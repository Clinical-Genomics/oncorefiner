#!/usr/bin/env python3
"""
Reformat ESVEE structural-variants from PURPLE VCF.

Supported SVTYPE values
-----------------------
The input VCF may contain only the following SVTYPE values:

* SGL
* BND
* DEL
* DUP
* INV
* INS

DEL, DUP, INV and INS are reformatable SVTYPE values. SGL variants are
reformatted as BND variants, while existing BND variants are retained unchanged.
Any missing or unsupported SVTYPE causes the program to terminate with an error.

Reformatting rules
------------------
* SVTYPE=SGL is changed to SVTYPE=BND. Its original ALT is retained.
* Existing SVTYPE=BND variants are retained unchanged.
* DEL, DUP, INV and INS variants are reformatted independently:
  - no variant lookup or mate pairing is performed;
  - the variant's own breakend ALT is parsed to obtain the remote coordinate;
  - ALT becomes <DEL>, <DUP>, <INV> or <INS>;
  - END is set to the remote coordinate encoded in that same ALT;
  - sequence retained from the original ALT column is stored in ALTCOLUMNSEQ.
* DEL, DUP, INV and INS variants are expected to satisfy all reformatting
  requirements. Unexpected variant structure causes the program to terminate
  with an error.

Every valid input variant produces exactly one output variant. Paired variants
therefore remain as two separate output variants even when their ID/MATEID
fields link them.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import click
import pysam


SVTYPE_CONFIG = {
    "SGL": {"action": "relabel"},
    "BND": {"action": "retain"},
    "DEL": {"action": "reformat"},
    "DUP": {"action": "reformat"},
    "INV": {"action": "reformat"},
    "INS": {"action": "reformat"},
}

BREAKEND_ALT_RE = re.compile(
    r"^(?P<left>[^\[\]]*)"
    r"(?P<bracket>[\[\]])"
    r"(?P<chrom>[^:\[\]]+):(?P<pos>\d+)"
    r"(?P=bracket)"
    r"(?P<right>[^\[\]]*)$"
)


@dataclass(frozen=True)
class ParsedBreakendAlt:
    """Parsed components of one breakend-style ALT allele.

    A breakend ALT combines sequence from the local VCF variant with a bracketed
    reference to the genomic position joined to that variant. Here, "remote"
    refers to that bracketed partner locus, while left_sequence and
    right_sequence are the ordinary bases surrounding it in the ALT column.
    """

    remote_chrom: str
    remote_pos: int
    left_sequence: str
    right_sequence: str

    @property
    def local_sequence(self) -> str:
        """Return all ordinary sequence surrounding the bracketed remote locus.

        This sequence still includes the local REF allele used to anchor the
        variant, so it is not yet the final ALTCOLUMNSEQ value.
        """
        return self.left_sequence + self.right_sequence

    @property
    def remote_first(self) -> bool:
        """Return True when the bracketed partner locus appears first in ALT.

        For example, ]chr1:200]CA has the remote locus first, whereas
        AC[chr1:200[ has local sequence first. This orientation determines which
        end of the local sequence contains the REF allele that must be removed.
        """
        return not self.left_sequence


@dataclass
class Report:
    """Collect summary statistics and optional per-variant-record reformatting details."""

    input_variant_records: int = 0
    output_variant_records: int = 0
    sgl_reformatted_to_bnd: int = 0
    reformatted_variant_records: int = 0
    reformatted_variant_records_with_altcolumnseq: int = 0
    unchanged_variant_records: int = 0
    input_svtypes: Counter[str] = field(default_factory=Counter)
    output_svtypes: Counter[str] = field(default_factory=Counter)
    reformatted_by_svtype: Counter[str] = field(default_factory=Counter)
    details: list[tuple[str, str, str]] = field(default_factory=list)

    def action(self, variant_record_label: str, message: str) -> None:
        """Store a per-variant-record action for an optional verbose report."""
        self.details.append(("action", variant_record_label, message))


def variant_record_label(variant_record: pysam.VariantRecord) -> str:
    """Return a concise label identifying one VCF variant record."""
    if variant_record.id is not None:
        return f"ID={variant_record.id}"
    return f"{variant_record.contig}:{variant_record.pos}"


def get_svtype(variant_record: pysam.VariantRecord) -> str:
    """Return and validate one variant record's required INFO/SVTYPE value."""
    svtype = variant_record.info.get("SVTYPE")

    if svtype is None:
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: missing required INFO/SVTYPE"
        )

    if not isinstance(svtype, str):
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: expected INFO/SVTYPE to contain one string, "
            f"but found {svtype!r}"
        )

    if svtype not in SVTYPE_CONFIG:
        supported = ", ".join(sorted(SVTYPE_CONFIG))
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: unsupported SVTYPE {svtype!r}. "
            f"Supported values are: {supported}"
        )

    return svtype


def parse_breakend_alt(alt: str) -> ParsedBreakendAlt | None:
    """Parse one breakend ALT allele, returning None when the syntax is invalid.

    Breakend notation uses matching square brackets around a chromosome and
    position, for example AC[chr1:200[ or ]chr1:200]CA. The chromosome and
    position describe the genomic locus connected to the current VCF variant record.
    """
    match = BREAKEND_ALT_RE.fullmatch(alt)
    if match is None:
        return None

    return ParsedBreakendAlt(
        remote_chrom=match.group("chrom"),
        remote_pos=int(match.group("pos")),
        left_sequence=match.group("left"),
        right_sequence=match.group("right"),
    )


def get_single_alt(variant_record: pysam.VariantRecord) -> str | None:
    """Return the single ALT allele, or None when exactly one ALT is not present."""
    if variant_record.alts is None or len(variant_record.alts) != 1:
        return None
    return variant_record.alts[0]


def extract_alt_column_sequence(
    variant_record: pysam.VariantRecord,
    parsed: ParsedBreakendAlt,
) -> str | None:
    """Extract ALT-column sequence after removing the local REF allele.

    In VCF, the ALT allele normally contains the REF allele from the current
    variant record so that the event is anchored to a concrete reference position.
    That anchoring sequence is representation rather than additional breakpoint
    sequence, so it is removed before storing ALTCOLUMNSEQ.
    """
    sequence = parsed.local_sequence
    ref = variant_record.ref

    # The local REF allele can occur on either side of the bracketed partner
    # locus. For example, with REF=A, AC[chr1:200[ starts with the REF allele,
    # whereas ]chr1:200]CA ends with it. Remove that REF sequence and retain only
    # the remaining bases for ALTCOLUMNSEQ.
    if parsed.remote_first:
        if not sequence.endswith(ref):
            return None
        return sequence[:-len(ref)] if ref else sequence

    if not sequence.startswith(ref):
        return None
    return sequence[len(ref):]


def reformat_sv_variant_record(variant_record: pysam.VariantRecord, svtype: str) -> pysam.VariantRecord:
    """Reformat one DEL/DUP/INV/INS variant record without consulting its mate."""
    if SVTYPE_CONFIG[svtype]["action"] != "reformat":
        raise ValueError(
            f"Internal error: reformat_sv_variant_record received SVTYPE {svtype!r} "
            "without action='reformat'"
        )

    alt_string = get_single_alt(variant_record)
    if alt_string is None:
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: SVTYPE={svtype} must contain exactly one "
            "ALT allele"
        )

    parsed = parse_breakend_alt(alt_string)
    if parsed is None:
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: SVTYPE={svtype} has an ALT allele that "
            f"cannot be parsed as a breakend: {alt_string!r}"
        )

    # In this ESVEE representation, DEL/DUP/INV/INS are interval events whose
    # two breakpoints lie on the same chromosome. The bracketed locus in ALT is
    # therefore expected to point back to the current chromosome. A different
    # chromosome would describe a translocation-like breakend instead.
    if parsed.remote_chrom != variant_record.contig:
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: SVTYPE={svtype} is expected to be "
            f"intrachromosomal, but ALT points to "
            f"{parsed.remote_chrom}:{parsed.remote_pos}"
        )

    altcolumn_sequence = extract_alt_column_sequence(variant_record, parsed)
    if altcolumn_sequence is None:
        raise click.ClickException(
            f"{variant_record_label(variant_record)}: SVTYPE={svtype} has ALT {alt_string!r}, "
            f"but the REF anchor {variant_record.ref!r} could not be identified"
        )

    reformatted = variant_record.copy()
    reformatted.alts = (f"<{svtype}>",)

    # pysam represents the VCF END coordinate as VariantRecord.stop rather than
    # as an ordinary INFO value. Assigning stop is therefore the pysam-supported
    # way of updating INFO/END in the variant record written to the output VCF.
    reformatted.stop = parsed.remote_pos

    if altcolumn_sequence:
        reformatted.info["ALTCOLUMNSEQ"] = altcolumn_sequence
    elif "ALTCOLUMNSEQ" in reformatted.info:
        # Remove a pre-existing value when this ALT contains no additional
        # sequence after the REF anchor has been stripped.
        del reformatted.info["ALTCOLUMNSEQ"]

    return reformatted


def reformat_variant_record(
    variant_record: pysam.VariantRecord,
    report: Report,
) -> pysam.VariantRecord:
    """Reformat exactly one variant record independently of every other variant record."""
    svtype = get_svtype(variant_record)
    action = SVTYPE_CONFIG[svtype]["action"]

    if action == "relabel":
        # An SGL variant record already carries a breakend-style ALT, but ESVEE labels it
        # as a single breakend. Downstream we want the generic VCF BND type, so
        # only SVTYPE changes; the ALT itself is deliberately left untouched.
        reformatted = variant_record.copy()
        reformatted.info["SVTYPE"] = "BND"
        report.sgl_reformatted_to_bnd += 1
        report.action(variant_record_label(variant_record), "reformatted SVTYPE=SGL to SVTYPE=BND")
        return reformatted

    if action == "retain":
        # Existing BND variant records are already expressed in standard breakend
        # notation, so no ALT, END, or INFO-field reformatting is required.
        report.unchanged_variant_records += 1
        report.action(variant_record_label(variant_record), "retained existing SVTYPE=BND unchanged")
        return variant_record

    if action != "reformat":
        # This guards against an invalid action being introduced into
        # SVTYPE_CONFIG in a future code change.
        raise ValueError(
            f"Internal error: unsupported action {action!r} configured for "
            f"SVTYPE {svtype!r}"
        )

    reformatted = reformat_sv_variant_record(variant_record, svtype)
    report.reformatted_variant_records += 1
    report.reformatted_by_svtype[svtype] += 1

    if reformatted.info.get("ALTCOLUMNSEQ"):
        report.reformatted_variant_records_with_altcolumnseq += 1

    report.action(
        variant_record_label(variant_record),
        f"reformatted independently to <{svtype}> with END={reformatted.stop}",
    )
    return reformatted


def process_variant_records(
    variant_records: list[pysam.VariantRecord],
    report: Report,
) -> list[pysam.VariantRecord]:
    """Reformat every variant record independently while preserving input order."""
    report.input_variant_records = len(variant_records)
    report.input_svtypes.update(get_svtype(variant_record) for variant_record in variant_records)

    # Reformat every variant record in memory before opening the output file. Because
    # this script treats malformed SV variant records as fatal errors, doing the work
    # first prevents an apparently valid but incomplete VCF from being left
    # behind if a later variant record fails validation.
    reformatted_variant_records = [reformat_variant_record(variant_record, report) for variant_record in variant_records]

    report.output_variant_records = len(reformatted_variant_records)
    report.output_svtypes.update(
        get_svtype(variant_record) for variant_record in reformatted_variant_records
    )
    return reformatted_variant_records


def add_required_headers(
    header: pysam.VariantHeader,
) -> pysam.VariantHeader:
    """Return a header copy containing INFO definitions used by this script."""
    output_header = header.copy()

    if "END" not in output_header.info:
        output_header.info.add(
            "END",
            number=1,
            type="Integer",
            description="Remote coordinate parsed from the original breakend ALT",
        )

    if "ALTCOLUMNSEQ" not in output_header.info:
        output_header.info.add(
            "ALTCOLUMNSEQ",
            number=1,
            type="String",
            description="Sequence retained from this ESVEE breakend ALT allele column",
        )

    return output_header


def read_variant_records(
    input_path: Path | None,
) -> tuple[pysam.VariantHeader, list[pysam.VariantRecord]]:
    """Read the input VCF and copy all variant records onto the augmented output header."""
    filename = "-" if input_path is None else str(input_path)

    try:
        with pysam.VariantFile(filename, "r") as reader:
            output_header = add_required_headers(reader.header)
            variant_records: list[pysam.VariantRecord] = []

            for variant_record in reader:
                # Each VCF variant record (`pysam.VariantRecord`) is tied to the header it was read
                # with. Because output_header contains newly added INFO
                # definitions such as ALTCOLUMNSEQ, copied variant records must first be
                # translated to that header before those fields can be assigned.
                copied_variant_record = variant_record.copy()
                copied_variant_record.translate(output_header)
                variant_records.append(copied_variant_record)

            return output_header, variant_records
    except (OSError, ValueError) as error:
        source = "stdin" if input_path is None else str(input_path)
        raise click.ClickException(
            f"Failed to read VCF from {source}: {error}"
        ) from error


def output_mode(path: Path | None) -> str:
    """Return the pysam write mode for stdout, plain VCF, or compressed VCF."""
    if path is None:
        return "w"
    return "wz" if path.suffix == ".gz" else "w"


def write_variant_records(
    output_path: Path | None,
    header: pysam.VariantHeader,
    variant_records: Iterable[pysam.VariantRecord],
) -> None:
    """Write reformatted variant records to stdout, VCF, or BGZF-compressed VCF."""
    filename = "-" if output_path is None else str(output_path)

    try:
        with pysam.VariantFile(
            filename,
            output_mode(output_path),
            header=header,
        ) as writer:
            for variant_record in variant_records:
                writer.write(variant_record)
    except (OSError, ValueError) as error:
        destination = "stdout" if output_path is None else str(output_path)
        raise click.ClickException(
            f"Failed to write VCF to {destination}: {error}"
        ) from error


def write_report(
    path: Path,
    report: Report,
    include_details: bool,
) -> None:
    """Write reformatting summary statistics and optional per-variant-record details."""
    rows: list[tuple[str, str, str]] = [
        ("summary", "input_variant_records", str(report.input_variant_records)),
        ("summary", "output_variant_records", str(report.output_variant_records)),
        (
            "summary",
            "variant_records_removed",
            str(report.input_variant_records - report.output_variant_records),
        ),
        (
            "summary",
            "sgl_reformatted_to_bnd",
            str(report.sgl_reformatted_to_bnd),
        ),
        ("summary", "reformatted_variant_records", str(report.reformatted_variant_records)),
        (
            "summary",
            "reformatted_variant_records_with_altcolumnseq",
            str(report.reformatted_variant_records_with_altcolumnseq),
        ),
        ("summary", "unchanged_variant_records", str(report.unchanged_variant_records)),
    ]

    for svtype, count in sorted(report.input_svtypes.items()):
        rows.append(("input_svtype", svtype, str(count)))

    for svtype, count in sorted(report.output_svtypes.items()):
        rows.append(("output_svtype", svtype, str(count)))

    for svtype, count in sorted(report.reformatted_by_svtype.items()):
        rows.append(("reformatted_by_svtype", svtype, str(count)))

    try:
        with path.open("w") as handle:
            handle.write("section\tname\tvalue\n")
            for section, name, value in rows:
                handle.write(f"{section}\t{name}\t{value}\n")

            if include_details:
                handle.write("\nkind\tvariant_record\tdescription\n")
                for kind, variant_record, description in report.details:
                    handle.write(f"{kind}\t{variant_record}\t{description}\n")
    except OSError as error:
        raise click.ClickException(
            f"Failed to write report to {path}: {error}"
        ) from error


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument(
    "input_vcf",
    required=False,
    type=click.Path(
        path_type=Path,
        exists=True,
        dir_okay=False,
        readable=True,
    ),
)
@click.option(
    "-o",
    "--output",
    "output_vcf",
    type=click.Path(path_type=Path, dir_okay=False),
    help="Output .vcf or BGZF-compressed .vcf.gz. Defaults to stdout.",
)
@click.option(
    "--report",
    "report_path",
    type=click.Path(path_type=Path, dir_okay=False),
    help="Write a TSV reformatting summary to this path.",
)
@click.option(
    "--verbose-report",
    is_flag=True,
    help="Append per-variant-record reformatting actions to the report. Requires --report.",
)
def main(
    input_vcf: Path | None,
    output_vcf: Path | None,
    report_path: Path | None,
    verbose_report: bool,
) -> None:
    """Reformat ESVEE SV variant records independently without merging mate pairs.

    INPUT_VCF may be plain text or BGZF/gzip-compressed. When omitted, input is
    read from stdin. Missing, unsupported, or structurally unexpected SV variant records
    terminate the program with an error.
    """
    if verbose_report and report_path is None:
        raise click.UsageError("--verbose-report requires --report")

    header, variant_records = read_variant_records(input_vcf)

    report = Report()
    reformatted_variant_records = process_variant_records(variant_records, report)
    write_variant_records(output_vcf, header, reformatted_variant_records)

    if report_path is not None:
        write_report(report_path, report, verbose_report)


if __name__ == "__main__":
    main()
