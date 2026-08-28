#!/usr/bin/env python3
"""
Reformat ESVEE structural-variant records independently.

Supported SVTYPE values
-----------------------
The input VCF may contain only the following SVTYPE values:

* SGL
* BND
* DEL
* DUP
* INV
* INS

DEL, DUP, INV and INS are reformatable SVTYPE values. SGL records are
reformatted as BND records, while existing BND records are retained unchanged.
Any missing or unsupported SVTYPE causes the program to terminate with an error.

Reformatting rules
------------------
* SVTYPE=SGL is changed to SVTYPE=BND. Its original ALT is retained.
* Existing SVTYPE=BND records are retained unchanged.
* DEL, DUP, INV and INS records are reformatted independently:
  - no record lookup or mate pairing is performed;
  - the record's own breakend ALT is parsed to obtain the remote coordinate;
  - ALT becomes <DEL>, <DUP>, <INV> or <INS>;
  - END is set to the remote coordinate encoded in that same ALT;
  - sequence retained from the original ALT column is stored in ALTCOLUMNSEQ.
* DEL, DUP, INV and INS records are expected to satisfy all reformatting
  requirements. Unexpected record structure causes the program to terminate
  with an error.

Every valid input record produces exactly one output record. Paired records
therefore remain as two separate output records even when their ID/MATEID
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
    """Parsed components of one breakend-style ALT allele."""

    remote_chrom: str
    remote_pos: int
    left_sequence: str
    right_sequence: str

    @property
    def local_sequence(self) -> str:
        """Return all non-breakend sequence from the ALT allele."""
        return self.left_sequence + self.right_sequence

    @property
    def remote_first(self) -> bool:
        """Return True when the bracketed remote locus precedes local sequence."""
        return not self.left_sequence


@dataclass
class Report:
    """Collect summary statistics and optional per-record reformatting details."""

    input_records: int = 0
    output_records: int = 0
    sgl_reformatted_to_bnd: int = 0
    reformatted_records: int = 0
    reformatted_with_altcolumnseq: int = 0
    unchanged_records: int = 0
    input_svtypes: Counter[str] = field(default_factory=Counter)
    output_svtypes: Counter[str] = field(default_factory=Counter)
    reformatted_by_svtype: Counter[str] = field(default_factory=Counter)
    details: list[tuple[str, str, str]] = field(default_factory=list)

    def action(self, record_label: str, message: str) -> None:
        """Store a per-record action for an optional verbose report."""
        self.details.append(("action", record_label, message))


def record_label(record: pysam.VariantRecord) -> str:
    """Return a concise record label for errors and report messages."""
    if record.id is not None:
        return f"ID={record.id}"
    return f"{record.contig}:{record.pos}"


def get_svtype(record: pysam.VariantRecord) -> str:
    """Return and validate the record's required single-valued INFO/SVTYPE."""
    svtype = record.info.get("SVTYPE")

    if svtype is None:
        raise click.ClickException(
            f"{record_label(record)}: missing required INFO/SVTYPE"
        )

    if not isinstance(svtype, str):
        raise click.ClickException(
            f"{record_label(record)}: expected INFO/SVTYPE to contain one string, "
            f"but found {svtype!r}"
        )

    if svtype not in SVTYPE_CONFIG:
        supported = ", ".join(sorted(SVTYPE_CONFIG))
        raise click.ClickException(
            f"{record_label(record)}: unsupported SVTYPE {svtype!r}. "
            f"Supported values are: {supported}"
        )

    return svtype


def parse_breakend_alt(alt: str) -> ParsedBreakendAlt | None:
    """Parse one breakend ALT allele, returning None when the syntax is invalid."""
    match = BREAKEND_ALT_RE.fullmatch(alt)
    if match is None:
        return None

    return ParsedBreakendAlt(
        remote_chrom=match.group("chrom"),
        remote_pos=int(match.group("pos")),
        left_sequence=match.group("left"),
        right_sequence=match.group("right"),
    )


def get_single_alt(record: pysam.VariantRecord) -> str | None:
    """Return the single ALT allele, or None when exactly one ALT is not present."""
    if record.alts is None or len(record.alts) != 1:
        return None
    return record.alts[0]


def extract_alt_column_sequence(
    record: pysam.VariantRecord,
    parsed: ParsedBreakendAlt,
) -> str | None:
    """Extract sequence added around the breakend after removing the REF anchor."""
    sequence = parsed.local_sequence
    ref = record.ref

    # In breakend notation the local REF anchor may occur either before or after
    # the bracketed remote locus. Which side it occupies determines where it must
    # be removed from the retained ALT-column sequence.
    if parsed.remote_first:
        if not sequence.endswith(ref):
            return None
        return sequence[:-len(ref)] if ref else sequence

    if not sequence.startswith(ref):
        return None
    return sequence[len(ref):]


def reformat_sv_record(record: pysam.VariantRecord, svtype: str) -> pysam.VariantRecord:
    """Reformat one DEL/DUP/INV/INS record without consulting its mate."""
    if SVTYPE_CONFIG[svtype]["action"] != "reformat":
        raise ValueError(
            f"Internal error: reformat_sv_record received SVTYPE {svtype!r} "
            "without action='reformat'"
        )

    alt_string = get_single_alt(record)
    if alt_string is None:
        raise click.ClickException(
            f"{record_label(record)}: SVTYPE={svtype} must contain exactly one "
            "ALT allele"
        )

    parsed = parse_breakend_alt(alt_string)
    if parsed is None:
        raise click.ClickException(
            f"{record_label(record)}: SVTYPE={svtype} has an ALT allele that "
            f"cannot be parsed as a breakend: {alt_string!r}"
        )

    # DEL/DUP/INV/INS describe an interval on one chromosome in this input
    # format. An interchromosomal remote coordinate therefore indicates input
    # that does not satisfy the assumptions of this reformatter.
    if parsed.remote_chrom != record.contig:
        raise click.ClickException(
            f"{record_label(record)}: SVTYPE={svtype} is expected to be "
            f"intrachromosomal, but ALT points to "
            f"{parsed.remote_chrom}:{parsed.remote_pos}"
        )

    altcolumn_sequence = extract_alt_column_sequence(record, parsed)
    if altcolumn_sequence is None:
        raise click.ClickException(
            f"{record_label(record)}: SVTYPE={svtype} has ALT {alt_string!r}, "
            f"but the REF anchor {record.ref!r} could not be identified"
        )

    reformatted = record.copy()
    reformatted.alts = (f"<{svtype}>",)

    # pysam/HTSlib exposes INFO/END through VariantRecord.stop, so assigning
    # stop updates END in the serialized VCF record.
    reformatted.stop = parsed.remote_pos

    if altcolumn_sequence:
        reformatted.info["ALTCOLUMNSEQ"] = altcolumn_sequence
    elif "ALTCOLUMNSEQ" in reformatted.info:
        # Remove a pre-existing value when this ALT contains no additional
        # sequence after the REF anchor has been stripped.
        del reformatted.info["ALTCOLUMNSEQ"]

    return reformatted


def reformat_record(
    record: pysam.VariantRecord,
    report: Report,
) -> pysam.VariantRecord:
    """Reformat exactly one record independently of every other record."""
    svtype = get_svtype(record)
    action = SVTYPE_CONFIG[svtype]["action"]

    if action == "relabel":
        # SGL already uses breakend-style ALT notation. Only its SVTYPE label is
        # changed so the resulting record is represented as BND.
        reformatted = record.copy()
        reformatted.info["SVTYPE"] = "BND"
        report.sgl_reformatted_to_bnd += 1
        report.action(record_label(record), "reformatted SVTYPE=SGL to SVTYPE=BND")
        return reformatted

    if action == "retain":
        # Existing BND records already have the desired representation.
        report.unchanged_records += 1
        report.action(record_label(record), "retained existing SVTYPE=BND unchanged")
        return record

    if action != "reformat":
        # This guards against an invalid action being introduced into
        # SVTYPE_CONFIG in a future code change.
        raise ValueError(
            f"Internal error: unsupported action {action!r} configured for "
            f"SVTYPE {svtype!r}"
        )

    reformatted = reformat_sv_record(record, svtype)
    report.reformatted_records += 1
    report.reformatted_by_svtype[svtype] += 1

    if reformatted.info.get("ALTCOLUMNSEQ"):
        report.reformatted_with_altcolumnseq += 1

    report.action(
        record_label(record),
        f"reformatted independently to <{svtype}> with END={reformatted.stop}",
    )
    return reformatted


def reformat_records(
    records: list[pysam.VariantRecord],
    report: Report,
) -> list[pysam.VariantRecord]:
    """Reformat every input record independently while preserving input order."""
    report.input_records = len(records)
    report.input_svtypes.update(get_svtype(record) for record in records)

    # Build the complete output list before opening the output VCF. This means
    # any unsupported or malformed input record aborts the run before a partial
    # output file can be written.
    reformatted_records = [reformat_record(record, report) for record in records]

    report.output_records = len(reformatted_records)
    report.output_svtypes.update(
        get_svtype(record) for record in reformatted_records
    )
    return reformatted_records


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


def read_records(
    input_path: Path | None,
) -> tuple[pysam.VariantHeader, list[pysam.VariantRecord]]:
    """Read the input VCF and copy all records onto the augmented output header."""
    filename = "-" if input_path is None else str(input_path)

    try:
        with pysam.VariantFile(filename, "r") as reader:
            output_header = add_required_headers(reader.header)
            records: list[pysam.VariantRecord] = []

            for record in reader:
                # Records must be translated onto the copied header before new
                # INFO fields such as ALTCOLUMNSEQ can safely be assigned.
                copied = record.copy()
                copied.translate(output_header)
                records.append(copied)

            return output_header, records
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


def write_records(
    output_path: Path | None,
    header: pysam.VariantHeader,
    records: Iterable[pysam.VariantRecord],
) -> None:
    """Write reformatted records to stdout, VCF, or BGZF-compressed VCF."""
    filename = "-" if output_path is None else str(output_path)

    try:
        with pysam.VariantFile(
            filename,
            output_mode(output_path),
            header=header,
        ) as writer:
            for record in records:
                writer.write(record)
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
    """Write reformatting summary statistics and optional per-record details."""
    rows: list[tuple[str, str, str]] = [
        ("summary", "input_records", str(report.input_records)),
        ("summary", "output_records", str(report.output_records)),
        (
            "summary",
            "records_removed",
            str(report.input_records - report.output_records),
        ),
        (
            "summary",
            "sgl_reformatted_to_bnd",
            str(report.sgl_reformatted_to_bnd),
        ),
        ("summary", "reformatted_records", str(report.reformatted_records)),
        (
            "summary",
            "reformatted_records_with_altcolumnseq",
            str(report.reformatted_with_altcolumnseq),
        ),
        ("summary", "unchanged_records", str(report.unchanged_records)),
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
                handle.write("\nkind\trecord\tdescription\n")
                for kind, record, description in report.details:
                    handle.write(f"{kind}\t{record}\t{description}\n")
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
    help="Append per-record reformatting actions to the report. Requires --report.",
)
def main(
    input_vcf: Path | None,
    output_vcf: Path | None,
    report_path: Path | None,
    verbose_report: bool,
) -> None:
    """Reformat ESVEE SV records independently without merging mate pairs.

    INPUT_VCF may be plain text or BGZF/gzip-compressed. When omitted, input is
    read from stdin. Missing, unsupported, or structurally unexpected SV records
    terminate the program with an error.
    """
    if verbose_report and report_path is None:
        raise click.UsageError("--verbose-report requires --report")

    header, records = read_records(input_vcf)

    report = Report()
    reformatted_records = reformat_records(records, report)
    write_records(output_vcf, header, reformatted_records)

    if report_path is not None:
        write_report(report_path, report, verbose_report)


if __name__ == "__main__":
    main()
