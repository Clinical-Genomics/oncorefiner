#!/usr/bin/env python3
"""
Convert paired ESVEE breakend-style SV records into symbolic VCF records.

Conversion rules
----------------
* SVTYPE=SGL is changed to SVTYPE=BND. Its original ALT is retained.
* Existing SVTYPE=BND records are retained unchanged.
* Paired DEL, DUP, INV and INS records are collapsed into one record:
  - records are paired through reciprocal ID/MATEID links;
  - both ALT alleles must point to the mate coordinates;
  - the lower coordinate is used as POS;
  - ALT becomes <DEL>, <DUP>, <INV> or <INS>;
  - END is set to the partner coordinate;
  - MATEID is removed;
  - SVID and the canonical record's other annotations are retained;
  - non-anchor sequence from the canonical ESVEE breakend ALT is stored
    in JUNCTIONSEQ. This sequence may contain reference-derived assembly
    sequence and must not be interpreted as necessarily representing a
    novel insertion.
* Structurally unsafe pairs are retained unchanged and reported.

The program reads the whole VCF before transforming it. This permits reliable
mate pairing and complete conversion statistics.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import click
import pysam


COLLAPSIBLE_TYPES = {"DEL", "DUP", "INV", "INS"}

BREAKEND_ALT_RE = re.compile(
    r"^(?P<left>[^\[\]]*)"
    r"(?P<bracket>[\[\]])"
    r"(?P<chrom>[^:\[\]]+):(?P<pos>\d+)"
    r"(?P=bracket)"
    r"(?P<right>[^\[\]]*)$"
)


@dataclass(frozen=True)
class ParsedBreakendAlt:
    remote_chrom: str
    remote_pos: int
    left_sequence: str
    right_sequence: str

    @property
    def local_sequence(self) -> str:
        return self.left_sequence + self.right_sequence

    @property
    def remote_first(self) -> bool:
        """Return True when the bracketed remote locus precedes local sequence."""
        return not self.left_sequence


@dataclass
class Report:
    input_records: int = 0
    output_records: int = 0
    sgl_to_bnd: int = 0
    merged_pairs: int = 0
    symbolic_with_junctionseq: int = 0
    unchanged_records: int = 0
    input_svtypes: Counter[str] = field(default_factory=Counter)
    output_svtypes: Counter[str] = field(default_factory=Counter)
    merged_by_svtype: Counter[str] = field(default_factory=Counter)
    warnings_by_reason: Counter[str] = field(default_factory=Counter)
    details: list[tuple[str, str, str]] = field(default_factory=list)

    @property
    def records_removed_by_merging(self) -> int:
        return self.input_records - self.output_records

    @property
    def warning_count(self) -> int:
        return sum(self.warnings_by_reason.values())

    def warn(self, reason: str, record_label: str, message: str) -> None:
        self.warnings_by_reason[reason] += 1
        self.details.append(("warning", record_label, message))
        click.echo(f"Warning: {record_label}: {message}", err=True)

    def action(self, record_label: str, message: str) -> None:
        self.details.append(("action", record_label, message))


def normalize_svtype(value: object | None) -> str:
    return str(value) if value is not None else "MISSING"


def scalar_info(record: pysam.VariantRecord, key: str) -> str | None:
    """
    Return a single INFO value as a string.

    MATEID and SVTYPE are expected to be scalar in the ESVEE VCF. A one-item
    tuple is accepted for robustness, while multi-valued fields are rejected.
    """
    value = record.info.get(key)
    if value is None:
        return None

    if isinstance(value, tuple):
        if len(value) != 1:
            return None
        value = value[0]

    return str(value)


def record_label(record: pysam.VariantRecord) -> str:
    return f"ID={record.id}"


def parse_breakend_alt(alt: str) -> ParsedBreakendAlt | None:
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
    if record.alts is None or len(record.alts) != 1:
        return None
    return record.alts[0]


def extract_non_anchor_sequence(
    record: pysam.VariantRecord,
) -> str | None:
    """
    Extract the non-anchor sequence from a breakend ALT.

    VCF breakend placement determines which side contains the local REF anchor:

    * local sequence before the bracket, e.g. A[chr1:100[
      -> REF anchor is at the beginning;
    * bracketed remote locus before local sequence, e.g. ]chr1:100]ACGT
      -> REF anchor is at the end.

    The returned sequence is retained as JUNCTIONSEQ. It may include
    reference-derived ESVEE assembly sequence and is not asserted to be a
    novel insertion.
    """
    alt = get_single_alt(record)
    if alt is None:
        return None

    parsed = parse_breakend_alt(alt)
    if parsed is None:
        return None

    sequence = parsed.local_sequence
    ref = record.ref

    if parsed.remote_first:
        if not sequence.endswith(ref):
            return None
        return sequence[:-len(ref)] if ref else sequence

    if not sequence.startswith(ref):
        return None
    return sequence[len(ref):]


def choose_junction_sequence(
    canonical: pysam.VariantRecord,
) -> tuple[str | None, str | None]:
    """
    Extract JUNCTIONSEQ from the canonical record only.

    The mate sequence is deliberately not compared. ESVEE may represent or
    assemble the two breakend alleles differently; this converter preserves
    one deterministic sequence annotation rather than validating agreement
    between mates.
    """
    junction_sequence = extract_non_anchor_sequence(canonical)
    if junction_sequence is None:
        return None, "anchor_not_identified"
    return junction_sequence, None


def validate_pair(
    first: pysam.VariantRecord,
    second: pysam.VariantRecord,
) -> tuple[bool, str, str]:
    svtype1 = scalar_info(first, "SVTYPE")
    svtype2 = scalar_info(second, "SVTYPE")

    if svtype1 != svtype2:
        return False, "different_svtype", (
            f"different SVTYPE values ({svtype1!r} and {svtype2!r})"
        )

    if svtype1 not in COLLAPSIBLE_TYPES:
        return False, "not_collapsible", f"SVTYPE {svtype1!r} is not collapsible"

    if first.contig != second.contig:
        return False, "interchromosomal_pair", "mates are on different chromosomes"

    if (
        scalar_info(first, "MATEID") != second.id
        or scalar_info(second, "MATEID") != first.id
    ):
        return False, "nonreciprocal_mateid", "MATEID links are not reciprocal"

    alt1_string = get_single_alt(first)
    alt2_string = get_single_alt(second)
    if alt1_string is None or alt2_string is None:
        return False, "alt_count_invalid", (
            "one or both records do not have exactly one ALT allele"
        )

    alt1 = parse_breakend_alt(alt1_string)
    alt2 = parse_breakend_alt(alt2_string)
    if alt1 is None or alt2 is None:
        return False, "alt_parse_failure", (
            "one or both ALT alleles are not parseable breakends"
        )

    if alt1.remote_chrom != second.contig or alt1.remote_pos != second.pos:
        return False, "alt_mate_mismatch", "first ALT does not point to its mate"

    if alt2.remote_chrom != first.contig or alt2.remote_pos != first.pos:
        return False, "alt_mate_mismatch", "second ALT does not point to its mate"

    return True, "", ""


def collapse_pair(
    first: pysam.VariantRecord,
    second: pysam.VariantRecord,
) -> tuple[pysam.VariantRecord | None, str | None, str | None]:
    valid, reason, message = validate_pair(first, second)
    if not valid:
        return None, reason, message

    canonical, partner = sorted(
        (first, second),
        key=lambda record: (record.contig, record.pos),
    )

    junction_sequence, sequence_error = choose_junction_sequence(canonical)
    if sequence_error:
        return (
            None,
            sequence_error,
            "the local REF anchor could not be identified in the canonical ALT allele",
        )

    svtype = scalar_info(canonical, "SVTYPE")
    assert svtype is not None

    output = canonical.copy()
    output.alts = (f"<{svtype}>",)

    if "MATEID" in output.info:
        del output.info["MATEID"]

    # pysam/HTSlib exposes INFO/END through VariantRecord.stop rather than
    # VariantRecord.info. Assigning the one-based partner position here causes
    # the VCF writer to emit END=<partner.pos>.
    output.stop = partner.pos

    if junction_sequence:
        output.info["JUNCTIONSEQ"] = junction_sequence
    elif "JUNCTIONSEQ" in output.info:
        del output.info["JUNCTIONSEQ"]

    return output, None, None


def build_record_index(
    records: Iterable[pysam.VariantRecord],
    report: Report,
) -> dict[str, pysam.VariantRecord]:
    by_id: dict[str, pysam.VariantRecord] = {}
    duplicate_ids: set[str] = set()

    for record in records:
        record_id = record.id
        if record_id in {None, "", "."}:
            continue

        if record_id in by_id:
            duplicate_ids.add(record_id)
        else:
            by_id[record_id] = record

    for duplicate_id in duplicate_ids:
        by_id.pop(duplicate_id, None)
        report.warn(
            "duplicate_record_id",
            f"ID={duplicate_id}",
            "duplicate record ID; records with this ID cannot be paired safely",
        )

    return by_id


def transform_records(
    records: list[pysam.VariantRecord],
    report: Report,
) -> list[pysam.VariantRecord]:
    """
    Transform each record or reciprocal mate pair exactly once.

    Record object identity is used to track consumed input records because
    pysam VariantRecord objects do not carry the original line index.
    """
    report.input_records = len(records)
    report.input_svtypes.update(
        normalize_svtype(scalar_info(record, "SVTYPE")) for record in records
    )

    by_id = build_record_index(records, report)
    input_order = {id(record): index for index, record in enumerate(records)}
    consumed: set[int] = set()
    output_with_order: list[tuple[int, pysam.VariantRecord]] = []

    for record in records:
        record_key = id(record)
        if record_key in consumed:
            continue

        svtype = scalar_info(record, "SVTYPE")

        if svtype == "SGL":
            converted = record.copy()
            converted.info["SVTYPE"] = "BND"
            output_with_order.append((input_order[record_key], converted))
            consumed.add(record_key)
            report.sgl_to_bnd += 1
            report.action(
                record_label(record),
                "converted SVTYPE=SGL to SVTYPE=BND",
            )
            continue

        if svtype == "BND" or svtype not in COLLAPSIBLE_TYPES:
            output_with_order.append((input_order[record_key], record))
            consumed.add(record_key)
            report.unchanged_records += 1
            continue

        mate_id = scalar_info(record, "MATEID")
        mate = by_id.get(mate_id) if mate_id else None

        if mate is None:
            report.warn(
                "mate_not_available",
                record_label(record),
                f"mate {mate_id!r} was not available",
            )
            output_with_order.append((input_order[record_key], record))
            consumed.add(record_key)
            report.unchanged_records += 1
            continue

        mate_key = id(mate)
        if mate_key in consumed:
            report.warn(
                "mate_already_consumed",
                record_label(record),
                f"mate {mate_id!r} was already processed",
            )
            output_with_order.append((input_order[record_key], record))
            consumed.add(record_key)
            report.unchanged_records += 1
            continue

        collapsed, reason, message = collapse_pair(record, mate)
        pair_label = f"ID={record.id}/{mate.id}"
        pair_order = min(input_order[record_key], input_order[mate_key])

        if collapsed is None:
            assert reason is not None and message is not None
            report.warn(reason, pair_label, message)
            output_with_order.extend(
                sorted(
                    (
                        (input_order[record_key], record),
                        (input_order[mate_key], mate),
                    ),
                    key=lambda item: item[0],
                )
            )
            report.unchanged_records += 2
        else:
            output_with_order.append((pair_order, collapsed))
            merged_svtype = normalize_svtype(
                scalar_info(collapsed, "SVTYPE")
            )
            report.merged_pairs += 1
            report.merged_by_svtype[merged_svtype] += 1
            if scalar_info(collapsed, "JUNCTIONSEQ"):
                report.symbolic_with_junctionseq += 1
            report.action(
                pair_label,
                f"merged as {merged_svtype} at "
                f"{collapsed.contig}:{collapsed.pos}-{collapsed.stop}",
            )

        consumed.add(record_key)
        consumed.add(mate_key)

    output_with_order.sort(key=lambda item: item[0])
    output = [record for _, record in output_with_order]

    report.output_records = len(output)
    report.output_svtypes.update(
        normalize_svtype(scalar_info(record, "SVTYPE")) for record in output
    )
    return output


def add_required_headers(
    header: pysam.VariantHeader,
) -> pysam.VariantHeader:
    output_header = header.copy()

    if "END" not in output_header.info:
        output_header.info.add(
            "END",
            number=1,
            type="Integer",
            description="End coordinate of the symbolic structural variant",
        )

    if "JUNCTIONSEQ" not in output_header.info:
        output_header.info.add(
            "JUNCTIONSEQ",
            number=1,
            type="String",
            description=(
                "Non-anchor sequence retained from the canonical ESVEE "
                "breakend ALT allele; may include reference-derived assembly "
                "sequence and does not necessarily represent a novel insertion"
            ),
        )

    return output_header


def read_records(
    input_path: Path | None,
) -> tuple[pysam.VariantHeader, list[pysam.VariantRecord]]:
    filename = "-" if input_path is None else str(input_path)

    try:
        with pysam.VariantFile(filename, "r") as reader:
            output_header = add_required_headers(reader.header)
            records: list[pysam.VariantRecord] = []

            for record in reader:
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
    if path is None:
        return "w"
    return "wz" if path.suffix == ".gz" else "w"


def write_records(
    output_path: Path | None,
    header: pysam.VariantHeader,
    records: Iterable[pysam.VariantRecord],
) -> None:
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
    rows: list[tuple[str, str, str]] = [
        ("summary", "input_records", str(report.input_records)),
        ("summary", "output_records", str(report.output_records)),
        (
            "summary",
            "records_removed_by_merging",
            str(report.records_removed_by_merging),
        ),
        ("summary", "merged_pairs", str(report.merged_pairs)),
        ("summary", "sgl_converted_to_bnd", str(report.sgl_to_bnd)),
        (
            "summary",
            "symbolic_records_with_junctionseq",
            str(report.symbolic_with_junctionseq),
        ),
        ("summary", "unchanged_records", str(report.unchanged_records)),
        ("summary", "warnings", str(report.warning_count)),
    ]

    for svtype, count in sorted(report.input_svtypes.items()):
        rows.append(("input_svtype", svtype, str(count)))

    for svtype, count in sorted(report.output_svtypes.items()):
        rows.append(("output_svtype", svtype, str(count)))

    for svtype, count in sorted(report.merged_by_svtype.items()):
        rows.append(("merged_pairs_by_svtype", svtype, str(count)))

    for reason, count in sorted(report.warnings_by_reason.items()):
        rows.append(("warnings_by_reason", reason, str(count)))

    try:
        with path.open("w") as handle:
            handle.write("section\tname\tvalue\n")
            for section, name, value in rows:
                handle.write(f"{section}\t{name}\t{value}\n")

            if include_details:
                handle.write("\nkind\trecords\tdescription\n")
                for kind, records, description in report.details:
                    handle.write(f"{kind}\t{records}\t{description}\n")
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
    help="Write a TSV conversion summary to this path.",
)
@click.option(
    "--verbose-report",
    is_flag=True,
    help=(
        "Append per-pair actions and warnings to the report. "
        "Requires --report."
    ),
)
def main(
    input_vcf: Path | None,
    output_vcf: Path | None,
    report_path: Path | None,
    verbose_report: bool,
) -> None:
    """
    Collapse paired ESVEE DEL/DUP/INV/INS breakends into symbolic VCF records.

    INPUT_VCF may be plain text or BGZF/gzip-compressed. When omitted, input
    is read from stdin.
    """
    if verbose_report and report_path is None:
        raise click.UsageError("--verbose-report requires --report")

    header, records = read_records(input_vcf)

    report = Report()
    transformed = transform_records(records, report)
    write_records(output_vcf, header, transformed)

    if report_path is not None:
        write_report(report_path, report, verbose_report)


if __name__ == "__main__":
    main()
