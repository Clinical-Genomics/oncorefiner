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
  - SVID and the template record's other annotations are retained;
  - non-reference junction sequence is extracted from the canonical
    lower-coordinate record and stored in INSSEQ.
* Structurally unsafe pairs are retained unchanged and reported.

The program reads the whole VCF before transforming it. This permits reliable
mate pairing and complete conversion statistics.
"""

from __future__ import annotations

import gzip
import re
import sys
from collections import Counter
from contextlib import nullcontext
from dataclasses import dataclass, field
from pathlib import Path
from typing import ContextManager, Iterable, TextIO

import click


COLLAPSIBLE_TYPES = {"DEL", "DUP", "INV", "INS"}

BREAKEND_ALT_RE = re.compile(
    r"^(?P<left>[^\[\]]*)"
    r"(?P<bracket>[\[\]])"
    r"(?P<chrom>[^:\[\]]+):(?P<pos>\d+)"
    r"(?P=bracket)"
    r"(?P<right>[^\[\]]*)$"
)


@dataclass
class VcfRecord:
    fields: list[str]
    input_index: int

    @classmethod
    def from_line(cls, line: str, input_index: int) -> "VcfRecord":
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise click.ClickException(
                f"Malformed VCF record {input_index + 1}: expected at least "
                f"8 columns, found {len(fields)}"
            )
        try:
            int(fields[1])
        except ValueError as error:
            raise click.ClickException(
                f"Malformed POS in record {input_index + 1}: {fields[1]!r}"
            ) from error
        return cls(fields=fields, input_index=input_index)

    def copy(self) -> "VcfRecord":
        return VcfRecord(self.fields.copy(), self.input_index)

    def to_line(self) -> str:
        return "\t".join(self.fields) + "\n"

    @property
    def chrom(self) -> str:
        return self.fields[0]

    @property
    def pos(self) -> int:
        return int(self.fields[1])

    @property
    def record_id(self) -> str:
        return self.fields[2]

    @property
    def ref(self) -> str:
        return self.fields[3]

    @property
    def alt(self) -> str:
        return self.fields[4]

    @alt.setter
    def alt(self, value: str) -> None:
        self.fields[4] = value

    @property
    def info(self) -> list[str]:
        value = self.fields[7]
        return [] if value in {"", "."} else value.split(";")

    @info.setter
    def info(self, entries: Iterable[str]) -> None:
        items = list(entries)
        self.fields[7] = ";".join(items) if items else "."

    def info_value(self, key: str) -> str | None:
        prefix = f"{key}="
        for entry in self.info:
            if entry.startswith(prefix):
                return entry[len(prefix):]
        return None

    def set_info(self, key: str, value: str) -> None:
        prefix = f"{key}="
        entries = [
            entry
            for entry in self.info
            if entry != key and not entry.startswith(prefix)
        ]
        entries.append(f"{key}={value}")
        self.info = entries

    def remove_info(self, key: str) -> None:
        prefix = f"{key}="
        self.info = [
            entry
            for entry in self.info
            if entry != key and not entry.startswith(prefix)
        ]


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
        """True when the bracketed remote locus appears before local sequence."""
        return not self.left_sequence


@dataclass
class Report:
    input_records: int = 0
    output_records: int = 0
    sgl_to_bnd: int = 0
    merged_pairs: int = 0
    symbolic_with_insseq: int = 0
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


def normalize_svtype(value: str | None) -> str:
    return value if value else "MISSING"


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


def extract_non_anchor_sequence(record: VcfRecord) -> str | None:
    """
    Extract sequence beyond the local REF anchor from a breakend ALT.

    VCF breakend placement determines which side contains the local anchor:

    * local sequence before the bracket, e.g. A[chr1:100[
      -> REF anchor is at the beginning;
    * bracketed remote locus before local sequence, e.g. ]chr1:100]ACGT
      -> REF anchor is at the end.

    Using bracket placement avoids ambiguity when the local sequence happens
    to begin and end with the same base.
    """
    parsed = parse_breakend_alt(record.alt)
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


def choose_insert_sequence(
    canonical: VcfRecord,
) -> tuple[str | None, str | None]:
    """
    Extract INSSEQ from the canonical record only.

    The mate's junction sequence is deliberately not compared. ESVEE may
    represent or assemble the two breakend alleles differently, and this
    converter is intended to preserve a deterministic sequence annotation,
    not validate assembly agreement between mates.
    """
    sequence = extract_non_anchor_sequence(canonical)
    if sequence is None:
        return None, "anchor_not_identified"
    return sequence, None


def validate_pair(first: VcfRecord, second: VcfRecord) -> tuple[bool, str, str]:
    svtype1 = first.info_value("SVTYPE")
    svtype2 = second.info_value("SVTYPE")

    if svtype1 != svtype2:
        return False, "different_svtype", (
            f"different SVTYPE values ({svtype1!r} and {svtype2!r})"
        )

    if svtype1 not in COLLAPSIBLE_TYPES:
        return False, "not_collapsible", f"SVTYPE {svtype1!r} is not collapsible"

    if first.chrom != second.chrom:
        return False, "interchromosomal_pair", "mates are on different chromosomes"

    if (
        first.info_value("MATEID") != second.record_id
        or second.info_value("MATEID") != first.record_id
    ):
        return False, "nonreciprocal_mateid", "MATEID links are not reciprocal"

    alt1 = parse_breakend_alt(first.alt)
    alt2 = parse_breakend_alt(second.alt)
    if alt1 is None or alt2 is None:
        return False, "alt_parse_failure", (
            "one or both ALT alleles are not parseable breakends"
        )

    if alt1.remote_chrom != second.chrom or alt1.remote_pos != second.pos:
        return False, "alt_mate_mismatch", "first ALT does not point to its mate"

    if alt2.remote_chrom != first.chrom or alt2.remote_pos != first.pos:
        return False, "alt_mate_mismatch", "second ALT does not point to its mate"

    return True, "", ""


def collapse_pair(
    first: VcfRecord,
    second: VcfRecord,
) -> tuple[VcfRecord | None, str | None, str | None]:
    valid, reason, message = validate_pair(first, second)
    if not valid:
        return None, reason, message

    canonical, partner = sorted(
        (first, second),
        key=lambda record: (record.chrom, record.pos, record.input_index),
    )

    insert_sequence, sequence_error = choose_insert_sequence(canonical)
    if sequence_error:
        return (
            None,
            sequence_error,
            "the local REF anchor could not be identified in the canonical ALT allele",
        )
    svtype = canonical.info_value("SVTYPE")
    assert svtype is not None

    output = canonical.copy()
    output.input_index = min(first.input_index, second.input_index)
    output.alt = f"<{svtype}>"
    output.remove_info("MATEID")
    output.set_info("END", str(partner.pos))

    if insert_sequence:
        output.set_info("INSSEQ", insert_sequence)
    else:
        output.remove_info("INSSEQ")

    return output, None, None


def build_record_index(
    records: list[VcfRecord],
    report: Report,
) -> dict[str, VcfRecord]:
    by_id: dict[str, VcfRecord] = {}
    duplicate_ids: set[str] = set()

    for record in records:
        record_id = record.record_id
        if record_id in {"", "."}:
            continue
        if record_id in by_id:
            duplicate_ids.add(record_id)
        else:
            by_id[record_id] = record

    for record_id in duplicate_ids:
        by_id.pop(record_id, None)
        report.warn(
            "duplicate_record_id",
            f"ID={record_id}",
            "duplicate record ID; records with this ID cannot be paired safely",
        )

    return by_id


def transform_records(records: list[VcfRecord], report: Report) -> list[VcfRecord]:
    """
    Two-pass transformation.

    Pass 1 builds the ID index. Pass 2 classifies and transforms each record or
    mate pair exactly once.
    """
    report.input_records = len(records)
    report.input_svtypes.update(
        normalize_svtype(record.info_value("SVTYPE")) for record in records
    )

    by_id = build_record_index(records, report)
    consumed: set[int] = set()
    output: list[VcfRecord] = []

    for record in records:
        if record.input_index in consumed:
            continue

        svtype = record.info_value("SVTYPE")

        if svtype == "SGL":
            converted = record.copy()
            converted.set_info("SVTYPE", "BND")
            output.append(converted)
            consumed.add(record.input_index)
            report.sgl_to_bnd += 1
            report.action(
                f"ID={record.record_id}",
                "converted SVTYPE=SGL to SVTYPE=BND",
            )
            continue

        if svtype == "BND" or svtype not in COLLAPSIBLE_TYPES:
            output.append(record)
            consumed.add(record.input_index)
            report.unchanged_records += 1
            continue

        mate_id = record.info_value("MATEID")
        mate = by_id.get(mate_id) if mate_id else None

        if mate is None:
            report.warn(
                "mate_not_available",
                f"ID={record.record_id}",
                f"mate {mate_id!r} was not available",
            )
            output.append(record)
            consumed.add(record.input_index)
            report.unchanged_records += 1
            continue

        if mate.input_index in consumed:
            report.warn(
                "mate_already_consumed",
                f"ID={record.record_id}",
                f"mate {mate_id!r} was already processed",
            )
            output.append(record)
            consumed.add(record.input_index)
            report.unchanged_records += 1
            continue

        collapsed, reason, message = collapse_pair(record, mate)
        pair_label = f"ID={record.record_id}/{mate.record_id}"

        if collapsed is None:
            assert reason is not None and message is not None
            report.warn(reason, pair_label, message)
            output.extend(sorted((record, mate), key=lambda item: item.input_index))
            report.unchanged_records += 2
        else:
            output.append(collapsed)
            merged_svtype = normalize_svtype(collapsed.info_value("SVTYPE"))
            report.merged_pairs += 1
            report.merged_by_svtype[merged_svtype] += 1
            if collapsed.info_value("INSSEQ"):
                report.symbolic_with_insseq += 1
            report.action(
                pair_label,
                f"merged as {merged_svtype} at "
                f"{collapsed.chrom}:{collapsed.pos}-{collapsed.info_value('END')}",
            )

        consumed.add(record.input_index)
        consumed.add(mate.input_index)

    output.sort(key=lambda record: record.input_index)
    report.output_records = len(output)
    report.output_svtypes.update(
        normalize_svtype(record.info_value("SVTYPE")) for record in output
    )
    return output


def open_text_input(path: Path | None) -> ContextManager[TextIO]:
    if path is None:
        return nullcontext(sys.stdin)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def open_text_output(path: Path | None) -> ContextManager[TextIO]:
    if path is None:
        return nullcontext(sys.stdout)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w")


def add_required_headers(headers: list[str]) -> list[str]:
    has_end = any(line.startswith("##INFO=<ID=END,") for line in headers)
    has_insseq = any(line.startswith("##INFO=<ID=INSSEQ,") for line in headers)

    additions: list[str] = []
    if not has_end:
        additions.append(
            '##INFO=<ID=END,Number=1,Type=Integer,'
            'Description="End coordinate of the symbolic structural variant">\n'
        )
    if not has_insseq:
        additions.append(
            '##INFO=<ID=INSSEQ,Number=1,Type=String,'
            'Description="Non-reference junction sequence extracted from the '
            'original ESVEE breakend ALT allele">\n'
        )

    result: list[str] = []
    inserted = False
    for line in headers:
        if line.startswith("#CHROM") and not inserted:
            result.extend(additions)
            inserted = True
        result.append(line)

    if not inserted:
        raise click.ClickException("VCF header does not contain a #CHROM line")

    return result


def write_report(path: Path, report: Report, include_details: bool) -> None:
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
            "symbolic_records_with_insseq",
            str(report.symbolic_with_insseq),
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

    with path.open("w") as handle:
        handle.write("section\tname\tvalue\n")
        for section, name, value in rows:
            handle.write(f"{section}\t{name}\t{value}\n")

        if include_details:
            handle.write("\nkind\trecords\tdescription\n")
            for kind, records, description in report.details:
                handle.write(f"{kind}\t{records}\t{description}\n")


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
    help="Output .vcf or .vcf.gz. Defaults to stdout.",
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

    INPUT_VCF may be plain text or gzip-compressed. When omitted, input is
    read from stdin.
    """
    if verbose_report and report_path is None:
        raise click.UsageError("--verbose-report requires --report")

    headers: list[str] = []
    records: list[VcfRecord] = []
    seen_column_header = False

    with open_text_input(input_vcf) as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                headers.append(line)
                if line.startswith("#CHROM"):
                    seen_column_header = True
                continue

            if not seen_column_header:
                raise click.ClickException(
                    "Encountered a VCF record before the #CHROM header"
                )

            if line.strip():
                records.append(VcfRecord.from_line(line, len(records)))

    report = Report()
    transformed = transform_records(records, report)
    output_headers = add_required_headers(headers)

    with open_text_output(output_vcf) as output_handle:
        output_handle.writelines(output_headers)
        for record in transformed:
            output_handle.write(record.to_line())

    if report_path is not None:
        write_report(report_path, report, verbose_report)


if __name__ == "__main__":
    main()
