#!/usr/bin/env python3
"""
Convert ESVEE breakend-style SV records independently.

Conversion rules
----------------
* SVTYPE=SGL is changed to SVTYPE=BND. Its original ALT is retained.
* Existing SVTYPE=BND records are retained unchanged.
* DEL, DUP, INV and INS records are converted independently:
  - no record lookup or mate pairing is performed;
  - the record's own breakend ALT is parsed to obtain the remote coordinate;
  - ALT becomes <DEL>, <DUP>, <INV> or <INS>;
  - END is set to the remote coordinate encoded in that same ALT;
  - non-anchor sequence from that record's own ALT is stored in JUNCTIONSEQ.
* Records that cannot be converted safely are retained unchanged and reported.

Every input record produces exactly one output record. Paired records therefore
remain as two separate output records even when their ID/MATEID fields link them.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Iterable

import click
import pysam


CONVERTIBLE_TYPES = {"DEL", "DUP", "INV", "INS"}

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
    symbolic_converted: int = 0
    symbolic_with_junctionseq: int = 0
    unchanged_records: int = 0
    input_svtypes: Counter[str] = field(default_factory=Counter)
    output_svtypes: Counter[str] = field(default_factory=Counter)
    converted_by_svtype: Counter[str] = field(default_factory=Counter)
    warnings_by_reason: Counter[str] = field(default_factory=Counter)
    details: list[tuple[str, str, str]] = field(default_factory=list)
    end_overrides: dict[tuple[str, int, str | None], int] = field(default_factory=dict)

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
    """Return a single INFO value as a string."""
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


def record_key(record: pysam.VariantRecord) -> tuple[str, int, str | None]:
    """Stable key used to carry a literal END value through serialization."""
    return (record.contig, record.pos, record.id)


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
    parsed: ParsedBreakendAlt,
) -> str | None:
    """
    Extract non-anchor sequence from this record's own breakend ALT.

    If local sequence is before the bracket, the REF anchor is expected at the
    beginning. If the remote locus comes first, the REF anchor is expected at
    the end.
    """
    sequence = parsed.local_sequence
    ref = record.ref

    if parsed.remote_first:
        if not sequence.endswith(ref):
            return None
        return sequence[:-len(ref)] if ref else sequence

    if not sequence.startswith(ref):
        return None
    return sequence[len(ref):]


def convert_symbolic_record(
    record: pysam.VariantRecord,
) -> tuple[pysam.VariantRecord | None, int | None, str | None, str | None]:
    """Convert one DEL/DUP/INV/INS record without consulting its mate."""
    svtype = scalar_info(record, "SVTYPE")
    if svtype not in CONVERTIBLE_TYPES:
        return None, None, "not_convertible", f"SVTYPE {svtype!r} is not convertible"

    alt_string = get_single_alt(record)
    if alt_string is None:
        return None, None, "alt_count_invalid", "record does not have exactly one ALT allele"

    parsed = parse_breakend_alt(alt_string)
    if parsed is None:
        return None, None, "alt_parse_failure", "ALT allele is not a parseable breakend"

    # DEL/DUP/INV/INS are expected to describe an intrachromosomal interval.
    # Interchromosomal breakends are left untouched rather than turning them
    # into misleading symbolic interval records.
    if parsed.remote_chrom != record.contig:
        return (
            None,
            None,
            "interchromosomal_record",
            f"ALT points to {parsed.remote_chrom}:{parsed.remote_pos}",
        )

    junction_sequence = extract_non_anchor_sequence(record, parsed)
    if junction_sequence is None:
        return (
            None,
            None,
            "anchor_not_identified",
            "the local REF anchor could not be identified in the ALT allele",
        )

    output = record.copy()
    output.alts = (f"<{svtype}>",)

    # MATEID is deliberately retained. It no longer controls transformation,
    # but keeping it preserves the original relationship/provenance between
    # independently emitted ESVEE records.
    #
    # Do not assign VariantRecord.stop here. pysam treats stop as an interval
    # endpoint and cannot faithfully represent a backward-pointing END < POS.
    # The literal remote coordinate is carried separately and forced into the
    # INFO/END field when the record is serialized.

    if junction_sequence:
        output.info["JUNCTIONSEQ"] = junction_sequence
    elif "JUNCTIONSEQ" in output.info:
        del output.info["JUNCTIONSEQ"]

    return output, parsed.remote_pos, None, None


def transform_record(
    record: pysam.VariantRecord,
    report: Report,
) -> pysam.VariantRecord:
    """Transform exactly one record, independently of every other record."""
    svtype = scalar_info(record, "SVTYPE")

    if svtype == "SGL":
        converted = record.copy()
        converted.info["SVTYPE"] = "BND"
        report.sgl_to_bnd += 1
        report.action(record_label(record), "converted SVTYPE=SGL to SVTYPE=BND")
        return converted

    if svtype == "BND" or svtype not in CONVERTIBLE_TYPES:
        report.unchanged_records += 1
        return record

    converted, end_override, reason, message = convert_symbolic_record(record)
    if converted is None:
        assert reason is not None and message is not None
        report.warn(reason, record_label(record), message)
        report.unchanged_records += 1
        return record

    assert end_override is not None
    report.end_overrides[record_key(converted)] = end_override

    report.symbolic_converted += 1
    report.converted_by_svtype[svtype] += 1
    if scalar_info(converted, "JUNCTIONSEQ"):
        report.symbolic_with_junctionseq += 1

    report.action(
        record_label(record),
        f"converted independently to <{svtype}> with END={end_override}",
    )
    return converted


def transform_records(
    records: list[pysam.VariantRecord],
    report: Report,
) -> list[pysam.VariantRecord]:
    """Transform each input record independently and preserve input order."""
    report.input_records = len(records)
    report.input_svtypes.update(
        normalize_svtype(scalar_info(record, "SVTYPE")) for record in records
    )

    output = [transform_record(record, report) for record in records]

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
            description="Remote coordinate parsed from the original breakend ALT",
        )

    if "JUNCTIONSEQ" not in output_header.info:
        output_header.info.add(
            "JUNCTIONSEQ",
            number=1,
            type="String",
            description=(
                "Non-anchor sequence retained from this ESVEE breakend ALT allele; "
                "may include reference-derived assembly sequence and does not "
                "necessarily represent a novel insertion"
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


def force_info_end(record_line: str, end_value: int) -> str:
    """Replace or append INFO/END in one serialized VCF record line."""
    newline = "\n" if record_line.endswith("\n") else ""
    fields = record_line.rstrip("\n").split("\t")
    if len(fields) < 8:
        raise ValueError("serialized VCF record has fewer than 8 columns")

    info_fields = [] if fields[7] in {"", "."} else fields[7].split(";")
    replaced = False
    for index, item in enumerate(info_fields):
        if item.startswith("END="):
            info_fields[index] = f"END={end_value}"
            replaced = True
            break

    if not replaced:
        info_fields.append(f"END={end_value}")

    fields[7] = ";".join(info_fields) if info_fields else f"END={end_value}"
    return "\t".join(fields) + newline


def write_records(
    output_path: Path | None,
    header: pysam.VariantHeader,
    records: Iterable[pysam.VariantRecord],
    end_overrides: dict[tuple[str, int, str | None], int],
) -> None:
    """Write records while preserving literal END values, including END < POS."""
    try:
        serialized_records: list[str] = []
        for record in records:
            line = str(record)
            end_override = end_overrides.get(record_key(record))
            if end_override is not None:
                line = force_info_end(line, end_override)
            serialized_records.append(line)

        content = str(header) + "".join(serialized_records)

        if output_path is None:
            sys.stdout.write(content)
        elif output_path.suffix == ".gz":
            with pysam.BGZFile(str(output_path), "w") as handle:
                handle.write(content.encode())
        else:
            output_path.write_text(content)
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
        ("summary", "records_removed", str(report.input_records - report.output_records)),
        ("summary", "sgl_converted_to_bnd", str(report.sgl_to_bnd)),
        ("summary", "symbolic_records_converted", str(report.symbolic_converted)),
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

    for svtype, count in sorted(report.converted_by_svtype.items()):
        rows.append(("converted_by_svtype", svtype, str(count)))

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
    help="Append per-record actions and warnings to the report. Requires --report.",
)
def main(
    input_vcf: Path | None,
    output_vcf: Path | None,
    report_path: Path | None,
    verbose_report: bool,
) -> None:
    """
    Reformat ESVEE SV records independently without merging mate pairs.

    INPUT_VCF may be plain text or BGZF/gzip-compressed. When omitted, input
    is read from stdin.
    """
    if verbose_report and report_path is None:
        raise click.UsageError("--verbose-report requires --report")

    header, records = read_records(input_vcf)

    report = Report()
    transformed = transform_records(records, report)
    write_records(output_vcf, header, transformed, report.end_overrides)

    if report_path is not None:
        write_report(report_path, report, verbose_report)


if __name__ == "__main__":
    main()
