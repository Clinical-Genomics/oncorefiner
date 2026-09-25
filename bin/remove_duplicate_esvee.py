#!/usr/bin/env python3

"""
Collapse paired ESVEE interval-SV records into one canonical VCF record.

This script is intended to run after the ESVEE SV reformatting step.

For DEL, DUP, INV and INS records:
    - records are grouped by INFO/SVID;
    - reciprocal mate pairs are validated;
    - the leftmost breakpoint is retained as POS;
    - END is set to the rightmost breakpoint;
    - MATEID is removed;
    - the removed mate ID is stored as COLLAPSED_MATEID;
    - breakpoint-specific INFO values from both records are retained in BP_*
      INFO fields, ordered as POS,END;
    - sample-level FORMAT evidence from both breakpoints is retained in
      BP_SAMPLE_* INFO fields, also ordered as POS,END;
    - the original FORMAT/sample column from the POS record is retained.

Existing BND records are never collapsed.

The script is deliberately conservative. Structurally inconsistent pairs cause
an error rather than being silently merged.

The transfer of FORMAT values into INFO currently requires exactly one sample,
because INFO is variant-level and cannot unambiguously represent multiple
samples without encoding sample identity.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Iterable

import click
import pysam


__version__ = "1.1.0"


###############################################################################
# Configuration
###############################################################################

COLLAPSIBLE_SVTYPES = {
    "DEL",
    "DUP",
    "INV",
    "INS",
}


# INFO fields describing evidence local to an individual breakpoint.
#
# These are retained from both original records in BP_* fields.
#
# For scalar fields:
#
#     BP_REF=left_value,right_value
#
# For fields that are already arrays, values from the two breakpoints are
# concatenated. For example:
#
#     CIPOS=-3,4
#     CIPOS=-3,4
#
# becomes:
#
#     BP_CIPOS=-3,4,-3,4
#
# Flag fields such as LINE are represented as 0/1 for each breakpoint:
#
#     BP_LINE=0,1
#
BREAKPOINT_INFO_FIELDS = (
    "ASMID",
    "ASMLEN",
    "ASMSEG",
    "AVGLEN",
    "BEAOR",
    "BEAPOS",
    "BEOR",
    "CIPOS",
    "DF",
    "IHOMPOS",
    "LINE",
    "MLR",
    "REF",
    "REFPAIR",
    "SEGALEN",
    "SEGID",
    "SEGMAPQ",
    "SEGRL",
    "SEGSCO",
    "SF",
    "UFP",
    "VF",
)


# INFO fields that already contain values for both breakpoints.
#
# In the original ESVEE records these are normally represented in opposite
# order between the two mates. For example:
#
#     record 1: PURPLE_AF=0.327,0.201
#     record 2: PURPLE_AF=0.201,0.327
#
# The script validates this relationship and retains the ordering from the
# canonical POS record.
#
PAIRED_INFO_FIELDS = (
    "PURPLE_AF",
    "PURPLE_CN",
    "PURPLE_CN_CHANGE",
)


# Scalar FORMAT/sample fields to preserve from both breakpoints.
#
# These are transferred to INFO as:
#
#     BP_SAMPLE_AF=value_at_POS,value_at_END
#
BREAKPOINT_SAMPLE_FIELDS = (
    "AF",
    "DF",
    "DP",
    "REF",
    "REFPAIR",
    "SB",
    "SF",
    "VF",
)


# VCF type for each sample field when transferred into INFO.
BREAKPOINT_SAMPLE_FIELD_TYPES = {
    "AF": "Float",
    "DF": "Integer",
    "DP": "Integer",
    "REF": "Integer",
    "REFPAIR": "Integer",
    "SB": "Float",
    "SF": "Integer",
    "VF": "Integer",
}


###############################################################################
# Exceptions
###############################################################################

class PairError(click.ClickException):
    """Raised when two records cannot safely be collapsed."""


###############################################################################
# Basic record helpers
###############################################################################

def record_label(record: pysam.VariantRecord) -> str:
    """Return a concise human-readable identifier for one VCF record."""

    record_id = record.id if record.id is not None else "."

    return f"{record.contig}:{record.pos}:{record_id}"


def get_svtype(record: pysam.VariantRecord) -> str | None:
    """Return INFO/SVTYPE when it contains a single string."""

    value = record.info.get("SVTYPE")

    return value if isinstance(value, str) else None


def get_svid(record: pysam.VariantRecord) -> str | None:
    """Return INFO/SVID when it contains a single string."""

    value = record.info.get("SVID")

    return value if isinstance(value, str) else None


def get_mateid(record: pysam.VariantRecord) -> str | None:
    """Return INFO/MATEID when it contains a single string."""

    value = record.info.get("MATEID")

    return value if isinstance(value, str) else None


###############################################################################
# Value helpers
###############################################################################

def normalise_value(value):
    """
    Convert array-like pysam INFO values into ordinary tuples.

    Scalar values are returned unchanged.
    """

    if value is None:
        return None

    if isinstance(value, tuple):
        return tuple(value)

    if isinstance(value, list):
        return tuple(value)

    return value


def values_equal(a, b) -> bool:
    """Compare INFO values after normalising array-like values."""

    return normalise_value(a) == normalise_value(b)


def reversed_value(value):
    """Reverse an array-like INFO value."""

    value = normalise_value(value)

    if not isinstance(value, tuple):
        return None

    return tuple(reversed(value))


def flatten_breakpoint_values(left_value, right_value):
    """
    Combine values from the POS and END breakpoints.

    Examples
    --------

    Scalar + scalar:

        103
        166

    becomes:

        (103, 166)

    Tuple + tuple:

        (-3, 4)
        (-3, 4)

    becomes:

        (-3, 4, -3, 4)
    """

    result = []

    for value in (left_value, right_value):

        if value is None:
            continue

        value = normalise_value(value)

        if isinstance(value, tuple):
            result.extend(value)
        else:
            result.append(value)

    if not result:
        return None

    return tuple(result)


###############################################################################
# Pair validation
###############################################################################

def validate_pair(
    first: pysam.VariantRecord,
    second: pysam.VariantRecord,
) -> None:
    """Validate that two records form one reciprocal interval-SV pair."""

    first_svtype = get_svtype(first)

    if first.contig != second.contig:
        raise PairError(
            f"SVID={get_svid(first)} is SVTYPE={first_svtype}, but its "
            f"records occur on different chromosomes: "
            f"{record_label(first)} and {record_label(second)}"
        )

    if first.id is None or second.id is None:
        raise PairError(
            f"SVID={get_svid(first)} cannot be collapsed because both "
            f"records must have an ID"
        )

    first_mate = get_mateid(first)
    second_mate = get_mateid(second)

    if first_mate != second.id:
        raise PairError(
            f"{record_label(first)} has MATEID={first_mate!r}; "
            f"expected {second.id!r}"
        )

    if second_mate != first.id:
        raise PairError(
            f"{record_label(second)} has MATEID={second_mate!r}; "
            f"expected {first.id!r}"
        )

    first_svid = get_svid(first)
    second_svid = get_svid(second)

    if first_svid != second_svid:
        raise PairError(
            f"Records {record_label(first)} and "
            f"{record_label(second)} have different SVID values"
        )


def validate_paired_fields(
    left: pysam.VariantRecord,
    right: pysam.VariantRecord,
) -> None:
    """
    Validate INFO fields that already describe both breakpoints.

    ESVEE/PURPLE normally writes these fields in opposite order on the two
    mate records.
    """

    for field in PAIRED_INFO_FIELDS:

        left_value = left.info.get(field)
        right_value = right.info.get(field)

        if left_value is None and right_value is None:
            continue

        if left_value is None or right_value is None:
            raise PairError(
                f"SVID={get_svid(left)}: INFO/{field} is present on only "
                f"one member of the pair"
            )

        expected = reversed_value(right_value)

        if expected is None:
            raise PairError(
                f"SVID={get_svid(left)}: INFO/{field} was expected to "
                f"contain paired breakpoint values, but found "
                f"{normalise_value(left_value)!r} and "
                f"{normalise_value(right_value)!r}"
            )

        if not values_equal(left_value, expected):
            raise PairError(
                f"SVID={get_svid(left)}: INFO/{field} values are "
                f"inconsistent between mates: "
                f"{normalise_value(left_value)!r} versus "
                f"{normalise_value(right_value)!r}"
            )


###############################################################################
# Header construction
###############################################################################

def add_output_headers(
    header: pysam.VariantHeader,
) -> pysam.VariantHeader:
    """Create the output VCF header with fields used by this script."""

    output_header = header.copy()

    ###########################################################################
    # Pair provenance
    ###########################################################################

    if "COLLAPSED_MATEID" not in output_header.info:
        output_header.info.add(
            "COLLAPSED_MATEID",
            number=1,
            type="String",
            description=(
                "ID of the mate breakend collapsed into this interval SV record"
            ),
        )

    ###########################################################################
    # Breakpoint-specific INFO fields
    ###########################################################################

    for field in BREAKPOINT_INFO_FIELDS:

        if field not in output_header.info:
            continue

        bp_field = f"BP_{field}"

        if bp_field in output_header.info:
            continue

        original = output_header.info[field]

        # VCF Flag fields do not contain a value. Store presence/absence
        # at the POS and END breakpoints as 1/0.
        if original.type == "Flag":

            number = 2
            field_type = "Integer"

            description = (
                f"Presence (1) or absence (0) of INFO/{field} at the POS "
                f"and END breakpoints, in that order"
            )

        # Scalar original fields become exactly two values:
        #
        #     POS,END
        #
        elif original.number == 1:

            number = 2
            field_type = original.type

            description = (
                f"Original INFO/{field} value at the POS and END "
                f"breakpoints, in that order"
            )

        # Fields such as CIPOS already contain multiple values at each
        # breakpoint. Those arrays are concatenated.
        else:

            number = "."
            field_type = original.type

            description = (
                f"Original INFO/{field} values from the POS breakpoint "
                f"followed by values from the END breakpoint"
            )

        output_header.info.add(
            bp_field,
            number=number,
            type=field_type,
            description=description,
        )

    ###########################################################################
    # Sample FORMAT values transferred into INFO
    ###########################################################################

    for field in BREAKPOINT_SAMPLE_FIELDS:

        info_field = f"BP_SAMPLE_{field}"

        if info_field in output_header.info:
            continue

        output_header.info.add(
            info_field,
            number=2,
            type=BREAKPOINT_SAMPLE_FIELD_TYPES[field],
            description=(
                f"Original sample FORMAT/{field} value at the POS and END "
                f"breakpoints, in that order"
            ),
        )

    ###########################################################################
    # FORMAT/AD
    ###########################################################################

    # AD is Number=R and therefore already contains REF and ALT depth.
    # Preserve the two components independently to avoid an ambiguous
    # four-value field.

    if "BP_SAMPLE_AD_REF" not in output_header.info:

        output_header.info.add(
            "BP_SAMPLE_AD_REF",
            number=2,
            type="Integer",
            description=(
                "Reference allele depth from FORMAT/AD at the POS and END "
                "breakpoints, in that order"
            ),
        )

    if "BP_SAMPLE_AD_ALT" not in output_header.info:

        output_header.info.add(
            "BP_SAMPLE_AD_ALT",
            number=2,
            type="Integer",
            description=(
                "Alternative allele depth from FORMAT/AD at the POS and END "
                "breakpoints, in that order"
            ),
        )

    return output_header


###############################################################################
# INFO transfer
###############################################################################

def transfer_breakpoint_info(
    collapsed: pysam.VariantRecord,
    left: pysam.VariantRecord,
    right: pysam.VariantRecord,
) -> None:
    """
    Transfer breakpoint-specific INFO values from both original records.

    Values always follow genomic interval order:

        POS,END

    VCF Flag fields are represented as 0/1 for each breakpoint.
    """

    for field in BREAKPOINT_INFO_FIELDS:

        bp_field = f"BP_{field}"

        if bp_field not in collapsed.header.info:
            continue

        original = collapsed.header.info[field]

        # For VCF Flag fields, retain whether the flag was present at each
        # breakpoint as 0/1.
        if original.type == "Flag":
            collapsed.info[bp_field] = (
                int(field in left.info),
                int(field in right.info),
            )
            continue

        left_value = left.info.get(field)
        right_value = right.info.get(field)

        combined = flatten_breakpoint_values(
            left_value,
            right_value,
        )

        if combined is not None:
            collapsed.info[bp_field] = combined


###############################################################################
# Sample FORMAT -> INFO transfer
###############################################################################

def transfer_sample_info(
    collapsed: pysam.VariantRecord,
    left: pysam.VariantRecord,
    right: pysam.VariantRecord,
) -> None:
    """
    Transfer sample-level FORMAT evidence from both breakpoints into INFO.

    The resulting fields are ordered:

        POS,END

    For example:

        BP_SAMPLE_AF=0.27,0.177
        BP_SAMPLE_DP=200,305

    FORMAT/AD is separated into:

        BP_SAMPLE_AD_REF
        BP_SAMPLE_AD_ALT

    The original FORMAT/sample column on the retained POS record remains
    unchanged.

    Because INFO is variant-level, this representation requires exactly one
    sample.
    """

    samples = list(collapsed.header.samples)

    if not samples:
        return

    if len(samples) != 1:
        raise PairError(
            f"SVID={get_svid(collapsed)}: cannot transfer FORMAT values "
            f"into INFO because the VCF contains {len(samples)} samples. "
            f"This representation requires exactly one sample."
        )

    sample_name = samples[0]

    left_sample = left.samples[sample_name]
    right_sample = right.samples[sample_name]

    ###########################################################################
    # Scalar FORMAT fields
    ###########################################################################

    for field in BREAKPOINT_SAMPLE_FIELDS:

        if field not in collapsed.header.formats:
            continue

        left_value = left_sample.get(field)
        right_value = right_sample.get(field)

        if left_value is None and right_value is None:
            continue

        if left_value is None or right_value is None:
            raise PairError(
                f"SVID={get_svid(collapsed)}: FORMAT/{field} is present "
                f"at only one breakpoint for sample {sample_name!r}"
            )

        # These fields are expected to be scalar.
        if isinstance(left_value, (tuple, list)):
            raise PairError(
                f"SVID={get_svid(collapsed)}: FORMAT/{field} at POS was "
                f"expected to be scalar but found {left_value!r}"
            )

        if isinstance(right_value, (tuple, list)):
            raise PairError(
                f"SVID={get_svid(collapsed)}: FORMAT/{field} at END was "
                f"expected to be scalar but found {right_value!r}"
            )

        collapsed.info[f"BP_SAMPLE_{field}"] = (
            left_value,
            right_value,
        )

    ###########################################################################
    # FORMAT/AD
    ###########################################################################

    if "AD" not in collapsed.header.formats:
        return

    left_ad = left_sample.get("AD")
    right_ad = right_sample.get("AD")

    if left_ad is None and right_ad is None:
        return

    if left_ad is None or right_ad is None:
        raise PairError(
            f"SVID={get_svid(collapsed)}: FORMAT/AD is present at only "
            f"one breakpoint for sample {sample_name!r}"
        )

    if len(left_ad) != 2:
        raise PairError(
            f"SVID={get_svid(collapsed)}: expected FORMAT/AD at POS to "
            f"contain exactly REF,ALT depths, but found {left_ad!r}"
        )

    if len(right_ad) != 2:
        raise PairError(
            f"SVID={get_svid(collapsed)}: expected FORMAT/AD at END to "
            f"contain exactly REF,ALT depths, but found {right_ad!r}"
        )

    collapsed.info["BP_SAMPLE_AD_REF"] = (
        left_ad[0],
        right_ad[0],
    )

    collapsed.info["BP_SAMPLE_AD_ALT"] = (
        left_ad[1],
        right_ad[1],
    )


###############################################################################
# Pair collapsing
###############################################################################

def collapse_pair(
    first: pysam.VariantRecord,
    second: pysam.VariantRecord,
) -> pysam.VariantRecord:
    """Collapse two reciprocal interval-SV records into one record."""

    validate_pair(first, second)

    ###########################################################################
    # Determine canonical POS -> END ordering
    ###########################################################################

    if first.pos <= second.pos:
        left = first
        right = second
    else:
        left = second
        right = first

    validate_paired_fields(
        left,
        right,
    )

    ###########################################################################
    # Start with the left/POS record
    ###########################################################################

    collapsed = left.copy()

    # Explicitly canonicalise the interval.
    collapsed.pos = left.pos
    collapsed.stop = right.pos

    ###########################################################################
    # Remove obsolete mate relationship
    ###########################################################################

    if "MATEID" in collapsed.info:
        del collapsed.info["MATEID"]

    collapsed.info["COLLAPSED_MATEID"] = right.id

    ###########################################################################
    # Preserve breakpoint-specific INFO evidence
    ###########################################################################

    transfer_breakpoint_info(
        collapsed=collapsed,
        left=left,
        right=right,
    )

    ###########################################################################
    # Preserve sample-level evidence from both breakpoints
    ###########################################################################

    transfer_sample_info(
        collapsed=collapsed,
        left=left,
        right=right,
    )

    return collapsed


###############################################################################
# Grouping
###############################################################################

def group_records(
    records: list[pysam.VariantRecord],
) -> tuple[
    dict[str, list[pysam.VariantRecord]],
    list[pysam.VariantRecord],
]:
    """
    Separate collapsible interval-SV records from records that should pass
    through unchanged.
    """

    groups = defaultdict(list)
    passthrough = []

    for record in records:

        svtype = get_svtype(record)

        # BND and any other non-interval type are retained unchanged.
        if svtype not in COLLAPSIBLE_SVTYPES:
            passthrough.append(record)
            continue

        # All collapsible interval SVs are expected to have an SVID.
        svid = get_svid(record)

        groups[svid].append(record)

    return groups, passthrough


###############################################################################
# Main record processing
###############################################################################

def process_records(
    records: list[pysam.VariantRecord],
) -> tuple[list[pysam.VariantRecord], list[str]]:

    groups, passthrough = group_records(records)

    output_records = list(passthrough)
    report_lines = []

    for svid, group in groups.items():

        #######################################################################
        # Unexpected variant group size
        #######################################################################

        if len(group) != 2:

            labels = ",".join(
                record_label(record)
                for record in group
            )

            raise PairError(
                f"SVID={svid} occurs {len(group)} times; expected exactly "
                f"two records. Records: {labels}"
            )

        #######################################################################
        # Collapse SV pair to one record
        #######################################################################

        collapsed = collapse_pair(
            group[0],
            group[1],
        )

        output_records.append(collapsed)

        report_lines.append(
            f"collapsed\t{svid}\t"
            f"{record_label(group[0])},{record_label(group[1])}\t"
            f"retained={record_label(collapsed)}"
        )

    ###########################################################################
    # Restore genomic ordering
    ###########################################################################

    if output_records:

        # Map each chromosome to its order in the VCF header.
        chrom_order = {
            contig: index
            for index, contig in enumerate(output_records[0].header.contigs)
        }

        # Sort SVs first by chromosome, then by genomic position.
        output_records.sort(
            key=lambda sv: (
                chrom_order.get(sv.contig, len(chrom_order)),
                sv.pos,
            )
        )

    return output_records, report_lines


###############################################################################
# VCF I/O
###############################################################################

def output_mode(path: Path | None) -> str:
    """Return pysam output mode."""

    if path is None:
        return "w"

    return "wz" if path.suffix == ".gz" else "w"


def read_vcf(
    path: Path,
) -> tuple[
    pysam.VariantHeader,
    list[pysam.VariantRecord],
]:
    """Read VCF and translate records onto the augmented output header."""

    try:

        with pysam.VariantFile(str(path), "r") as reader:

            output_header = add_output_headers(
                reader.header
            )

            records = []

            for record in reader:

                copied = record.copy()

                # VariantRecords are tied to their original header. Translate
                # onto the augmented header before assigning new INFO fields.
                copied.translate(output_header)

                records.append(copied)

            return output_header, records

    except (OSError, ValueError) as error:

        raise click.ClickException(
            f"Failed to read VCF {path}: {error}"
        ) from error


def write_vcf(
    path: Path | None,
    header: pysam.VariantHeader,
    records: Iterable[pysam.VariantRecord],
) -> None:
    """Write VCF or BGZF-compressed VCF."""

    filename = "-" if path is None else str(path)

    try:

        with pysam.VariantFile(
            filename,
            output_mode(path),
            header=header,
        ) as writer:

            for record in records:
                writer.write(record)

    except (OSError, ValueError) as error:

        raise click.ClickException(
            f"Failed to write VCF {filename}: {error}"
        ) from error


###############################################################################
# Report
###############################################################################

def write_report(
    path: Path,
    input_count: int,
    output_count: int,
    report_lines: list[str],
) -> None:
    """Write a TSV summary of pair collapsing."""

    collapsed_count = sum(
        line.startswith("collapsed\t")
        for line in report_lines
    )

    try:

        with path.open("w") as handle:

            handle.write(
                "section\tname\tvalue\n"
            )

            handle.write(
                f"summary\tinput_records\t"
                f"{input_count}\n"
            )

            handle.write(
                f"summary\toutput_records\t"
                f"{output_count}\n"
            )

            handle.write(
                f"summary\trecords_removed\t"
                f"{input_count - output_count}\n"
            )

            handle.write(
                f"summary\tcollapsed_pairs\t"
                f"{collapsed_count}\n"
            )

            handle.write(
                "\naction\tSVID\trecords\tdescription\n"
            )

            for line in report_lines:
                handle.write(line + "\n")

    except OSError as error:

        raise click.ClickException(
            f"Failed to write report {path}: {error}"
        ) from error


###############################################################################
# CLI
###############################################################################

@click.command(
    context_settings={
        "help_option_names": ["-h", "--help"]
    }
)
@click.version_option(
    version=__version__
)
@click.argument(
    "input_vcf",
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
    type=click.Path(
        path_type=Path,
        dir_okay=False,
    ),
    help=(
        "Output .vcf or BGZF-compressed .vcf.gz. "
        "Defaults to stdout."
    ),
)
@click.option(
    "--report",
    "report_path",
    type=click.Path(
        path_type=Path,
        dir_okay=False,
    ),
    help=(
        "Write a TSV report describing collapsed pairs."
    ),
)
def main(
    input_vcf: Path,
    output_vcf: Path | None,
    report_path: Path | None,
) -> None:
    """
    Collapse paired ESVEE DEL/DUP/INV/INS records.

    INPUT_VCF should be the output from the ESVEE SV reformatting script.

    Existing BND records are retained independently.
    """

    header, records = read_vcf(
        input_vcf
    )

    output_records, report_lines = process_records(
        records
    )

    write_vcf(
        output_vcf,
        header,
        output_records,
    )

    if report_path is not None:

        write_report(
            report_path,
            input_count=len(records),
            output_count=len(output_records),
            report_lines=report_lines,
        )


if __name__ == "__main__":
    main()