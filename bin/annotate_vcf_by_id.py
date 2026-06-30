#!/usr/bin/env python3

from collections import defaultdict

import click
import pysam
import csv
import ast

"""
This script aims to annotate the info field of a vcf based on the ID field, using:
1. a VCF file following the VCF specifications, where the ID field will be used to annotate entries
2. a TSV file where the ids to annotate is always the first column and the remaining are new INFO fields: ex.: ID, TAG1, TAG2, TAG3
3. a TXT file with the header lines of the tags added, following the VCF specifications (i.e. ##INFO=<ID=ID,Number=number,Type=type,Description="description")

NOTE: the script adapts the annotation style to the Types in the header.
This means that the ‘Flag’ type indicates that the INFO field does not contain a Value entry, and hence the Number should be
0 and the flag will be added as is to the INFO field based on 0 (false) or 1 (true).
"""

def prepare_header(
    vcf_in: pysam.VariantFile, header_file_path: click.Path
) -> tuple[pysam.VariantHeader, dict[str, str]]:
    """ Adds new header lines to original header from header file, and returns the merged header"""

    # add new header lines to original header from header file
    with open(header_file_path, "r") as header_file:
        new_header_lines = header_file.readlines()

    for line in new_header_lines:
        vcf_in.header.add_line(line)
    merged_header = vcf_in.header

    # create a dictionary to store the type of each header ID
    header_id_type = {}
    for name, meta in merged_header.info.items():
        header_id_type[name] = meta.type

    return merged_header, header_id_type

def string_to_object(string: str) -> str | int | float:
    """Convert a string to the appropriate variable type (int, float, bool, etc.). If the string cannot be converted, return the original string."""
    try:
        return ast.literal_eval(string)
    except:
        return string

def format_annotation_dict(raw_annotation_dict: dict) -> dict:
    vcf_id_key: str = next(iter(raw_annotation_dict))
    vcf_id: str = raw_annotation_dict.pop(vcf_id_key)
    info_fields_dict: dict = {key: string_to_object(value) for key, value in raw_annotation_dict.items()}
    return {vcf_id: info_fields_dict}

"""def merge_dict(dicts, header_id_type: dict[str, str]):
    grouped = {}
    for dict in dicts:
        for key in dict:
            if key not in grouped:
                grouped[key] = {field: [value] for field, value in dict[key].items()}
            else:
                grouped[key] = {field: grouped[key].get(field, []) + [value] for field, value in dict[key].items() }

    for key, fields in grouped.items():
        for field, values in fields.items():
            if len(values) == 1:
                continue
            else:
                if header_id_type.get(field) == "Flag":
                    grouped[key][field] = max(values)
                else:
                    grouped[key][field] = tuple(values)

    return grouped"""

def merge_dict_by_field(dicts: list[dict], field_type: dict[str, str]):
    """Merges a list of dictionaries into a single dictionary, grouping values by key, based on field type. If a field is of type "Flag", the maximum
    value is kept; otherwise, all values are stored in a tuple.

    Keyword arguments:
    dicts -- list of dictionaries to merge i.e. [{vcf_id: {field: value, ...}}, {vcf_id: {field: value, ...}}, ...]
    field_type -- dictionary mapping field names to their types i.e. {field_name: "Flag" | "Other", ...}
    """

    merged_dict = {}
    for dict in dicts:
        vcf_id: str = next(iter(dict))

        if vcf_id not in merged_dict:
            merged_dict[vcf_id] = dict[vcf_id]
        else:
            for field, value in dict[vcf_id].items():
                if field_type[field] == "Flag":
                    merged_dict[vcf_id][field] = max(merged_dict[vcf_id][field], value)
                else:
                    # Append to existing list or create list with existing value and new value
                    if isinstance(merged_dict[vcf_id][field], list):
                        merged_dict[vcf_id][field].append(value)
                    else:
                        merged_dict[vcf_id][field] = [merged_dict[vcf_id][field], value]
    return merged_dict


def get_dict_from_tsv(tsv_file_path: click.Path, header_id_type: dict[str, str]) -> dict:
    """Parses a tsv file into a dictionary with format: {id: [{field: value, ...},
    {field: value, ...}, ...]},
    compatible with annotation with pysam"""
    with open(tsv_file_path, newline="") as f:
        raw_annotation_dict_list: list[dict] = [row for row in csv.DictReader(f, delimiter="\t")]
        annotation_dict_list: list[dict] = [format_annotation_dict(raw_annotation_dict) for raw_annotation_dict in raw_annotation_dict_list]
        converted_dict: dict = merge_dict_by_field(annotation_dict_list, header_id_type)
        print(converted_dict)

        return converted_dict

def annotate_vcf(
    vcf_in: pysam.VariantFile,
    tsv_dict: dict[str, dict[str, str | int | float | tuple]],
    vcf_out: pysam.VariantFile,
) :
    """Annotate VCF file using the tsv_dict and header_file, and save the annotated VCF to output_file."""

    # need to loop through records to get entries with correct ID
    for record in vcf_in:
        if record.id in tsv_dict:
            record.info.update(tsv_dict[record.id])
        vcf_out.write(record)

    return vcf_out

@click.command()

@click.option(
    "-v",
    "--vcf_file",
    type=click.Path(exists=True),
    help="The VCF file to annotate",
    required=True,
)
@click.option(
    "-t",
    "--tsv_file",
    type=click.Path(exists=True),
    help="The TSV file with annotations",
    required=True,
)
@click.option(
    "-h",
    "--header_file",
    type=click.Path(exists=True),
    help="A txt file with new header lines for the annotations",
    required=True,
)
@click.option(
    "-o",
    "--output_file",
    type=click.STRING,
    help="Output name or path for the annotated VCF file",
    required=True,
)

# run
def main(
    vcf_file: click.Path,
    header_file: click.Path,
    tsv_file: click.Path,
    output_file: click.Path,
) -> None:

    # prepare header
    vcf_in: pysam.VariantFile = pysam.VariantFile(vcf_file, "rb")
    merged_header, header_id_type = prepare_header(vcf_in, header_file_path=header_file)

    # parse tsv file
    tsv_dict: dict = get_dict_from_tsv(tsv_file_path=tsv_file, header_id_type=header_id_type)

    # annotate vcf
    # pysam auto-detects compression from file extension
    vcf_out: pysam.VariantFile = pysam.VariantFile(
        output_file, "w", header=merged_header
    )
    annotate_vcf(vcf_in, tsv_dict, vcf_out)

    vcf_in.close()
    vcf_out.close()

    # write index file for annotated vcf
    output_index: pysam.TabixFile = pysam.tabix_index(
        output_file, preset="vcf", force=True
    )


# main
if __name__ == "__main__":
    main()
