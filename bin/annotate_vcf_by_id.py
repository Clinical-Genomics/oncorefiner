#!/usr/bin/env python3

import click
import pysam
import csv

"""
This script aims to annotate the info field of a vcf based on the ID field, using:
1. a VCF file following the VCF specifications, where the ID field will be used to annotate entries
2. a TSV file where the ids to annotate is always the first column and the remaining are new INFO fields: ex.: ID, TAG1, TAG2, TAG3
3. a TXT file with the header lines of the tags added, following the VCF specifications (i.e. ##INFO=<ID=ID,Number=number,Type=type,Description="description")

NOTE: the script adapts the annotation style to the Types in the header.
This means that the ‘Flag’ type indicates that the INFO field does not contain a Value entry, and hence the Number should be
0 and the flag will be added as is to the INFO field based on 0 (false) or 1 (true).
"""


def convert_data_types(tsv_row_dict: dict[str, str]) -> dict[str, str | int | float ]:
    """ Convert values to float, int or string - pysam requires correct types for annotation"""
    for key, value in tsv_row_dict.items():
        try:
            tsv_row_dict[key] = int(value)
        except ValueError:
            try:
                tsv_row_dict[key] = float(value)
            except ValueError:
                tsv_row_dict[key] = str(value)

    return tsv_row_dict

def convert_entry_to_list(
    tsv_annotations_dict: dict[str, dict[str, str | int | float]],
    vcf_id: str,
    annotation_tag: str,
    annotation_value: str | int | float
) -> dict[str, dict[str, str | int | float]]:
    """ Converts a single annotation entry to a list when multiple values exist for the same ID and tag"""

    # current entry for column
    current_annotation_value: str | int | float | list = tsv_annotations_dict[vcf_id][annotation_tag]

    if isinstance(current_annotation_value, list):
        current_annotation_value.append(annotation_value)
    else:
        tsv_annotations_dict[vcf_id][annotation_tag] = [
            current_annotation_value,
            annotation_value,
        ]

    return tsv_annotations_dict

def add_entry_to_existing_id(
    vcf_id: str,
    vcf_id_annotations: dict[str, str | int | float ],
    tsv_annotations_dict: dict[str, dict[str, str | int | float]],
) -> dict[str, dict[str, str | int | float ]]:
    """ Adds new annotation entries to existing ID in the annotation dictionary,
    converting to list if multiple values exist for the same tag (using function above)"""

    for annotation_tag, annotation_value in vcf_id_annotations.items():
        # if new annotation is the same as existing annotation, skip
        if tsv_annotations_dict[vcf_id][annotation_tag] == annotation_value:
            continue

        tsv_annotations_dict = convert_entry_to_list(tsv_annotations_dict, vcf_id, annotation_tag, annotation_value)

    return tsv_annotations_dict


def add_entry_to_annotation_dict(
    vcf_id: str,
    vcf_id_annotations: dict[str, str | int | float ],
    tsv_annotations_dict: dict[str, dict[str, str | int | float]]
) -> dict[str, dict[str, str | int | float ]]:
    """ Adds annotation entries to the annotation dictionary,
    avoiding overwriting existing annotations for the same ID by using the function above"""

    # avoid overwriting if id has multiple annotations
    if vcf_id not in tsv_annotations_dict:
        tsv_annotations_dict[vcf_id] = vcf_id_annotations

    else:
        # append new annotation value to existing ID annotations
        tsv_annotations_dict = add_entry_to_existing_id(vcf_id, vcf_id_annotations, tsv_annotations_dict )

    return tsv_annotations_dict


def convert_lists_to_tuples(
    tsv_annotations_dict: dict[str, dict[str, str | int | float]]
) -> dict[str, dict[str, str | int | float]]:
    """ Converts lists to tuples in the annotation dictionary when multiple values exist for the same ID and tag,
    to be compatible with pysam specifications"""

    for vcf_id, vcf_id_annotations in tsv_annotations_dict.items():
        for annotation_tag, annotation_value in vcf_id_annotations.items():
            if isinstance(annotation_value, list):
                vcf_id_annotations[annotation_tag] = tuple(annotation_value)

    return tsv_annotations_dict


def update_vcf_header_with_new_lines(
    vcf_in: pysam.VariantFile, header_file_path: click.Path
) -> pysam.VariantHeader:
    """ Adds new header lines to original header from header file, and returns the merged header"""

    # add new header lines to original header from header file
    with open(header_file_path, "r") as hf:
        new_header_lines = hf.read().splitlines()

    for line in new_header_lines:
        vcf_in.header.add_line(line)
    merged_header = vcf_in.header

    return merged_header


def update_vcf_info_field(
    record: pysam.VariantRecord, tsv_dict: dict[str, dict[str, str | int | float | tuple]]
) -> pysam.VariantRecord:
    """Updates the INFO field of a VCF record with annotations from the tsv_dict"""

    if record.id in tsv_dict:
        annotations: dict[str, str | int | float | tuple] = tsv_dict[record.id]
        record.info.update(annotations)
    return record


# parse tsv file
def get_dict_from_tsv(tsv_file_path: click.Path) -> dict:
    """Parses the tsv file into a dictionary (id as key and the remaining columns as values in a nested dict),
    to be compatible with annotation with pysam"""

    tsv_annotations_dict = {}
    with open(tsv_file_path, newline="") as f:
        # converts each tsv row to dict
        tsv_reader: csv.DictReader = csv.DictReader(f, delimiter="\t")

        vcf_id_column_header: str = tsv_reader.fieldnames[0]

        for tsv_row_dict in tsv_reader:
            # removes id column from dict, returns id
            vcf_id: str = tsv_row_dict.pop(vcf_id_column_header) # maybe split up

            # convert values to ints, float or keep as string
            vcf_id_annotations: dict[str, str | int | float] = convert_data_types(
                tsv_row_dict
            )

            # add row to dictionary, avoid overwriting if id has multiple annotations
            tsv_annotations_dict: dict[str, dict[str, str | int | float]] = (
                add_entry_to_annotation_dict(
                    vcf_id, vcf_id_annotations, tsv_annotations_dict
                )
            )

    # convert list to tuple if multiple values for same field
    tsv_annotations_dict: dict[str, dict[str, str | int | float | tuple]] = (
        convert_lists_to_tuples(tsv_annotations_dict)
    )
    return tsv_annotations_dict


def annotate_vcf(
    vcf_file: click.Path,
    header_file: click.Path,
    tsv_dict: dict[str, dict[str, str | int | float | tuple]],
    output_file: click.Path,
) -> tuple:
    """Annotate VCF file using the tsv_dict and header_file, and save the annotated VCF to output_file.
    Returns the path to the annotated VCF and its index file"""

    vcf_in: pysam.VariantFile = pysam.VariantFile(vcf_file, "rb")

    merged_header: pysam.VariantHeader = update_vcf_header_with_new_lines(
        vcf_in, header_file_path=header_file
    )

    # pysam auto-detects compression from file extension
    vcf_out: pysam.VariantFile = pysam.VariantFile(
        output_file, "w", header=merged_header
    )

    # need to loop through records to get entries with correct ID
    for record in vcf_in:
        record: pysam.VariantRecord = update_vcf_info_field(record, tsv_dict)
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

    # write index file for annotated vcf
    output_index: pysam.TabixFile = pysam.tabix_index(
        output_file, preset="vcf", force=True
    )

    return output_file, output_index


@click.command()

# get input files
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
    type=click.Path(),  # change to string
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

    # parse tsv file
    tsv_dict: dict = get_dict_from_tsv(tsv_file_path=tsv_file)

    # annotate vcf
    annotate_vcf(vcf_file, header_file, tsv_dict, output_file)


# main
if __name__ == "__main__":
    main()
