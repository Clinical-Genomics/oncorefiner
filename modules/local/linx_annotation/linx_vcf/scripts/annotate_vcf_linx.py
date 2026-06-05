import click
import pysam
import csv

"""
This script aims to annotate a vcf using:
1. a VCF file following the VCF specifications, where the ID field will be used to annoate entries
2. a TSV file where the ids to annotate is always the first column and the remaining are new INFO fields: ex.: ID, TAG1, TAG2, TAG3
3. a TXT file with the header lines of the tags added, following the VCF specifications (i.e. ##INFO=<ID=ID,Number=number,Type=type,Description="description")

NOTE: the script adapts the annotation style to the Types in the header.
This means that the ‘Flag’ type indicates that the INFO field does not contain a Value entry, and hence the Number should be
0 and the flag will be added as is to the INFO field based on 0 (false) or 1 (true).
"""


# parse tsv file
def parse_tsv(tsv_file: str) -> dict:

    """Parses the tsv file into a dictionary with the id as key and the remaining columns as values in a nested dict"""

    tsv_dict = {}
    with open(tsv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter="\t") # convert each tsv row to dict

        id_column = reader.fieldnames[0]
        for row in reader:
            key = row.pop(id_column) # removes id column from dict, returns id (value)

            # convert numeric strings to ints
            values = { k: int(v) if v.isdigit() else v
                      for k, v in row.items() }

            # avoid overwriting if id has multiple annotations
            if key not in tsv_dict:
                tsv_dict[key] = values

            else: # append new value to existing annotations
                for field, value in values.items():
                    if tsv_dict[key][field] == value:
                        continue
                    current_entry = tsv_dict[key][field] # current entry for column
                    if isinstance(current_entry, list):
                        current_entry.append(value)
                    else:
                        tsv_dict[key][field] = [current_entry, value]

    # convert to tuple if multiple values for same field, to be compatible with pysam specifications
    for key, fields in tsv_dict.items():
        for field, value in fields.items():
            if isinstance(value, list):
                fields[field] = tuple(value)

    return tsv_dict

def annotate_vcf(vcf_file: str, header_file: str, tsv_dict: dict, output_file: str) -> tuple:

    """ Annotate VCF file using the tsv_dict and header_file, and save the annotated VCF to output_file. Returns the path to the annotated VCF and its index file"""

    vcf_in = pysam.VariantFile(vcf_file, 'rb')

    # add new header lines to original header from header file
    with open(header_file, 'r') as hf:
        new_header_lines = hf.read().splitlines()

    for line in new_header_lines:
        vcf_in.header.add_line(line)
    merged_header = vcf_in.header

    vcf_out = pysam.VariantFile( output_file, "w", header=merged_header ) # auto-detects compression from file extension

    for record in vcf_in:
        if record.id in list(tsv_dict.keys()):
            annotations = tsv_dict[record.id]
            record.info.update(annotations)
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    output_index = pysam.tabix_index(output_file, preset="vcf", force=True)

    return output_file, output_index

@click.command()

# get input files
@click.option("-v", "--vcf_file", type=click.Path(exists=True), help="The VCF file to annotate", required=True)
@click.option("-t", "--tsv_file", type=click.Path(exists=True), help="The TSV file with annotations", required=True)
@click.option("-h", "--header_file", type=click.Path(exists=True), help="A txt file with new header lines for the annotations", required=True)
@click.option("-o", "--output_file", type=click.Path(), help="name of output annotated VCF file", required=True)

# run
def main(vcf_file: str, header_file: str, tsv_file: str, output_file: str) -> None:

    # parse tsv file
    tsv_dict = parse_tsv(tsv_file)

    # annotate vcf
    annotate_vcf(vcf_file, header_file, tsv_dict, output_file)

# main
if __name__ == "__main__":
    main()
