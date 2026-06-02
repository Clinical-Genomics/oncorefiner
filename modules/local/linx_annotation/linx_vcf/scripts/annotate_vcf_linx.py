import sys
import pysam
import csv

"""
This script aims to annotate a vcf using:
1. a VCF file following the specifications, where the ID field will be used to annoate entries
2. a TSV file where the ids to annotate is always the first column and the remaining are new INFO fields: ex.: ID, TAG1, TAG2, TAG3
3. a TXT file with the header lines of the tags added, following the VCF specifications (i.e. ##INFO=<ID=ID,Number=number,Type=type,Description="description")

NOTE: the script adapts the annotation style to the Types in the header.
This means that the ‘Flag’ type indicates that the INFO field does not contain a Value entry, and hence the Number should be
0 and the flag will be added as is to the INFO field based on 0 (false) or 1 (true).
"""

# get input files from command line arguments, TODO add click
vcf_file = sys.argv[1]
tsv_file = sys.argv[2]
header_file = sys.argv[3]
output_file = sys.argv[4]


# parse new header
def new_header(header_file, header_lines):
    f= open(header_file, 'rt')
    for line in f:
        line = line.strip()

        # append line to list to be added to new header
        header_lines.append(line)

    return header_lines

# parse tsv file
def parse_tsv(tsv_file):
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

    # convert to tuple if multiple values for same field, to be compatible with VCF specifications
    for key, fields in tsv_dict.items():
        for field, value in fields.items():
            if isinstance(value, list):
                fields[field] = tuple(value)

    return tsv_dict

def annotate_vcf(vcf, new_header_lines, tsv_dict):
    vcf_in = pysam.VariantFile(vcf, 'rb')

    # add new header lines)
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

# run

# parse new header
new_header_lines = []
new_header_lines = new_header(header_file, new_header_lines)

# parse tsv file
tsv_dict = parse_tsv(tsv_file)

# annotate vcf
annotate_vcf(vcf_file, new_header_lines, tsv_dict)
