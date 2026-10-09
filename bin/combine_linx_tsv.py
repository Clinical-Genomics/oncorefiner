#!/usr/bin/env python3

import click
import pandas as pd

"""
This script aims to combine the three linx files (fusions, breakends, svs) into one merged tsv file
Expects three linx tsv files with a header line and tab-delimited columns:
- Fusion TSV file
- Breakend TSV file

At least the following columns are expected in the input files:
- Fusion TSV file: fivePrimeBreakendId, threePrimeBreakendId, name, reported
- Breakend TSV file: id, svId

The output is a merged tsv file with the following columns:
- vcfId: the ID of the SV in the VCF file, from the vcfId column in the breakends TSV file
- FUSION_NAME: the name of the fusion, from the name column in the fusion TSV file
- REPORTED: whether the fusion is reported, from the reported column in the fusion TSV file, converted to 0/1 for compatibility with pysam specifications
"""

__version__ = "2.0"

def merge_linx_dataframes(
    fusions: pd.DataFrame
) -> pd.DataFrame:
    """Merges input tsv files and returns a dataframe with the relevant columns for annotation of VCF file"""

    # rename columns to desired header names in subsequenct VCF annotation step
    fusions.rename(columns={"name": "FUSION_NAME", "reported": "REPORTED"}, inplace=True)

    # change reported column to 0/1 to adhere to pysam specifications
    fusions["REPORTED"] = fusions["REPORTED"].replace({"false": 0, "true": 1})

    result: pd.DataFrame = fusions[
            ["fivePrimeVcfId",  "threePrimeVcfId", "FUSION_NAME", "REPORTED"]
        ]

    """
    # merge fusions and breakends on 'fivePrimeBreakendId' & 'threePrimeBreakendId' and 'id' columns
    fusion_fivebreakend: pd.DataFrame = fusions.merge(breakends, left_on="fivePrimeBreakendId", right_on="id", how="left")
    fusion_threebreakend: pd.DataFrame = fusions.merge(breakends, left_on="threePrimeBreakendId", right_on="id", how="left")

    # concatenate fusion_fivebreakend and fusion_threebreakend to get all fusions
    fusion_breakend_merge: pd.DataFrame = pd.concat(
        [fusion_fivebreakend, fusion_threebreakend], ignore_index=True
    )


    # keep only the relevant columns for the final output - for debug: add columns svId, fivePrimeBreakendId, threePrimeBreakendId
    result: pd.DataFrame = fusion_breakend_merge[
        ["vcfIdfivePrimeVcfId",  "threePrimeVcfId", "FUSION_NAME", "REPORTED"]
    ]

    # remove duplicates
    result = result.drop_duplicates()
    """
    return result


@click.command()
@click.option(
    "-f",
    "--fusion_file",
    type=click.Path(exists=True),
    help="The fusion linx tsv file",
    required=True,
)
@click.option(
    "-o",
    "--output_file",
    type=click.STRING,
    help="name of output tsv file",
    required=True,
)
@click.version_option(__version__, prog_name="combine_linx_tsv")

def combine_linx_files(
    fusion_file: click.Path,
    output_file: click.STRING,
) -> None:

    # load tsv into pandas df
    fusions: pd.DataFrame = pd.read_csv(fusion_file, sep="\t", dtype=str)

    # merge dataframes
    merged_df: pd.DataFrame = merge_linx_dataframes(fusions)

    # save the final merged dataframe to a tsv file
    merged_df.to_csv(output_file, sep="\t", index=False)


# main
if __name__ == "__main__":
    combine_linx_files()
