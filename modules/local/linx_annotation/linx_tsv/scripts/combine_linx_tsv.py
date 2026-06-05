import click
import pandas as pd

"""
This script aims to combine the three linx files (fusions, breakends, svs) into one merged tsv file
Expects three linx tsv files with a header line and tab-delimited columns:
- Fusion TSV file
- Breakend TSV file
- SV TSV file

At least the following columns are expected in the input files:
- Fusion TSV file: fivePrimeBreakendId, threePrimeBreakendId, name, reported
- Breakend TSV file: id, svId
- SV TSV file: svId, vcfId
"""

@click.command()
@click.option("-f", "--fusion_file", type=click.Path(exists=True), help="The fusion linx tsv file", required=True)
@click.option("-b", "--breakend_file", type=click.Path(exists=True), help="The breakends linx tsv file", required=True)
@click.option("-sv", "--sv_file", type=click.Path(exists=True), help="The svs linx tsv file", required=True)
@click.option("-o", "--output_file", type=click.Path(), help="name of output tsv file", required=True)

def merge_linx_files(fusion_file: str, breakend_file: str, sv_file: str, output_file: str) -> None:

    # load tsv into pandas df
    fusions   = pd.read_csv(fusion_file, sep='\t', dtype=str)
    breakends = pd.read_csv(breakend_file, sep='\t', dtype=str)
    svs       = pd.read_csv(sv_file, sep='\t', dtype=str)

    fusions.rename(columns={'name': 'FUSION_NAME', 'reported': 'REPORTED'}, inplace=True) # rename columns to desired header names in subsequenct VCF annotation step
    fusions['REPORTED'] = fusions['REPORTED'].replace({'false': 0, 'true': 1}) # to adhere to pysam

    # subset dataframes
    fusions = fusions[['fivePrimeBreakendId', 'threePrimeBreakendId', 'FUSION_NAME', 'REPORTED']]
    breakends = breakends[['id', 'svId']]

    # merge fusions and breakends on 'fivePrimeBreakendId' & 'threePrimeBreakendId' and 'id' columns
    fusion_fivebreakend = fusions.merge(breakends, left_on='fivePrimeBreakendId', right_on='id', how='left')
    fusion_threebreakend = fusions.merge(breakends, left_on='threePrimeBreakendId', right_on='id', how='left')

    # merge fusion_fivebreakend and fusion_threebreakend to get all fusions
    fusion_breakend_merge = pd.concat([fusion_fivebreakend, fusion_threebreakend], ignore_index=True)

    # merge fusion_breakend_merge with svs on 'svId' column in sv file
    fusion_breakend_svs_merge = fusion_breakend_merge.merge(svs, on='svId', how='left')

    # keep only the relevant columns for the final output - for debug: add columns svId, fivePrimeBreakendId, threePrimeBreakendId
    result = fusion_breakend_svs_merge[
        [
            "vcfId",
            "FUSION_NAME",
            "REPORTED"
        ]
    ]

    # remove duplicates
    result = result.drop_duplicates()

    # save the final merged dataframe to a tsv file
    result.to_csv(output_file, sep='\t', index=False)

# main
if __name__ == "__main__":
    merge_linx_files()
