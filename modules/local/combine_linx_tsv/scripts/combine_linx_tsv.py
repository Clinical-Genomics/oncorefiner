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

The output is a merged tsv file with the following columns:
- vcfId: the ID of the SV in the VCF file, from the vcfId column in the SV TSV file
- FUSION_NAME: the name of the fusion, from the name column in the fusion TSV file
- REPORTED: whether the fusion is reported, from the reported column in the fusion TSV file, converted to 0/1 for compatibility with pysam specifications
- Columns svId, fivePrimeBreakendId, threePrimeBreakendId: Linx IDs for SVs, 5' and 3' breakends. Added for debugging purpose.
"""


def rename_columns(df: pd.DataFrame, rename_dict: dict[str, str]) -> pd.DataFrame:
    """Renames columns in a dataframe based on a provided mapping dictionary"""

    df.rename(columns=rename_dict, inplace=True)

    return df


def replace_values(
    df: pd.DataFrame, column: str, replace_dict: dict[str, int]
) -> pd.DataFrame:
    """Replaces values in a specified column of a dataframe based on a provided mapping dictionary"""

    df[column] = df[column].replace(replace_dict)

    return df


def merge_dataframes(
    df1: pd.DataFrame, df2: pd.DataFrame, left_on: str, right_on: str, how: str
) -> pd.DataFrame:
    """Merges two dataframes based on specified columns and merge method"""

    merged_df = df1.merge(df2, left_on=left_on, right_on=right_on, how=how)

    return merged_df


def merge_linx_files(
    fusion_file: click.Path, breakend_file: click.Path, sv_file: click.Path
) -> pd.DataFrame:
    """Merges input tsv files and returns a dataframe with the relevant columns for annotation of VCF file"""

    # load tsv into pandas df
    fusions: pd.DataFrame = pd.read_csv(fusion_file, sep="\t", dtype=str)
    breakends: pd.DataFrame = pd.read_csv(breakend_file, sep="\t", dtype=str)
    svs: pd.DataFrame = pd.read_csv(sv_file, sep="\t", dtype=str)

    # rename columns to desired header names in subsequenct VCF annotation step
    fusions = rename_columns(fusions, {"name": "FUSION_NAME", "reported": "REPORTED"})

    # change reported column to 0/1 to adhere to pysam specifications

    fusions = replace_values(fusions, "REPORTED", {"false": 0, "true": 1})

    # merge fusions and breakends on 'fivePrimeBreakendId' & 'threePrimeBreakendId' and 'id' columns
    fusion_fivebreakend: pd.DataFrame = merge_dataframes(
        fusions, breakends, left_on="fivePrimeBreakendId", right_on="id", how="left"
    )

    fusion_fivebreakend: pd.DataFrame = merge_dataframes(
        fusions, breakends, left_on="fivePrimeBreakendId", right_on="id", how="left"
    )

    fusion_threebreakend: pd.DataFrame = merge_dataframes(
        fusions, breakends, left_on="threePrimeBreakendId", right_on="id", how="left"
    )

    # concatenate fusion_fivebreakend and fusion_threebreakend to get all fusions
    fusion_breakend_merge: pd.DataFrame = pd.concat(
        [fusion_fivebreakend, fusion_threebreakend], ignore_index=True
    )

    # merge fusion_breakend_merge with svs on 'svId' column in sv file
    fusion_breakend_svs_merge: pd.DataFrame = merge_dataframes(
        fusion_breakend_merge, svs, left_on="svId", right_on="svId", how="left"
    )

    # keep only the relevant columns for the final output - for debug: add columns svId, fivePrimeBreakendId, threePrimeBreakendId
    result: pd.DataFrame = fusion_breakend_svs_merge[
        ["vcfId", "FUSION_NAME", "REPORTED"]
    ]

    # remove duplicates
    result = result.drop_duplicates()

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
    "-b",
    "--breakend_file",
    type=click.Path(exists=True),
    help="The breakends linx tsv file",
    required=True,
)
@click.option(
    "-sv",
    "--sv_file",
    type=click.Path(exists=True),
    help="The svs linx tsv file",
    required=True,
)
@click.option(
    "-o",
    "--output_file",
    type=click.STRING,
    help="name of output tsv file",
    required=True,
)
def combine_linx_files(
    fusion_file: click.Path,
    breakend_file: click.Path,
    sv_file: click.Path,
    output_file: click.STRING,
) -> None:

    # run main function
    merged_df: pd.DataFrame = merge_linx_files(fusion_file, breakend_file, sv_file)

    # save the final merged dataframe to a tsv file
    merged_df.to_csv(output_file, sep="\t", index=False)


# main
if __name__ == "__main__":
    combine_linx_files()
