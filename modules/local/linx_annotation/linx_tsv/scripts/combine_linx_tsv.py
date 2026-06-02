import pandas as pd
import argparse

"""
This script aims to combine the three linx files (fusions, breakends, svs) into one merged tsv file
"""

def merge_linx_files(fusions, breakends, svs, output_file):

    fusions.rename(columns={'name': 'FUSION_NAME', 'reported': 'REPORTED'}, inplace=True)
    fusions['REPORTED'] = fusions['REPORTED'].replace({'false': 0, 'true': 1}) # to adhere to pysam

    # merge fusions and breakends on 'fivePrimeBreakendId' & 'threePrimeBreakendId' and 'id' columns
    fusion_fivebreakend = fusions[['fivePrimeBreakendId', 'threePrimeBreakendId', 'FUSION_NAME', 'REPORTED']].merge(breakends[['id', 'svId']], left_on='fivePrimeBreakendId', right_on='id', how='left')
    fusion_threebreakend = fusions[['fivePrimeBreakendId', 'threePrimeBreakendId', 'FUSION_NAME', 'REPORTED']].merge(breakends[['id', 'svId']], left_on='threePrimeBreakendId', right_on='id', how='left')

    # merge the two
    fusion_breakend_merge = pd.concat([fusion_fivebreakend, fusion_threebreakend], ignore_index=True)


    # merge above with svs on 'svId' column in sv file
    fusion_breakend_svs_merge = fusion_breakend_merge.merge(svs, on='svId', how='left')

    # keep only the relevant columns for the final output - for debug - add columns svId, fivePrimeBreakendId, threePrimeBreakendId
    result = fusion_breakend_svs_merge[
        [
            "vcfId",
            "FUSION_NAME",
            "REPORTED"
        ]
    ]

    # remove duplicates (entries with same svId)
    result = result.drop_duplicates()

    # save the final merged dataframe to a tsv file
    result.to_csv(output_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description="Merge LINX files into one TSV file")

    parser.add_argument("-f", type=str, help="The fusion linx tsv file", required=True)
    parser.add_argument("-b", type=str, help="The breakends linx tsv file", required=True)
    parser.add_argument("-s", type=str, help="The svs linx tsv file", required=True)
    parser.add_argument("-o", type=str, help="name of output tsv file", required=True)

    args = parser.parse_args()

    # load tsv into pandas df
    fusions   = pd.read_csv(args.f, sep='\t', dtype=str)
    breakends = pd.read_csv(args.b, sep='\t', dtype=str)
    svs       = pd.read_csv(args.s, sep='\t', dtype=str)

    output_file = args.o

    merge_linx_files(fusions, breakends, svs, output_file)

if __name__ == "__main__":
    main()
