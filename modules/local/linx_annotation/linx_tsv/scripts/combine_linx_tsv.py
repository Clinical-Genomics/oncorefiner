import pandas as pd
import sys

"""
This script aim to first combine the three linx files (fusions, breakends, svs) into one merged tsv file
"""

# Load tsv files TODO: add click module
fusions   = pd.read_csv(sys.argv[1], sep='\t', dtype=str)
breakends = pd.read_csv(sys.argv[2], sep='\t', dtype=str)
svs       = pd.read_csv(sys.argv[3], sep='\t', dtype=str)

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
result.to_csv(sys.argv[4], sep='\t', index=False)
