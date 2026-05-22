import pandas as pd
import sys

# Load tsv files
fusions   = pd.read_csv(sys.argv[1], sep='\t', dtype=str)
breakends = pd.read_csv(sys.argv[2], sep='\t', dtype=str)
svs       = pd.read_csv(sys.argv[3], sep='\t', dtype=str)

fusions.rename(columns={'name': 'FUSIONID', 'reported': 'REPORTED_FUSION'}, inplace=True)

# merge fusions and breakends on 'fivePrimeBreakendId' ('threePrimeBreakendId') and 'id' columns
fusion_breakend_merge = fusions[['fivePrimeBreakendId', 'threePrimeBreakendId', 'FUSIONID', 'REPORTED_FUSION']].merge(breakends[['id', 'svId']], left_on='fivePrimeBreakendId', right_on='id', how='left')

# merge above with svs on 'svId' column
fusion_breakend_svs_merge = fusion_breakend_merge.merge(svs, on='svId', how='left')

# keep only the relevant columns for the final output
result = fusion_breakend_svs_merge[
    [
        "vcfId",
        "svId",
        "fivePrimeBreakendId",
        "threePrimeBreakendId",
        "FUSIONID",
        "REPORTED_FUSION"
    ]
]

# save the final merged dataframe to a tsv file
result.to_csv(sys.argv[4], sep='\t', index=False)
