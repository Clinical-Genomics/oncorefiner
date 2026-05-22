import pandas as pd
import sys

# Load files
vcfentries = pd.read_csv(sys.argv[1], sep='\t', dtype=str, header=None, names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
linxinfo   = pd.read_csv(sys.argv[2], sep='\t', dtype=str)

# Merge vcfentries and linxinfo on 'ID' and 'vcfId' columns
vcf_linx_merge = vcfentries.merge(linxinfo, left_on='ID', right_on='vcfId', how='left')

result = vcf_linx_merge[
    [
        "#CHROM",
        "POS",
        "ID",
        "FUSIONID",
        "REPORTED_FUSION"
    ]
]

result.to_csv(sys.argv[3], sep='\t', index=False)
