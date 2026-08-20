%%metro title: Clinical-Genomics/oncorefiner
%%metro style: dark
%%metro line: snv | PROCESS SNVs | #0570b0
%%metro line: sv | PROCESS SVs | #2db572
%%metro line: cnv | PROCESS CNVs | #f5c542
%%metro line: qc | QC | #e63946 | dashed
%%metro line: compression | File compression | #7b2d3b | dashed
%%metro legend: br

graph LR
    input[Oncoanalyser Output]
    snv_custom_annot[Custom annotation]
    snv_ensembl_vep[Ensemble VEP]
    cadd[CADD]
    snv_research_filter[Research Filter]
    snv_clinical_filter[Clinical Filter]
    snv_ranking[Variant ranking]
    sv_custom_annot[Custom annotation]
    sv_ensembl_vep[Ensemble VEP]
    sv_research_filter[Research Filter]
    sv_clinical_filter[Clinical Filter]
    sv_clinical_filter[Clinical Filter]
    sv_ranking[Variant ranking]
    cytosure_support[Cytosure Support]
    gens[GENS support]
    cnv_report[CNV report]
    qc_report[MultiQC]
    qc_thresholds[QC thresholds]
    bam2cram[BAM to CRAM]

    input -->|snv| snv_custom_annot
    snv_custom_annot -->|snv| snv_ensembl_vep
    snv_ensembl_vep -->|snv| cadd
    cadd -->|snv| snv_research_filter
    snv_research_filter -->|snv| snv_clinical_filter
    snv_clinical_filter -->|snv| snv_ranking
    input -->|sv| sv_custom_annot
    sv_custom_annot -->|sv| sv_ensembl_vep
    sv_ensembl_vep -->|sv| sv_research_filter
    sv_research_filter -->|sv| sv_clinical_filter
    sv_clinical_filter -->|sv| sv_ranking
    sv_ranking -->|sv| cytosure_support
    input -->|cnv| gens
    input -->|cnv| cnv_report
    input -->|qc| qc_report
    input -->|qc| qc_thresholds
    input -->|compression| bam2cram
