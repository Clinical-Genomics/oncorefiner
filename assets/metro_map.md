%%metro title: Clinical-Genomics/oncorefiner
%%metro style: dark
%%metro line: snv | PROCESS SNVs | #0570b0
%%metro line: sv | PROCESS SVs | #2db572
%%metro line: cnv | PROCESS CNVs | #f5c542
%%metro line: qc | QC | #e63946 | dashed
%%metro line: compression | File compression | #7b2d3b | dashed
%%metro legend: br

graph LR
    input[Oncoanalyser \n output]
    snv_vcfanno[vcfanno \n custom annotation]
    snv_vep[VEP]
    snv_cadd[CADD \n annotation]
    snv_research_bcftools[bcftools \n filter]
    snv_clinical_bcftools[bcftools \n clinical filter]
    snv_genmod[genmod ranking]
    sv_svdb[SVDB \n custom annotation]
    sv_vep[VEP]
    sv_research_bcftools[bcftools \n filter]
    sv_clinical_bcftools[bcftools \n clinical filter]
    sv_genmod[genmod ranking]
    sv_vcf2cytosure[vcf2cytosure]
    cnv_gens[GENS]
    cnv_report[CNV report]
    qc_report[MultiQC]
    bam2cram[BAM to CRAM]
    
    input -->|snv| snv_vcfanno
    snv_vcfanno -->|snv| snv_research_bcftools
    snv_research_bcftools -->|snv| snv_cadd
    snv_cadd -->|snv| snv_vep
    snv_vep -->|snv| snv_genmod
    snv_genmod -->|snv| snv_clinical_bcftools
    input -->|sv| sv_svdb
    sv_svdb -->|sv| sv_research_bcftools
    sv_research_bcftools -->|sv| sv_vep
    sv_vep -->|sv| sv_genmod
    sv_genmod -->|sv| sv_clinical_bcftools
    sv_clinical_bcftools -->|sv| sv_vcf2cytosure
    input -->|cnv| cnv_gens
    input -->|cnv| cnv_report
    input -->|qc| qc_report
    input -->|compression| bam2cram
