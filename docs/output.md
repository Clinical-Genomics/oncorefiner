# Clinical-Genomics/oncorefiner: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

### PROCESS_CNVS

This process generates an interactive CNV report and files compatible with visualization in GENS.

<details markdown="1">
<summary>Output files</summary>

- `cnv/`
  - `<meta.id>.cnv_report.html`: interactive CNV report.
  - `<meta.id>_CNV_BAF_GENS.<tumor or normal>.baf.zoom.tsv.gz`: a gzipped tsv containing b-allele frequency (BAF) for visualization using GENS, produced from Oncoanalysers AMBER files.
  - `<meta.id>_CNV_BAF_GENS.<tumor or normal>.baf.zoom.tsv.gz.tbi`: index file for the gzipped tsv with BAF.
  - `<meta.id>_<meta.sample_type>_CNV_COV_GENS.bed.gz`: a gzipped bed file with genomic coverage levels for visualization using GENS, produced from Oncoanalysers COBALT files.
  - `<meta.id>_<meta.sample_type>_CNV_COV_GENS.bed.gz.tbi`: index file for the gzipped bed with coverage levels.

</details>

[`RMARKDOWNNOTEBOOK`](https://github.com/rstudio/rmarkdown) - used to generate a custom interactive CNV report.
`PREPARE_AMBER_FOR_GENS` - Custom script that takes output from Oncoanalysers [`AMBER`](https://github.com/hartwigmedical/hmftools/tree/master/amber) VCF and converts the B-allele frequency levels into a GENS compatible zoom level.
`PREPARE_COBALT_FOR_GENS` - Custom script that takes output from Oncoanalysers [`COBALT`](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) VCF and converts the coverage levels into GENS compatible levels.

### PROCESS_SNVS

This process annotates, filters and ranks single nucleotide variants.

<details markdown="1">
<summary>Output files</summary>

- `snv`
  - `<meta.id>_SNV_annotated_vcfanno.vcf.gz`: a gzipped VCF containing annotated SNVs.
  - `<meta.id>_SNV_annotated_vcfanno.vcf.gz.tbi`: an index file for the gzipped VCF.
  - `<meta.id>_SNV_annotated_vep.vcf.gz`: a gzipped VCF from step 4 with annotated and filtered variants.
  - `<meta.id>_SNV_annotated_vep.vcf.gz.tbi`: an index file for the gzipped VCF.
  - `<meta.id>_SNV_annotated_vep.vcf.gz_summary.html`: a html summary file produced by VEP.
  - `<meta.id>_SNV_clinical_filtered_bcftools.vcf.gz`: a gzipped VCF from step 6 with annotated, ranked and clinically filtered variants.
  - `<meta.id>_SNV_clinical_filtered_bcftools.vcf.gz.tbi`: an index file for the gzipped VCF.
  - `<meta.id>_SNV_genmod_score.vcf.gz`: a gzipped VCF from step 5 with annotated, filtered and ranked variants. Only produced if GENMOD config is provided. If not provided, the corresponding file would be the `_SNV_annotated_vep.vcf.gz`.
  - `<meta.id>_SNV_genmod_score.vcf.gz.tbi`: an index file for the gzipped VCF.

</details>

[`Vcfanno`](https://github.com/brentp/vcfanno) annotates VCF files with a number of INFO fields from the VCFs or BED files provided.
[`Genmod score`](https://github.com/Clinical-Genomics/genmod) assigns a score to each variant based on the genmod score config file provided. The output VCF file is annotated with the following INFO fields which reflect the assigned score:

```
RankScore
RankScoreNormalized
RankScoreMinMax
RankResult
```

[`bcftools`](https://github.com/samtools/bcftools) - This tool can filter VCFs using custom settings. In step 2 it applies quality and population level filtering, whilst in step 6 applies clinically relevant filters as defined in the configuration settings.
[`CADD`](https://github.com/kircherlab/CADD-scripts/) - This tool will annotate indels with a deleteriousness score based on the tools resources and calculations.
[`Ensembl VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html) - Annotation using the Ensembk VEP resources.

### PROCESS_SVS

This process annotates, ranks and filters structural variants

<details markdown="1">
<summary>Output files</summary>

- `sv/`
  - `<meta.id>_SV_research_filtered_bcftools.vcf.gz`: a gzipped VCF containing LINX fusion and SVDB annotated variants filtered with bcftools from step 3 in README.md.
  - `<meta.id>_SV_research_filtered_bcftools.vcf.gz.tbi`: an index file for the gzipped VCF.
  - `<meta.id>_SV_annotated_vep.vcf.gz`: a gzipped VCF containing the variants from the file above, annotated with VEP from step 4.
  - `<meta.id>_SV_annotated_vep.vcf.gz.tbi`: an index file for the above gzipped VCF.
  - `<meta.id>_SV_annotated_vep.vcf.gz_summary.html`: a html summary file produced by VEP.
  - `<meta.id>_SV_clinical_filtered_bcftools.vcf.gz`: a gzipped VCF containing variants annoated by VEP, ranked using GENMOD and filtered using bcftools from step 6.
  - `<meta.id>_SV_clinical_filtered_bcftools.vcf.gz.tbi`: an index file for the above gzipped VCF.
  - `<meta.id>_[tumor/normal]_SV_vcf2cytosure.cgh`: cgh file produced from step 7 to use for visualization in Cytosure.

</details>

`vcf_annotate_linx_fusions` - Annotation of fusions from [`LINX`](https://github.com/hartwigmedical/hmftools/tree/master/linx) to the VCF using cutom scripts. The VCF from LINX in Oncoanalyser displays the SVs as two instances, where only one entry will be annotated with the LINX information in Oncorefiner. The two instances can be connected manually using the `SVID` in the INFO field.  
[`SVDB`](https://github.com/J35P312/SVDB) - The tool applies annotation from external databases to the VCF.
[`bcftools`](https://github.com/samtools/bcftools) - This tool can filter VCFs using custom settings. In step 3 above it applies quality and population level filtering, whilst in step 6, applies clinically relevant filters as defined in configuration settings.
[`Ensembl VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html) - Annotation using the Ensembk VEP resources.
[`Genmod score`](https://github.com/Clinical-Genomics/genmod) assigns a score to each variant based on the genmod score config file provided. The output VCF file is annotated with the following INFO fields which reflect the assigned score:

```
RankScore
RankScoreNormalized
RankScoreMinMax
RankResult
```

[`vcf2cytosure`](https://github.com/NBISweden/vcf2cytosure) - Generation of cytosure files for visualization in Cytosure.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### BAM to CRAM

<details markdown="1">
<summary>Output files</summary>

- `alignments/`
  - `<meta.sample_id>_<meta.sample_type>.cram`: CRAM file for the input BAM
  - `<meta.sample_id>_<meta.sample_type>.cram.crai`: CRAI file for the input BAI

</details>

[`SAMTOOLS_VIEW`](https://www.htslib.org/doc/samtools-view.html) - This module compresses the input BAM file(s) into CRAM for storage.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
