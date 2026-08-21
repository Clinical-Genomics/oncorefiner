# Clinical-Genomics/oncorefiner

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/Clinical-Genomics/oncorefiner)
[![GitHub Actions CI Status](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/linting.yml/badge.svg)](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/Clinical-Genomics/oncorefiner)

## Introduction

**Clinical-Genomics/oncorefiner** is a customizable post-processing and extension layer for [nf-core/Oncoanalyser](https://github.com/nf-core/oncoanalyser) that adapts its outputs according to clinical and operational needs, adds missing analyses, and ensures flexibility for evolving standards while retaining Oncoanalyser's robust core.

<!-- TODO: Add information about the processes and the output when the architecture of the pipeline is decided and the information is clear.


### Workflow diagram

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->

1. Process CNV
   1. Prepare files from Oncoanalysers [`AMBER`](https://github.com/hartwigmedical/hmftools/tree/master/amber) and [`COBALT`](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) for visualization using GENS using cutom scripts `PREPARE_AMBER_FOR_GENS` and `PREPARE_COBALT_FOR_GENS` which is wrapped in the subworkflow `PREPARE_AMBER_COBALT_GENS`.

1. Process SNV VCF files
   1. Annotate with [`Vcfanno`](https://github.com/brentp/vcfanno). Intended for local/custom annotation.
   1. Filter with [`bcftools`](https://github.com/samtools/bcftools). This step aims to apply quality, population-level filtering and/or other general criteria as defined in the configuration settings.
   1. Annotate with [`Ensembl VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html)
   1. Rank variants and annotate with [`Genmod score`](https://github.com/Clinical-Genomics/genmod). A score is assigned to each variant based on the genmod score config file provided. Annotation with the score is added to the output vcf file.
   1. Filter with [`bcftools`](https://github.com/samtools/bcftools). This step applies clinically relevant filters as defined in the configuration settings. For example, it may involve subsetting variants based on a a list of clinically relevant genes.

1. Process SV VCF files
   1. Add annotation of fusions from LINX files using custom script `vcf_annotate_linx_fusions`.
   1. Annotate VCF with external database (params) using [`SVDB`](https://github.com/J35P312/SVDB).
   1. Filter with [`bcftools`](https://github.com/samtools/bcftools). This step aims to apply quality, population-level filtering and/or other general criteria as defined in the configuration settings.
   1. Annotate with [`Ensembl VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html).
   1. Filter with [`bcftools`](https://github.com/samtools/bcftools). This step applies clinically relevant filters as defined in the configuration settings. For example, it may involve subsetting variants based on a a list of clinically relevant genes.

1. Process CNV TSV files
   1. Generate interactive CNV html report using [`rmarkdownnotebook`](https://rmarkdown.rstudio.com/), and the script `assets/cnv_report.Rmd`, based on the result files from `oncoanalyser`: ` *.purple.cnv.gene.tsv` and `*.purple.cnv.somatic.tsv`.
1. Present QC for raw reads ([`MultiQC`](http://multiqc.info/)).

For further information about the each step and output files, please refer to the [output documentation](.github/docs/output.md).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Now, you can run the pipeline using:

```bash
nextflow run Clinical-Genomics/oncorefiner \
   -profile <docker/singularity/.../institute> \
   -params-file <params.yaml/params.json> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](./docs/usage.md) and the [parameter documentation](./docs/parameters.md).

## Pipeline output

For more details about the output files and reports, please refer to the [output documentation](.github/docs/output.md).

## Credits

Clinical-Genomics/oncorefiner was originally written by Clinical Genomics Stockholm.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Beatriz Sá Vinhas](https://github.com/beatrizsavinhas)
- [Eva Caceres](https://github.com/fevac)
- [Kristine Bilgrav Sæther](https://github.com/kristinebilgrav)
- [Felix Lenner](https://github.com/fellen31)
- [Mathias Johansson](https://github.com/mathiasbio)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use Clinical-Genomics/oncorefiner for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
