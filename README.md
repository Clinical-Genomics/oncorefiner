# Clinical-Genomics/oncorefiner

<!-- TODO:
   Add logo here.
-->

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/Clinical-Genomics/oncorefiner)
[![GitHub Actions CI Status](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/linting.yml/badge.svg)](https://github.com/Clinical-Genomics/oncorefiner/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/Clinical-Genomics/oncorefiner)

## Introduction

**Clinical-Genomics/oncorefiner** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

### Workflow diagram

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

1. Process SNV VCF files
   1. Annotate with [`Vcfanno`](https://github.com/brentp/vcfanno)
   1. Filter according to call quality with [`bcftools`](https://github.com/samtools/bcftools)
   1. Filter according to user provided list of research relevant variant with [`bcftools`](https://github.com/samtools/bcftools)
   1. Annotate with [`VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html)
   1. Filter according to user provided list of clinically relevant variants with [`bcftools`](https://github.com/samtools/bcftools)

1. Process SV VCF files
   1. SVDB QUERY ???
   1. Filter according to call quality with ???
   1. Filter according to user provided list of research relevant variant with [`bcftools`](https://github.com/samtools/bcftools)
   1. Annotate with [`VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html)
   1. Filter according to user provided list of clinically relevant variants with [`bcftools`](https://github.com/samtools/bcftools)

1. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

### Summary of tools and version used in the pipeline

| Step                  | Tool          | Version |
| --------------------- | ------------- | ------- |
| Clinical Filtering    | bcftools      | 1.22    |
| Clinical Filtering SV | bcftools      | 1.22    |
| EnsemblVEP SNV        | ensemblvep    | 115.2   |
| EnsemblVEP SNV        | perl-math-cdf | 0.1     |
| EnsemblVEP SNV        | tabix         | 1.21    |
| EnsemblVEP SV         | ensemblvep    | 115.2   |
| EnsemblVEP SV         | perl-math-cdf | 0.1     |
| EnsemblVEP SV         | tabix         | 1.21    |
| Research Filtering    | bcftools      | 1.22    |
| Research Filtering SV | bcftools      | 1.22    |
| SVDB Query DB         | svdb          | 2.8.4   |
| Untar VEP Cache       | untar         | 1.34    |
| Vcfanno               | vcfanno       | 0.3.7   |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SAMPLE_ID,SAMPLE_ID_S1_L002_R1_001.fastq.gz,SAMPLE_ID_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end with only `fastq_1`) or a pair of fastq files (paired end with `fastq_1` and `fastq_2`).

Now, you can run the pipeline using:

```bash
nextflow run Clinical-Genomics/oncorefiner \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](./docs/usage.md) and the [parameter documentation](./docs/parameters.md).

## Pipeline output

For more details about the output files and reports, please refer to the [output documentation](.github/docs/output.md).

## Credits

Clinical-Genomics/oncorefiner was originally written by Clinical Genomics Stockholm.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Eva Caceres](https://github.com/fevac)
- [Kristine Bilgrav Sæther](https://github.com/kristinebilgrav)
- [Beatriz Sá Vinhas](https://github.com/beatrizsavinhas)
- [Mathias Johansson](https://github.com/mathiasbio)
- [Felix Lenner](https://github.com/fellen31)

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
