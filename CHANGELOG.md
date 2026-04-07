# Clinical-Genomics/oncorefiner: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of Clinical-Genomics/oncorefiner, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#1](https://github.com/Clinical-Genomics/oncorefiner/pull/1) Added Ensembl VEP annotation for SNV vcf file.
- [#2](https://github.com/Clinical-Genomics/oncorefiner/pull/2) Added VCFANNO annotation for SNV vcf file.
- [#3](https://github.com/Clinical-Genomics/oncorefiner/pull/3) Added filtering for SNV vcf file.
- [#4](https://github.com/Clinical-Genomics/oncorefiner/pull/4) Added annotation for SV vcf file.
- [#5](https://github.com/Clinical-Genomics/oncorefiner/pull/5) Added filtering for SV vcf file.
- [#8](https://github.com/Clinical-Genomics/oncorefiner/pull/8) Added small test profile. The related test dataset have been added as a branch called oncorefiner under [Clinical-Genomics/test-datasets](https://github.com/Clinical-Genomics/test-datasets/tree/oncorefiner).
- [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18) Added CI checks for `Conventional PR title`, `Updated changelog` and `Add PR checklist comment`.
- [#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25) Added parameters documentation.
- [#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25) Added pre-commit hook for automatic generation of parameters documentation.
- [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30) Added Nextflow strict syntax compatibility.
- [#61](https://github.com/Clinical-Genomics/oncorefiner/pull/61) Added bam and bai parameters for tumor and normal samples, and logic to create joint bam_bai channels.
- [#62](https://github.com/Clinical-Genomics/oncorefiner/pull/62) Added `sex` parameter.

### `Changed`

- [#24](https://github.com/Clinical-Genomics/oncorefiner/pull/24) Updated PR template, PR checklist, feature request template, bug report template and issue template chooser.
- [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30) Updated nf-schema to 2.6.1.
- [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30) Updated minimum Nextflow version to 25.10.0.
- [#37](https://github.com/Clinical-Genomics/oncorefiner/pull/37) Added wgs-cancer-pipeline projects list in the issue templates.
- [#56](https://github.com/Clinical-Genomics/oncorefiner/pull/56) Updated link to Contributing Guidelines in the PR checklist to point to the rendered version of the document in `dev`.

### `Fixed`

- [#9](https://github.com/Clinical-Genomics/oncorefiner/pull/9) Removed snv_vcf_tbi and sv_vcf_tbi parameter. VCF indexes are now automatically detected.
- [#10](https://github.com/Clinical-Genomics/oncorefiner/pull/10) Renamed pipeline from postprocessing to oncorefiner.
- [#20](https://github.com/Clinical-Genomics/oncorefiner/pull/20) Fixed linting issues.
- [#26](https://github.com/Clinical-Genomics/oncorefiner/pull/26) Fixed nf-test to run a functional default test, and generated a snapshot.
- [#32](https://github.com/Clinical-Genomics/oncorefiner/pull/32) Added missing description in bug report template.
- [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35) Updated template settings to set organisation to `Clinical-Genomics` and skip unused features `igenomes` and `fastqc`.
- [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35) Refactored `genome` parameter to have default value 'GRCh38' and no longer refer to igenomes.
- [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35) Updated linting config to fix linting issues and re-added/removed checks for files where nf-core file structure is no longer required.
- [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35) Updated template for nf-core/tools version 3.5.2 to apply updated settings and changes missed in previous template update ([#14](https://github.com/Clinical-Genomics/oncorefiner/pull/14)).
- [36](https://github.com/Clinical-Genomics/oncorefiner/pull/36) Fixed prepare_references config that was defined but not used.
- [#39](https://github.com/Clinical-Genomics/oncorefiner/pull/39) Fixed bug and formatting in feature request template.
- [#41](https://github.com/Clinical-Genomics/oncorefiner/pull/41) Fixed merge mistake in `.nf-core.yml` introduced in previous PR ([#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25)).
- [#42](https://github.com/Clinical-Genomics/oncorefiner/pull/42) Added necessary GITHUB_TOKEN permissions for action add_pr_checklist_comment.
- [#34](https://github.com/Clinical-Genomics/oncorefiner/pull/34) Updated all modules and removed deprecated `ch_versions` to implement latest nf-core changes that use the `versions` topic channel to collect software versions.
- [#45](https://github.com/Clinical-Genomics/oncorefiner/pull/45) Fixed settings for `add_pr_checklist_comment` to allow action to run on a PR originated from a fork.
- [#49](https://github.com/Clinical-Genomics/oncorefiner/pull/49) Added `species` parameter to provide information for annotation which was previously hardcoded.
- [#50](https://github.com/Clinical-Genomics/oncorefiner/pull/50) Added settings and moved ungrouped parameters to relevant groups.
- [#54](https://github.com/Clinical-Genomics/oncorefiner/pull/54) Fixed bug in `MULTIQC` input channel that prevented the step from running.
- [#57](https://github.com/Clinical-Genomics/oncorefiner/pull/57) Refactored subworkflow `PREPARE_REFERENCES` to include logic for untarring vep cache and be called in the main workflow, before `ONCOREFINER`.
- [#58](https://github.com/Clinical-Genomics/oncorefiner/pull/58) Fixed so that parameters are only accessed in `main.nf` and provided to subsequent workflows as `val_*`.
- [#64](https://github.com/Clinical-Genomics/oncorefiner/pull/64) Generalised description of `sex` parameter.
- [#68](https://github.com/Clinical-Genomics/oncorefiner/pull/68) Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/`.
- [#48](https://github.com/Clinical-Genomics/oncorefiner/pull/48) Updated documentation.

### `Dependencies`

### `Deprecated`

### `Removed`

- [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18) Removed CI checks `awstest` and `awsfulltest`.
- [#51](https://github.com/Clinical-Genomics/oncorefiner/pull/51) Removed unused parameter `custom_extra_files`.
