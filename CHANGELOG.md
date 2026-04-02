# Clinical-Genomics/oncorefiner: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of Clinical-Genomics/oncorefiner, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added Ensembl VEP annotation for SNV vcf file [#1](https://github.com/Clinical-Genomics/oncorefiner/pull/1)
- Added VCFANNO annotation for SNV vcf file [#2](https://github.com/Clinical-Genomics/oncorefiner/pull/2)
- Added filtering for SNV vcf file [#3](https://github.com/Clinical-Genomics/oncorefiner/pull/3)
- Added annotation for SV vcf file [#4](https://github.com/Clinical-Genomics/oncorefiner/pull/4)
- Added filtering for SV vcf file [#5](https://github.com/Clinical-Genomics/oncorefiner/pull/5)
- Added small test profile. The related test dataset have been added as a branch called oncorefiner under [Clinical-Genomics/test-datasets](https://github.com/Clinical-Genomics/test-datasets/tree/oncorefiner) [#8](https://github.com/Clinical-Genomics/oncorefiner/pull/8)
- Added CI checks for `Conventional PR title`, `Updated changelog` and `Add PR checklist comment` [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18)
- Added parameters documentation [#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25)
- Added pre-commit hook for automatic generation of parameters documentation [#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25)
- Added Nextflow strict syntax compatibility [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)
- Added bam and bai parameters for tumor and normal samples, and logic to create joint bam_bai channels [#61](https://github.com/Clinical-Genomics/oncorefiner/pull/61)
- Added `sex` parameter [#62](https://github.com/Clinical-Genomics/oncorefiner/pull/62)

### Changed

- Updated PR template, PR checklist, feature request template, bug report template and issue template chooser [#24](https://github.com/Clinical-Genomics/oncorefiner/pull/24)
- Updated nf-schema to 2.6.1 [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)
- Updated minimum Nextflow version to 25.10.0 [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)
- Added wgs-cancer-pipeline projects list in the issue templates [#37](https://github.com/Clinical-Genomics/oncorefiner/pull/37)
- Updated link to Contributing Guidelines in the PR checklist to point to the rendered version of the document in `dev` [#56](https://github.com/Clinical-Genomics/oncorefiner/pull/56)

### `Fixed`

- Removed snv_vcf_tbi and sv_vcf_tbi parameter. VCF indexes are now automatically detected [#9](https://github.com/Clinical-Genomics/oncorefiner/pull/9)
- Renamed pipeline from postprocessing to oncorefiner []()
- Fixed linting issues [#20](https://github.com/Clinical-Genomics/oncorefiner/pull/20)
- Fixed nf-test to run a functional default test, and generated a snapshot [#26](https://github.com/Clinical-Genomics/oncorefiner/pull/26)
- Added missing description to bug_report.yml [32](https://github.com/Clinical-Genomics/oncorefiner/pull/32)
- Updated template settings to set organisation to `Clinical-Genomics` and skip unused features `igenomes` and `fastqc` [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35)
- Refactored `genome` parameter to have default value 'GRCh38' and no longer refer to igenomes [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35)
- Updated linting config to fix linting issues and re-added/removed checks for files where nf-core file structure is no longer required [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35)
- Updated template for nf-core/tools version 3.5.2 to apply updated settings and changes missed in previous template update ([14](https://github.com/Clinical-Genomics/oncorefiner/pull/14)) [#35](https://github.com/Clinical-Genomics/oncorefiner/pull/35)
- Fixed prepare_references config that was defined but not used [36](https://github.com/Clinical-Genomics/oncorefiner/pull/36)
- Fixed bug and formatting in feature request template [#39](https://github.com/Clinical-Genomics/oncorefiner/pull/39)
- Fixed merge mistake introduced in [#25](https://github.com/Clinical-Genomics/oncorefiner/pull/25) [#41](https://github.com/Clinical-Genomics/oncorefiner/pull/41)
- Added necessary GITHUB_TOKEN permissions for action add_pr_checklist_comment [#42](https://github.com/Clinical-Genomics/oncorefiner/pull/42)
- Updated all modules and removed deprecated `ch_versions` to implement latest nf-core changes that use the `versions` topic channel to collect software versions [#34](https://github.com/Clinical-Genomics/oncorefiner/pull/34)
- Fixed settings for `add_pr_checklist_comment` to allow action to run on a PR originated from a fork [#45](https://github.com/Clinical-Genomics/oncorefiner/pull/45)
- Added `species` parameter to provide information for annotation which was previously hardcoded [#49](https://github.com/Clinical-Genomics/oncorefiner/pull/49)
- Added settings and moved ungrouped parameters to relevant groups [#50](https://github.com/Clinical-Genomics/oncorefiner/pull/50)
- Fixed bug in `MULTIQC` input channel that prevented the step from running [#54](https://github.com/Clinical-Genomics/oncorefiner/pull/54)
- Refactored subworkflow `PREPARE_REFERENCES` to include logic for untarring vep cache and be called in the main workflow, before `ONCOREFINER` [#57](https://github.com/Clinical-Genomics/oncorefiner/pull/57)
- Fixed so that parameters are only accessed in `main.nf` and provided to subsequent workflows as `val_*` [#58](https://github.com/Clinical-Genomics/oncorefiner/pull/58)
- Generalised description of `sex` parameter [#64](https://github.com/Clinical-Genomics/oncorefiner/pull/64)
- Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/` [#68](https://github.com/Clinical-Genomics/oncorefiner/pull/68)

### `Dependencies`

### `Deprecated`

### `Removed`

- Removed CI checks `awstest` and `awsfulltest` [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18)
- Removed unused parameter `custom_extra_files` [#51](https://github.com/Clinical-Genomics/oncorefiner/pull/51)
