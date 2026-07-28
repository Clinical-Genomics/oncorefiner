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
- [#60](https://github.com/Clinical-Genomics/oncorefiner/pull/60) Added `GENERATE_CYTOSURE_FILES` subworkflow and necessary nf-core modules `TIDDIT_COV` and `VCF2CYTOSURE`.
- [#70](https://github.com/Clinical-Genomics/oncorefiner/pull/70) Added `SAMTOOLS/VIEW` for bam to cram conversion in the `main.nf`.
- [#66](https://github.com/Clinical-Genomics/oncorefiner/pull/66) Added `PROCESS_SNVS` subworkflow.
- [#73](https://github.com/Clinical-Genomics/oncorefiner/pull/73) Added `PROCESS_SVS` subworkflow.
- [#59](https://github.com/Clinical-Genomics/oncorefiner/pull/59) Added `ANNOTATE_CADD` subworkflow with following test (stub only), for CADD scoring of InDels, used in `PROCESS_SNVS`.
- [#69](https://github.com/Clinical-Genomics/oncorefiner/pull/69) Added `tumor_normal` config file, used by the default test profile.
- [#69](https://github.com/Clinical-Genomics/oncorefiner/pull/69) Added `tumor_only` config file, profile and pipeline test and snapshot.
- [#101](https://github.com/Clinical-Genomics/oncorefiner/pull/101) Added parameters for arguments in processes `PROCESS_SNVS:BCFTOOLS_VIEW_*`, `PROCESS_*:ENSEMBLVEP_VEP` and `ANNOTATE_CADD:BCFTOOLS_ANNOTATE_INDELS`.
- [#111](https://github.com/Clinical-Genomics/oncorefiner/pull/111) Added hidden `genome_version_number` parameter, parsed by default from `params.genome`.
- [#112](https://github.com/Clinical-Genomics/oncorefiner/pull/112) Added parameters `amber_baf_tsv_gz`, `cobalt_ratio_pcf_tumor`, `cobalt_ratio_pcf_normal` and `sampletype` needed for GENS workflow and populated them in the test dataset.
- [#117](https://github.com/Clinical-Genomics/oncorefiner/pull/117) Added stub version for all pipeline and local subworkflow tests.
- [#133](https://github.com/Clinical-Genomics/oncorefiner/pull/133) Added pipeline stub test without `vep_plugin_files` parameter.
- [#137](https://github.com/Clinical-Genomics/oncorefiner/pull/137) Added nf-schema validation to VEP database inputs (`vep_plugin_files`) and support for TSV, JSON and YAML file formats in addition to CSV.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/103) Added `genmod_score_config` parameter.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/103) Added `genomic-medicine-sweden/vcf_annotate_score_genmod` subworkflow to `PROCESS_SNVS`.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/103) Added `tests/optional_inputs_stub.nf.test` to test running the pipeline without optional inputs. Includes stub test for running the pipeline without providing the `genmod_score_config` parameter.
- [#135](https://github.com/Clinical-Genomics/oncorefiner/pull/135) Added test config for `GENERATE_CYTOSURE_FILES`.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Added metadata parameters `case_id`, `sample_id_tumor` and `sample_id_normal`.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Added `channelFromMetaAndPath` custom channel factory.
- [#136](https://github.com/Clinical-Genomics/oncorefiner/pull/136) Added nf-schema validation to SVDB database inputs (`svdb_query_dbs`)
- []() Added the local subworkflow `GENS`
- [#120](https://github.com/Clinical-Genomics/oncorefiner/pull/120) Added local module `PREPARE_AMBER_FOR_GENS`, using script in `/bin`, along with tests.
- [#122](https://github.com/Clinical-Genomics/oncorefiner/pull/122) Added local module `prepare_cobalt_for_gens` using python script in `/bin`, along with tests.

### `Changed`

- [#24](https://github.com/Clinical-Genomics/oncorefiner/pull/24) Updated PR template, PR checklist, feature request template, bug report template and issue template chooser.
- [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30) Updated nf-schema to 2.6.1.
- [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30) Updated minimum Nextflow version to 25.10.0.
- [#37](https://github.com/Clinical-Genomics/oncorefiner/pull/37) Added wgs-cancer-pipeline projects list in the issue templates.
- [#56](https://github.com/Clinical-Genomics/oncorefiner/pull/56) Updated link to Contributing Guidelines in the PR checklist to point to the rendered version of the document in `dev`.
- [#67](https://github.com/Clinical-Genomics/oncorefiner/pull/67) Updated nf-core subworkflow `utils_nfschema_plugin`
- [#71](https://github.com/Clinical-Genomics/oncorefiner/pull/71) Updated `pipelines_testdata_base_path` and paths to test vcf files.
- [#66](https://github.com/Clinical-Genomics/oncorefiner/pull/66) Moved logic for processing SNV VCF files, previously in `workflows/oncorefiner.nf`, to `PROCESS_SNVS` subworkflow.
- [#73](https://github.com/Clinical-Genomics/oncorefiner/pull/73) Moved logic for processing SV VCF files, previously in `workflows/oncorefiner.nf`, to `PROCESS_SVS` subworkflow.
- [#69](https://github.com/Clinical-Genomics/oncorefiner/pull/69) Renamed and refactored test config to `test_base.config` to include common parameters and files used for all tests.
- [#69](https://github.com/Clinical-Genomics/oncorefiner/pull/69) Refactored default pipeline test and test profile to run a `tumor_normal` default test and updated snapshot.
- [#87](https://github.com/Clinical-Genomics/oncorefiner/pull/87) Updated all testdata file paths to GRCh38 files, updated snapshot.
- [#100](https://github.com/Clinical-Genomics/oncorefiner/pull/100) Template update for nf-core tools 4.0.2.
- [#100](https://github.com/Clinical-Genomics/oncorefiner/pull/100) Changed minimum required nextflow version to 25.10.4.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/108) Updated `bcftools/view` module, also included in `genomic-medicine-sweden/vcf_annotate_score_genmod` subworkflow.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/103) Changed `PROCESS_SNVS` output logic to return `ch_research_filtered_vcf/tbi` after scoring if this step is run. The file names for these outputs are then given by `VCF_ANNOTATE_SCORE_GENMOD:BCFTOOLS_VIEW`, i.e. with prefix `${meta.id}_genmod_score` as specified in the config file.
- [#106](https://github.com/Clinical-Genomics/oncorefiner/pull/106) Changed logic for `stable_path` in pipeline tests to exclude all vcf and index files by default. Removed individual entries for these files in `tests/.nftignore`.
- [#108](https://github.com/Clinical-Genomics/oncorefiner/pull/108) Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/` that includes the genmod score config file.
- [#81](https://github.com/Clinical-Genomics/oncorefiner/pull/81) Update publishing strategy for `PROCESS_SNVS` subworkflow to use output blocks. Removed publishing settings from the subworkflow config file.
- [#88](https://github.com/Clinical-Genomics/oncorefiner/pull/88) Update contributing guidelines to follow current conventions on running `prek` for pre-commit hooks, indentation and item order in code blocks and remove outdated channel publishing instructions.
- [#127](https://github.com/Clinical-Genomics/oncorefiner/pull/127) Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/oncorefiner` that removes `testdata/samplesheet_test.csv`.
- [#128](https://github.com/Clinical-Genomics/oncorefiner/pull/128) Upgraded nf-test version to 0.9.5 so that the CI check supports topic channel tests.
- [#138](https://github.com/Clinical-Genomics/oncorefiner/pull/138) Changed indentation to consistently match the nf-core template.
- [#139](https://github.com/Clinical-Genomics/oncorefiner/pull/139) Changed `tests/.nftignore` to not include `process_svs/*.{vcf}`, since uncompressed VCF files are no longer output
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/108) Updated `bcftools/view` module, also included in `genomic-medicine-sweden/vcf_annotate_score_genmod` subworkflow.
- [#103](https://github.com/Clinical-Genomics/oncorefiner/pull/103) Changed `PROCESS_SNVS` output logic to return `ch_research_filtered_vcf/tbi` after scoring if this step is run. The file names for these outputs are then given by `VCF_ANNOTATE_SCORE_GENMOD:BCFTOOLS_VIEW`, i.e. with prefix `${meta.id}_genmod_score` as specified in the config file.
- [#142](https://github.com/Clinical-Genomics/oncorefiner/pull/142) Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/oncorefiner` that adds LINX TSV files.
- [#141](https://github.com/Clinical-Genomics/oncorefiner/pull/141) Refactored `GENERATE_CYTOSURE_FILES` to emit output from `VCF2CYTOSURE`.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Refactored input channel initialisation to update meta according to new metadata parameters for: Case specific file channels - `ch_snv_vcf`, `ch_snv_vcf_tbi`, `ch_sv_vcf`, `ch_sv_vcf_tbi`; and Sample specific file channel - `ch_bam_bai_normal`, `ch_bam_bai_tumor`, `ch_genome_fasta`, `ch_genome_fai`.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Update configuration settings to use `meta.sample_type` instead of the previous `meta.type`.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Update `generate_cytosure_files` and `process_svs` test and snapshots to accommodate to changes in pipeline configuration.
- [#124](https://github.com/Clinical-Genomics/oncorefiner/pull/124) Update documentation in `docs/usage.md` and `README.md` on how to start a pipeline run given the necessary metadata input parameters.
- [#145](https://github.com/Clinical-Genomics/oncorefiner/pull/145) Small fixes to `GENERATE_CYTOSURE_FILES` tests from #135: updated config regex to include the fully qualified name of the process in order for the settings to be applied with the highest priority and removed parameter dependency from each test by moving the values to the config file instead.
- [#149](https://github.com/Clinical-Genomics/oncorefiner/pull/149) Updated test data base path to remove unused column in SVDB files

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
- [#36](https://github.com/Clinical-Genomics/oncorefiner/pull/36) Fixed prepare_references config that was defined but not used.
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
- [#67](https://github.com/Clinical-Genomics/oncorefiner/pull/67) `--help` parameter not working
- [#72](https://github.com/Clinical-Genomics/oncorefiner/pull/72) Added subworkflow test to `PREPARE_REFERENCES`
- [#82](https://github.com/Clinical-Genomics/oncorefiner/pull/82) Add index files for SNV clinical and research filtered vcfs
- [#90](https://github.com/Clinical-Genomics/oncorefiner/pull/90) Fixed bug in `GENERATE_CYTOSURE` and `PROCESS_SNVs`: updated test meta to be `subject_a` and input to contain tbi.
- [#95](https://github.com/Clinical-Genomics/oncorefiner/pull/95) Fixed so VEP annotates in the same order for tests in `PROCESS_SVs`.
- [#96](https://github.com/Clinical-Genomics/oncorefiner/pull/96) Moved custom test settings for `PROCESS_SNVS:BCFTOOLS_VIEW_*` to test configs.
- [#101](https://github.com/Clinical-Genomics/oncorefiner/pull/101) Refactored config settings for `PROCESS_SNVS:BCFTOOLS_VIEW_*`, `PROCESS_*:ENSEMBLVEP_VEP`, and `ANNOTATE_CADD:BCFTOOLS_ANNOTATE_INDELS` to take parameters and specified specific test settings for these in `test_base.config`
- [#110](https://github.com/Clinical-Genomics/oncorefiner/pull/110) Added `genome` argument for `VCF2CYTOSURE`, given by `params.genome_version_number`, to ensure the reference genome used for the process matches the given input data.
- [#116](https://github.com/Clinical-Genomics/oncorefiner/pull/116) Updated `GENERATE_CYTOSURE_FILES` subworkflow test snapshot in sequence of [#110](https://github.com/Clinical-Genomics/oncorefiner/pull/110).
- [#123](https://github.com/Clinical-Genomics/oncorefiner/pull/123) Changed the subworkflow test when untarring is needed in `PREPARE_REFERENCES` to check `workflow.out.vep_resources` instead of `params.outdir`
- [#129](https://github.com/Clinical-Genomics/oncorefiner/pull/129) Changed the extra arguments `extra_args_*` to default to '' instead of null and removed empty `extra_args_snv_vep` in `test_base.config`
- [#130](https://github.com/Clinical-Genomics/oncorefiner/pull/130) Changed output logic for `PROCESS_SNVS` to output the `research_filtered_vcf` file after annotation with vep.
- [#131](https://github.com/Clinical-Genomics/oncorefiner/pull/131) Fixed strict Nextflow syntax compatibility in `ANNOTATE_CADD`
- [#137](https://github.com/Clinical-Genomics/oncorefiner/pull/134) Fixed VEP able to run without `vep_plugin_files` parameter.
- [#135](https://github.com/Clinical-Genomics/oncorefiner/pull/135) Fixed test for `GENERATE_CYTOSURE_FILES` so that it states the necessary parameters for the test, uses configurations from the test config and is independent from pipeline wide configuration.
- [#156](https://github.com/Clinical-Genomics/oncorefiner/pull/156) Update `pipelines_testdata_base_path` to reflect the latest commit to `Clinical-Genomics/test-datasets/` with addions of AMBER and COBALT files needed for GENS workflow.

### `Dependencies`

### `Deprecated`

### `Removed`

- [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18) Removed CI checks `awstest` and `awsfulltest`.
- [#51](https://github.com/Clinical-Genomics/oncorefiner/pull/51) Removed unused parameter `custom_extra_files`.
- [#125](https://github.com/Clinical-Genomics/oncorefiner/pull/125) Removed `input` parameter and samplesheet logic and documentation.
- [#139](https://github.com/Clinical-Genomics/oncorefiner/pull/139) Removed SVDB output file, since it's an intermediary file not useful to have in the outputs.
- [#152](https://github.com/Clinical-Genomics/oncorefiner/pull/152) Removed publishDir directive for `PREPARE_REFERENCES`.
