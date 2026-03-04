# nf-core/oncorefiner: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/oncorefiner, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added Ensembl VEP annotation for SNV vcf file [#1](https://github.com/Clinical-Genomics/oncorefiner/pull/1)
- Added VCFANNO annotation for SNV vcf file [#2](https://github.com/Clinical-Genomics/oncorefiner/pull/2)
- Added filtering for SNV vcf file [#3](https://github.com/Clinical-Genomics/oncorefiner/pull/3)
- Added annotation for SV vcf file [#4](https://github.com/Clinical-Genomics/oncorefiner/pull/4)
- Added filtering for SV vcf file [#5](https://github.com/Clinical-Genomics/oncorefiner/pull/5)
- Added small test profile. The related test dataset have been added as a branch called oncorefiner under [Clinical-Genomics/test-datasets](https://github.com/Clinical-Genomics/test-datasets/tree/oncorefiner) [#8](https://github.com/Clinical-Genomics/oncorefiner/pull/8)
- Added CI checks for `Conventional PR title`, `Updated changelog` and `Add PR checklist comment` [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18)
- Added Nextflow strict syntax compatibility [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)
- Updated nf-schema to 2.6.1 [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)

### `Fixed`

- Removed snv_vcf_tbi and sv_vcf_tbi parameter. VCF indexes are now automatically detected [#9](https://github.com/Clinical-Genomics/oncorefiner/pull/9)
- Renamed pipeline from postprocessing to oncorefiner []()
- Fixed linting issues [#20](https://github.com/Clinical-Genomics/oncorefiner/pull/20)
- Fixed minimum Nextflow version mismatches between GitHub CI tests and `nextflow.config` [#30](https://github.com/Clinical-Genomics/oncorefiner/pull/30)

### `Dependencies`

### `Deprecated`

### `Removed`

- Removed CI checks `awstest` and `awsfulltest` [#18](https://github.com/Clinical-Genomics/oncorefiner/pull/18)
