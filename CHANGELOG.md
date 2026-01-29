# nf-core/postprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/postprocessing, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added Ensembl VEP annotation for SNV vcf file [#1](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/1)
- Added VCFANNO annotation for SNV vcf file [#2](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/2)
- Added filtering for SNV vcf file [#3](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/3)
- Added annotation for SV vcf file [#4](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/4)
- Added filtering for SV vcf file [#5](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/5)
- Added small test profile. The related test dataset have been added as a branch called oncorefiner under [Clinical-Genomics/test-datasets](https://github.com/Clinical-Genomics/test-datasets/tree/oncorefiner) [#8](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/8)

### `Fixed`

- Removed snv_vcf_tbi and sv_vcf_tbi parameter. VCF indexes are now automatically detected [#9](https://github.com/Clinical-Genomics/nf-core-postprocessing/pull/9)

### `Dependencies`

### `Deprecated`
