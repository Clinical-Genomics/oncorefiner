//
// A subworkflow to annotate cadd
//

include { BCFTOOLS_ANNOTATE                        } from '../../../modules/nf-core/bcftools/annotate/main'
include { CADD                                     } from '../../../modules/nf-core/cadd/main'
include { TABIX_TABIX as TABIX_CADD                } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_ANNOTATE            } from '../../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_CADD {

    take:
        ch_snv_vcf         // channel: [mandatory] [ val(meta), path(vcfs), path(idx) ]
        ch_cadd_header     // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources  // channel: [mandatory] [ path(dir) ]
        ch_cadd_prescored_indels // channel: [mandatory] [ val(meta), path(dir) ]

    main:
        ch_versions = channel.empty()

        CADD(ch_snv_vcf, ch_cadd_resources, ch_cadd_prescored_indels)

        TABIX_CADD(CADD.out.tsv)

        ch_snv_vcf
            .join(CADD.out.tsv)
            .join(TABIX_CADD.out.tbi)
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in, ch_cadd_header )

        TABIX_ANNOTATE (BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(CADD.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_CADD.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_ANNOTATE.out.versions.first())

    emit:
        vcf = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi = TABIX_ANNOTATE.out.tbi
        versions = ch_versions
}
