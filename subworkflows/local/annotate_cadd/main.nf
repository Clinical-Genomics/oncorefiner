//
// A subworkflow to annotate cadd
//

include { BCFTOOLS_ANNOTATE as RENAME_CHR_CADD     } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as ANNOTATE_INDELS     } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW                            } from '../../../modules/nf-core/bcftools/view/main'
include { CADD                                     } from '../../../modules/nf-core/cadd/main'
include { GAWK as REFERENCE_TO_CADD_CHRNAMES       } from '../../../modules/nf-core/gawk/main'
include { GAWK as CADD_TO_REFERENCE_CHRNAMES       } from '../../../modules/nf-core/gawk/main'
include { TABIX_TABIX as TABIX_CADD                } from '../../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_CADD {

    take:
        ch_vcf                     // channel: [mandatory] [ val(meta), path(vcf), path(idx) ]
        val_genome                 // string:  [mandatory] GRCh37 or GRCh38
        ch_fai                     // channel: [mandatory] [ val(meta), path(fai) ]
        ch_header                  // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources          // channel: [mandatory] [ val(meta), path(dir) ]
        ch_cadd_prescored_indels   // channel: [mandatory] [ val(meta), path(dir) ]

    main:

        ch_rename_chrs_ref    = channel.value([[]])

        // Create files and rename chromosomes if reference is GRCh38
        if (val_genome.equals('GRCh38')) {

            // Create txt files for changing of chromosomes
            REFERENCE_TO_CADD_CHRNAMES ( ch_fai , [], false )

            REFERENCE_TO_CADD_CHRNAMES.out.output.map { _meta, txt -> txt }
                .set {ch_chrnames_cadd}

            CADD_TO_REFERENCE_CHRNAMES ( ch_fai , [], false )

            CADD_TO_REFERENCE_CHRNAMES.out.output.map { _meta, txt -> txt }
                .set { ch_rename_chrs_ref }

            ch_vcf
                .combine(ch_chrnames_cadd)
                .map { meta, vcf, tbi, txt -> tuple( meta, vcf, tbi, [], [], [], [], txt ) }
                .set {rename_chrnames_in}

            // Change chr names to CADD compatible names
            RENAME_CHR_CADD( rename_chrnames_in )

            RENAME_CHR_CADD.out.vcf
                .map {meta, vcf -> tuple( meta , vcf, [] )}
                .set { ch_vcf }
        }

        // Filter to extract indels
        BCFTOOLS_VIEW(ch_vcf, [], [], [])

        // CADD
        CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_resources, ch_cadd_prescored_indels)

        // Index CADD
        TABIX_CADD(CADD.out.tsv)

        // Change chr names back to desired naming and annotate original vcf with cadd results
        ch_vcf
            .join(CADD.out.tsv, failOnMismatch: true, failOnDuplicate: true)
            .join(TABIX_CADD.out.index, failOnMismatch: true, failOnDuplicate: true)
            .combine( ch_header )
            .combine( ch_rename_chrs_ref )
            .map { meta, vcf, tbi, annotations, annotations_index, header, txt -> tuple( meta, vcf, [], annotations, annotations_index, [], header, txt )  }
            .set { ch_annotate }


        ANNOTATE_INDELS( ch_annotate )

    emit:
        vcf = ANNOTATE_INDELS.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi = ANNOTATE_INDELS.out.tbi // channel: [ val(meta), path(tbi) ]
}
