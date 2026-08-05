//
// A subworkflow to annotate indels with CADD scores
//

include { BCFTOOLS_ANNOTATE as BCFTOOLS_RENAME_CHR_CADD } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_INDELS } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW                                 } from '../../../modules/nf-core/bcftools/view/main'
include { CADD                                          } from '../../../modules/nf-core/cadd/main'
include { GAWK as GAWK_REF_TO_CADD_CHRNAMES             } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_CADD_TO_REF_CHRNAMES             } from '../../../modules/nf-core/gawk/main'
include { TABIX_TABIX as TABIX_CADD                     } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_INPUT                    } from '../../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_CADD {

    take:
    ch_vcf                   // channel: [mandatory] [val(meta), path(vcf)]
    val_genome               // string:  [mandatory] GRCh37 or GRCh38
    ch_fai                   // channel: [mandatory] [val(meta), path(fai)]
    ch_header                // channel: [mandatory] [val(meta), path(header)]
    ch_cadd_resources        // channel: [mandatory] [val(meta), path(dir)]
    ch_cadd_prescored_indels // channel: [mandatory] [val(meta), path(dir)]

    main:

    ch_rename_chrs_ref    = channel.value([[]])

    TABIX_INPUT(ch_vcf) //Subworkflow needs tabix index

    ch_vcf_tbi = ch_vcf
        .join(TABIX_INPUT.out.index, failOnMismatch:true, failOnDuplicate:true)

    // Create files and rename chromosomes if reference is GRCh38
    if (val_genome.equals('GRCh38')) {

        // Create txt files for changing of chromosomes
        GAWK_REF_TO_CADD_CHRNAMES ( ch_fai , [], false )

        GAWK_REF_TO_CADD_CHRNAMES.out.output.map { _meta, txt -> txt }
            .set {ch_chrnames_cadd}

        GAWK_CADD_TO_REF_CHRNAMES ( ch_fai , [], false )

        GAWK_CADD_TO_REF_CHRNAMES.out.output.map { _meta, txt -> txt }
            .set { ch_rename_chrs_ref }

        rename_chrnames_in = ch_vcf_tbi
            .combine(ch_chrnames_cadd)
            .map { meta, vcf, tbi, txt -> tuple( meta, vcf, tbi, [], [], [], [], txt ) }

        // Change chr names to CADD compatible names
        BCFTOOLS_RENAME_CHR_CADD( rename_chrnames_in )

        ch_vcf_tbi = BCFTOOLS_RENAME_CHR_CADD.out.vcf
            .map {meta, vcf -> tuple( meta , vcf, [] )}
    }

    // Filter to extract indels
    BCFTOOLS_VIEW(ch_vcf_tbi, [], [], [])

    // CADD
    CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_resources, ch_cadd_prescored_indels)

    // Index CADD
    TABIX_CADD(CADD.out.tsv)

    // Change chr names back to desired naming and annotate original vcf with cadd results
    ch_annotate = ch_vcf_tbi
        .join(CADD.out.tsv, failOnMismatch: true, failOnDuplicate: true)
        .join(TABIX_CADD.out.index, failOnMismatch: true, failOnDuplicate: true)
        .combine( ch_header )
        .combine( ch_rename_chrs_ref )
        .map { meta, vcf, tbi, annotations, annotations_index, header, txt -> tuple( meta, vcf, tbi, annotations, annotations_index, [], header, txt )  } //THERE IS A TBI?


    BCFTOOLS_ANNOTATE_INDELS( ch_annotate )

    emit:
    vcf = BCFTOOLS_ANNOTATE_INDELS.out.vcf // channel: [val(meta), path(vcf)]
    tbi = BCFTOOLS_ANNOTATE_INDELS.out.tbi // channel: [val(meta), path(tbi)]
}
