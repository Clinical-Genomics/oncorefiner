/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// LOCAL MODULES
//
include { PREPARE_AMBER_FOR_GENS  } from '../../../modules/local/prepare_amber_for_gens/main'
include { PREPARE_COBALT_FOR_GENS } from '../../../modules/local/prepare_cobalt_for_gens/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN GENS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENS {

    take:
    ch_amber_baf_tsv     // channel: [required]  [val(meta), path(tsv)]
    ch_cobalt_pcf_normal // channel: [required]  [val(meta), path(cobalt_pcf_normal)]
    ch_cobalt_pcf_tumor  // channel: [required]  [val(meta), path(cobalt_pcf_tumor)]
    val_analysis_type    // string:  [required]  analysis type, e.g. "tumor_only" or "tumor_normal"

    main:

    PREPARE_AMBER_FOR_GENS (
        ch_amber_baf_tsv,
        val_analysis_type
    )

    def ch_prepare_cobalt_for_gens_in = channel.empty().mix(
        ch_cobalt_pcf_normal,
        ch_cobalt_pcf_tumor
    )

    PREPARE_COBALT_FOR_GENS (
        ch_prepare_cobalt_for_gens_in
    )

    emit:
    gens_baf_normal_tsv = PREPARE_AMBER_FOR_GENS.out.normal_tsv // channel: [val(meta), path(tsv.gz)]
    gens_baf_normal_tbi = PREPARE_AMBER_FOR_GENS.out.normal_tbi // channel: [val(meta), path(tbi)]
    gens_baf_tumor_tsv  = PREPARE_AMBER_FOR_GENS.out.tumor_tsv  // channel: [val(meta), path(tsv.gz)]
    gens_baf_tumor_tbi  = PREPARE_AMBER_FOR_GENS.out.tumor_tbi  // channel: [val(meta), path(tbi)]
    gens_cov_bed        = PREPARE_COBALT_FOR_GENS.out.bed       // channel: [val(meta), path(bed.gz)]
    gens_cov_bed_tbi    = PREPARE_COBALT_FOR_GENS.out.tbi       // channel: [val(meta), path(tbi)]

}
