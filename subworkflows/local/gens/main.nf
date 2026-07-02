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
    ch_amber_baf_tsv     // channel: [required]  [val(meta), path(amber_baf_tsv)]
    ch_cobalt_pcf_normal // channel: [required]  [val(meta), path(cobalt_pcf_normal)]
    ch_cobalt_pcf_tumor  // channel: [required]  [val(meta), path(cobalt_pcf_tumor)]
    val_sampletype       // string:  [required]  sample type, e.g. "tumor_only" or "tumor_normal"

    main:
    ch_amber_baf_tsv_sampletype = ch_amber_baf_tsv
        .combine(val_sampletype)

    PREPARE_AMBER_FOR_GENS (
        ch_amber_baf_tsv_sampletype
    )

    ch_cobalt_pcf_tumor_normal = channel.empty().mix(
        ch_cobalt_pcf_normal,
        ch_cobalt_pcf_tumor
    )
    ch_cobalt_pcf_in = ch_cobalt_pcf_tumor_normal
        .multiMap { meta, cobalt_pcf ->
            pcf: tuple(meta, cobalt_pcf)
        }

    PREPARE_COBALT_FOR_GENS (
        ch_cobalt_pcf_in.pcf
    )

    emit:
    gens_baf_tsv_normal = PREPARE_AMBER_FOR_GENS.out.tsv_normal
    gens_baf_tbi_normal = PREPARE_AMBER_FOR_GENS.out.tbi_normal
    gens_baf_tsv_tumor  = PREPARE_AMBER_FOR_GENS.out.tsv_tumor
    gens_baf_tbi_tumor  = PREPARE_AMBER_FOR_GENS.out.tbi_tumor
    gens_cov_bed        = PREPARE_COBALT_FOR_GENS.out.bed
    gens_cov_bed_tbi    = PREPARE_COBALT_FOR_GENS.out.bed_tbi

}
