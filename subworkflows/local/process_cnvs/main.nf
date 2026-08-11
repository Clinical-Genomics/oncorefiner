/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// LOCAL SUBWORKFLOWS
//

include { GENS } from '../../../subworkflows/local/gens/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_CNVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_CNVS {

    take:
    ch_amber_baf_tsv_gz        // channel: [optional]  [val(meta), path(tsv.gz)]
    ch_cobalt_ratio_pcf_normal // channel: [optional]  [val(meta), path(pcf)]
    ch_cobalt_ratio_pcf_tumor  // channel: [optional]  [val(meta), path(pcf)]
    val_analysis_type          // string:  [optional]  analysis type, e.g. "tumor_only" or "tumor_normal"

    main:

    GENS (
        ch_amber_baf_tsv_gz,
        ch_cobalt_ratio_pcf_normal,
        ch_cobalt_ratio_pcf_tumor,
        val_analysis_type
    )

    emit:
    gens_baf_normal_tsv = GENS.out.gens_baf_normal_tsv // channel: [val(meta), path(tsv.gz)]
    gens_baf_normal_tbi = GENS.out.gens_baf_normal_tbi // channel: [val(meta), path(tsv.gz.tbi)]
    gens_baf_tumor_tsv  = GENS.out.gens_baf_tumor_tsv  // channel: [val(meta), path(tsv.gz)]
    gens_baf_tumor_tbi  = GENS.out.gens_baf_tumor_tbi  // channel: [val(meta), path(tsv.gz.tbi)]
    gens_cov_bed        = GENS.out.gens_cov_bed        // channel: [val(meta), path(bed.gz)]
    gens_cov_bed_tbi    = GENS.out.gens_cov_bed_tbi    // channel: [val(meta), path(bed.gz.tbi)]
}
