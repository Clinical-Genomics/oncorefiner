/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// LOCAL SUBWORKFLOWS
//

include { PREPARE_AMBER_COBALT_GENS } from '../../../subworkflows/local/prepare_amber_cobalt_gens/main'

//
// MODULE: Installed directly from nf-core/modules
//
include { RMARKDOWNNOTEBOOK } from '../../../modules/nf-core/rmarkdownnotebook/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_CNVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_CNVS {

    take:
    ch_amber_baf_tsv_gz        // channel: [optional]  [val(meta), path(tsv.gz)]
    ch_cnv_gene_tsv            // channel: [optional]  [val(meta), path(tsv)]
    ch_cnv_segment_tsv         // channel: [optional]  [val(meta), path(tsv)]
    ch_cobalt_ratio_pcf_normal // channel: [optional]  [val(meta), path(pcf)]
    ch_cobalt_ratio_pcf_tumor  // channel: [optional]  [val(meta), path(pcf)]
    val_analysis_type          // string:  [optional]  analysis type, e.g. "tumor_only" or "tumor_normal"

    main:

    // Run the PREPARE_AMBER_COBALT_GENS subworkflow
    PREPARE_AMBER_COBALT_GENS (
        ch_amber_baf_tsv_gz,
        ch_cobalt_ratio_pcf_normal,
        ch_cobalt_ratio_pcf_tumor,
        val_analysis_type
    )

    // Path to report template (Rmd file)
    // Todo: eventually this should be passed as an input to the workflow, but for now we can hardcode it since it's part of our assets
    def cnv_report_notebook = file("${projectDir}/assets/cnv_report.Rmd", checkIfExists: true)

    // Join the two distinct file channels together based on the meta.id
    ch_rmarkdownnotebook_in = ch_cnv_gene_tsv
        .join(ch_cnv_segment_tsv)
        .multiMap { meta, cnv_gene_tsv, cnv_segment_tsv ->
        notebook   : [meta, cnv_report_notebook]
        parameters : [cnv_gene: cnv_gene_tsv.name, cnv_segment: cnv_segment_tsv.name]
        input_files: [cnv_gene_tsv, cnv_segment_tsv]
    }


    RMARKDOWNNOTEBOOK(
    ch_rmarkdownnotebook_in.notebook,
    ch_rmarkdownnotebook_in.parameters,
    ch_rmarkdownnotebook_in.input_files,
    )

    emit:
    gens_baf_normal_tsv = PREPARE_AMBER_COBALT_GENS.out.gens_baf_normal_tsv // channel: [val(meta), path(tsv.gz)]
    gens_baf_normal_tbi = PREPARE_AMBER_COBALT_GENS.out.gens_baf_normal_tbi // channel: [val(meta), path(tsv.gz.tbi)]
    gens_baf_tumor_tsv  = PREPARE_AMBER_COBALT_GENS.out.gens_baf_tumor_tsv  // channel: [val(meta), path(tsv.gz)]
    gens_baf_tumor_tbi  = PREPARE_AMBER_COBALT_GENS.out.gens_baf_tumor_tbi  // channel: [val(meta), path(tsv.gz.tbi)]
    gens_cov_bed        = PREPARE_AMBER_COBALT_GENS.out.gens_cov_bed        // channel: [val(meta), path(bed.gz)]
    gens_cov_bed_tbi    = PREPARE_AMBER_COBALT_GENS.out.gens_cov_bed_tbi    // channel: [val(meta), path(bed.gz.tbi)]
    html_report = RMARKDOWNNOTEBOOK.out.report         // channel: [val(meta), path(*.html)]
}
