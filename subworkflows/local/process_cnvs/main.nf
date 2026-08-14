/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    ch_cnv_gene_tsv    // channel: [optional]  [val(meta), path(tsv)]
    ch_cnv_segment_tsv // channel: [optional]  [val(meta), path(tsv)]

    main:

    // Path to report template (Rmd file)
    // Todo: eventually this should be passed as an input to the workflow, but for now we can hardcode it since it's part of our assets
    def cnv_report_notebook = file("${projectDir}/assets/cnv_report.Rmd", checkIfExists: true)

    // Join the two distinct file channels together based on the meta.id
    ch_rmarkdown_in = ch_cnv_gene_tsv
        .join(ch_cnv_segment_tsv)
        .multiMap { meta, cnv_gene_tsv, cnv_segment_tsv ->
        notebook   : [meta, cnv_report_notebook]
        parameters : [cnv_gene: cnv_gene_tsv.name, cnv_segment: cnv_segment_tsv.name]
        input_files: [cnv_gene_tsv, cnv_segment_tsv]
    }


    RMARKDOWNNOTEBOOK(
    ch_rmarkdown_in.notebook,
    ch_rmarkdown_in.parameters,
    ch_rmarkdown_in.input_files,
    )

    emit:
    html_report = RMARKDOWNNOTEBOOK.out.report   // channel: [ val(meta), path(*.html) ]
}
