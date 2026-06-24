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
    RUN PROCESS_SNVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_CNVS {

    take:
    ch_cnv_gene_tsv    // channel: [optional]  [val(meta), path(tsv)]
    ch_cnv_segment_tsv // channel: [optional]  [val(meta), path(tsv)]

    main:

    // Path to report template (Rmd file)
    // Todo: eventually this should be passed as an input to the workflow, but for now we can hardcode it since it's part of our assets
    def cnv_report_template = file("${projectDir}/assets/cnv_report.Rmd", checkIfExists: true)

    // Join the two distinct file channels together based on the meta.id
    ch_report_inputs = ch_cnv_gene_tsv
            .join(ch_cnv_segment_tsv)
            .map { meta, cnv_gene_tsv, cnv_segment_tsv ->
                // Define any internal variables to pass to the Rmd document params
                def r_params = [
                    cnv_gene   : cnv_gene_tsv.name,
                    cnv_segment: cnv_segment_tsv.name,
                ]

                // Group the 2 target data files into a single list element
                def data_files = [cnv_gene_tsv, cnv_segment_tsv]

                return [ [meta, cnv_report_template], r_params, data_files ]
            }

    RMARKDOWNNOTEBOOK (
        ch_report_inputs.map { it[0] }, // tuple val(meta), path(notebook)
        ch_report_inputs.map { it[1] }, // val parameters
        ch_report_inputs.map { it[2] }  // path input_files
    )

    emit:
    html_report = RMARKDOWNNOTEBOOK.out.report   // channel: [ val(meta), path(*.html) ]
}
