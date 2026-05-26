/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { RMARKDOWNNOTEBOOK                       } from '../../../modules/nf-core/rmarkdownnotebook/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_SNVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_CNVS {

    take:
    ch_cnv_gene_tsv                 // channel: [optional]  [val(meta), path(tsv)]
    ch_cnv_segment_tsv              // channel: [optional]  [val(meta), path(tsv)]  
    
    main:
    // Path to report template (Rmd file)
    def cnv_report_template = file("${projectDir}/assets/cnv_report.Rmd", checkIfExists: true)

    // Join the two distinct file channels together based on the meta.id
    ch_combined_inputs = ch_cnv_gene_tsv.join(ch_cnv_segment_tsv)

    // Prepare the exact tuple format expected by the module's input:
    // Input 1: tuple val(meta), path(notebook)
    // Input 2: val(parameters)
    // Input 3: path(input_files)
    ch_report_inputs = ch_combined_inputs.map { meta, cnv_gene_tsv, cnv_segment_tsv ->
        
        // Define any internal variables you want to pass to the Rmd document params
        def r_params = [
            cnv_gene   : cnv_gene_tsv.name,
            cnv_segment: cnv_segment_tsv.name,
            output_file: "${meta}_cnv_report"
        ]

        // Group the 2 target data files into a single list element
        def data_files = [cnv_gene_tsv, cnv_segment_tsv]

        return [ [meta, cnv_report_template], r_params, data_files ]
    }

    // -------------------------------------------------------------
    // 2. Invoke the Module
    // -------------------------------------------------------------
    RMARKDOWNNOTEBOOK (
        ch_report_inputs.map { it[0] }, // tuple val(meta), path(notebook)
        ch_report_inputs.map { it[1] }, // val parameters
        ch_report_inputs.map { it[2] }  // path input_files
    )

    emit:
    html_report = RMARKDOWNNOTEBOOK.out.report       // channel: [ val(meta), path(*.html) ]
    versions    = RMARKDOWNNOTEBOOK.out.versions     // channel: [ path(versions.yml) ]
}