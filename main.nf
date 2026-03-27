#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Clinical-Genomics/oncorefiner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Clinical-Genomics/oncorefiner
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ONCOREFINER             } from './workflows/oncorefiner'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CLINICALGENOMICS_ONCOREFINER {

    take:
    samplesheet                 // channel: [mandatory] samplesheet read in from --input
    val_snv_vcf                 // string:  [optional]  path to input SNV vcf file
    val_sv_vcf                  // string:  [optional]  path to input SV vcf file
    val_genome_fasta            // string:  [optional]  path to genome fasta file
    val_vep_cache               // string:  [optional] path to vep cache tar gzip file




    main:

    // Initialize input channels
    ch_snv_vcf              = channel.fromPath(val_snv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_snv_vcf_tbi          = channel.fromPath(val_snv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf               = channel.fromPath(val_sv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf_tbi           = channel.fromPath(val_sv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()

    // Reference files
    ch_genome_fasta         = channel.fromPath(val_genome_fasta).map { it -> [[id:it.simpleName], it] }.collect()

    // File channels for PREPARE_REFERENCES
    ch_vep_cache_unprocessed     = val_vep_cache           ? channel.fromPath(val_vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                           : channel.value([[],[]])



    //
    // Subworkflow: Prepare reference files
    //

    PREPARE_REFERENCES (
        params.vep_cache
        )

    //
    // WORKFLOW: Run pipeline
    //
    ONCOREFINER (
        samplesheet,
        ch_snv_vcf,
        ch_snv_vcf_tbi,
        ch_sv_vcf,
        ch_sv_vcf_tbi,
        ch_genome_fasta,
        PREPARE_REFERENCES.out.vep_resources
    )
    emit:
    multiqc_report = ONCOREFINER.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    CLINICALGENOMICS_ONCOREFINER (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.snv_vcf,
        params.sv_vcf,
        params.fasta,
        params.vep_cache
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        CLINICALGENOMICS_ONCOREFINER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
