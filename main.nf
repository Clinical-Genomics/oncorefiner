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
    val_genome                  // string:  [optional]  genome assembly (e.g. "GRCh38")
    val_genome_fasta            // string:  [optional]  path to genome fasta file
    val_vep_cache               // string:  [optional]  path to vep cache tar gzip file
    val_vep_plugin_files        // string:  [optional]  path to file containing paths to vep plugin files (one per line)
    val_vcfanno_resources       // string:  [optional]  path to file containing paths to vcfanno resources (one per line)
    val_vcfanno_lua             // string:  [optional]  path to vcfanno lua file
    val_vcfanno_toml            // string:  [optional]  path to vcfanno toml file
    val_vcfanno_extra           // string:  [optional]  path to file containing paths to extra files for vcfanno (one per line)
    val_svdb_query_dbs          // string:  [optional]  path to file containing paths to SVDB query databases and additional information (one per line)


    main:

    // Initialize input channels for oncorefiner
    ch_snv_vcf              = channel.fromPath(val_snv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_snv_vcf_tbi          = channel.fromPath(val_snv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf               = channel.fromPath(val_sv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf_tbi           = channel.fromPath(val_sv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_vep_extra_files      = channel.empty()
    ch_svdb_dbs             = channel.empty()

    // Reference files
    ch_genome_fasta         = channel.fromPath(val_genome_fasta).map { it -> [[id:it.simpleName], it] }.collect()

    // File channels for PREPARE_REFERENCES
    ch_vep_cache_unprocessed     = val_vep_cache           ? channel.fromPath(val_vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                           : channel.value([[],[]])

    // VEP: Parse paths in the file 'vep_plugin_files' to create ch_vep_extra_files
    ch_vep_extra_files_unsplit  = val_vep_plugin_files ? channel.fromPath(val_vep_plugin_files).collect() : channel.value([])
    if (val_vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
            .map { row ->
                def f = file(row.vep_files[0])
                if(f.isFile() || f.isDirectory()){
                    return [f]
                } else {
                    error("\nVep database file ${f} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    // Vcfanno: Initialize channels for vcfanno resources, lua and toml files, and extra files
    ch_vcfanno_resources = val_vcfanno_resources ? channel.fromPath(val_vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                 : channel.value([])
    ch_vcfanno_lua       = val_vcfanno_lua       ? channel.fromPath(val_vcfanno_lua).collect()
                                                 : channel.value([])
    ch_vcfanno_toml      = val_vcfanno_toml      ? channel.fromPath(val_vcfanno_toml).collect()
                                                 : channel.value([])
    ch_vcfanno_extra     = val_vcfanno_extra     ? channel.fromPath(val_vcfanno_extra).collect()
                                                 : []

    // SVDB: Initialize channel for SVDB query csv file
    ch_sv_dbs            = val_svdb_query_dbs    ? channel.fromPath(val_svdb_query_dbs)
                                                 : channel.empty()

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
        PREPARE_REFERENCES.out.vep_resources,
        ch_vep_extra_files,
        ch_vcfanno_resources,
        ch_vcfanno_lua,
        ch_vcfanno_toml,
        ch_vcfanno_extra,
        ch_sv_dbs,
        val_genome
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
        params.vep_cache,
        params.vep_plugin_files,
        params.vcfanno_resources,
        params.vcfanno_lua,
        params.vcfanno_toml,
        params.vcfanno_extra,
        params.svdb_query_dbs,
        params.genome
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
