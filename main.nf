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
include { SAMTOOLS_VIEW           } from './modules/nf-core/samtools/view/main'
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
    samplesheet                     // channel: [mandatory] samplesheet read in from --input
    val_bam_normal                  // string:  [optional]  path to BAM file for the normal sample
    val_bai_normal                  // string:  [optional]  path to BAI file for the normal sample
    val_bam_tumor                   // string:  [optional]  path to BAM file for the tumor sample
    val_bai_tumor                   // string:  [optional]  path to BAI file for the tumor sample
    val_cadd_prescored_indels       // string:  [optional]  path to CADD prescored indels file
    val_cadd_resources              // string:  [optional]  path to CADD resources directory
    val_genmod_score_config         // string:  [optional]  path to Genmod score config file
    val_genome                      // string:  [optional]  genome assembly (e.g. "GRCh38")
    val_genome_fasta                // string:  [optional]  path to genome fasta file
    val_genome_fai                  // string:  [optional]  path to genome fasta index file
    val_multiqc_config              // string:  [optional]  path to multiqc config file
    val_multiqc_logo                // string:  [optional]  path to image file to be used as logo in multiqc report
    val_multiqc_methods_description // string:  [optional]  path to text file containing methods description to be included in multiqc report
    val_outdir                      // string:  [mandatory] path to output directory (default: ./results)
    val_snv_vcf                     // string:  [optional]  path to input SNV vcf file
    val_species                     // string:  [optional]  species (e.g. "homo_sapiens")
    val_sv_vcf                      // string:  [optional]  path to input SV vcf file
    val_svdb_query_dbs              // string:  [optional]  path to file containing paths to SVDB query databases and additional information (one per line)
    val_vcfanno_extra               // string:  [optional]  path to file containing paths to extra files for vcfanno (one per line)
    val_vcfanno_lua                 // string:  [optional]  path to vcfanno lua file
    val_vcfanno_resources           // string:  [optional]  path to file containing paths to vcfanno resources (one per line)
    val_vcfanno_toml                // string:  [optional]  path to vcfanno toml file
    val_vep_cache                   // string:  [optional]  path to vep cache tar gzip file
    val_vep_cache_version           // string:  [optional]  version of vep cache to use (e.g. "107")
    val_vep_plugin_files            // string:  [optional]  path to file containing paths to vep plugin files (one per line)

    main:

    //
    // Subworkflow: Prepare reference files
    //
    PREPARE_REFERENCES (
        val_vep_cache
        )

    //
    // WORKFLOW: Run pipeline
    //

    // Input channels
    ch_snv_vcf         = channel.fromPath(val_snv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_snv_vcf_tbi     = channel.fromPath(val_snv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf          = channel.fromPath(val_sv_vcf).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_sv_vcf_tbi      = channel.fromPath(val_sv_vcf + '.tbi', checkIfExists: true).map { vcf -> [[id:vcf.simpleName], vcf] }.collect()
    ch_vep_extra_files = channel.empty()
    ch_svdb_dbs        = channel.empty()

    // Alignment files
    ch_bam_bai_normal  = channel.empty()

    if (val_bam_normal && val_bai_normal) {
        ch_bam_bai_normal = channel.fromPath(val_bam_normal)
                            .combine(channel.fromPath(val_bai_normal))
                            .map { bam, bai -> [[type:'normal'], bam, bai] }
    }

    ch_bam_bai_tumor = channel.empty()

    if (val_bam_tumor && val_bai_tumor) {
        ch_bam_bai_tumor = channel.fromPath(val_bam_tumor)
                            .combine(channel.fromPath(val_bai_tumor))
                            .map { bam, bai -> [[type:'tumor'], bam, bai] }
    }

    // Reference files
    ch_genome_fasta          = channel.fromPath(val_genome_fasta).map { it -> [[id:it.simpleName], it] }.collect()
    ch_genome_fai            = channel.fromPath(val_genome_fai).map { it -> [[id:it.simpleName], it] }.collect()
    ch_genome_fasta_fai      = ch_genome_fasta.join(ch_genome_fai, failOnMismatch: true, failOnDuplicate: true)

    // CADD input files
    ch_cadd_header           = channel.fromPath("$projectDir/assets/cadd_to_vcf_header.txt", checkIfExists: true).collect()
    ch_cadd_resources        = val_cadd_resources        ? channel.fromPath(val_cadd_resources).map { it -> [[id:'cadd_resources'], it] }.collect()
                                                         : channel.value([])

    ch_cadd_prescored_indels = val_cadd_prescored_indels ? channel.fromPath(val_cadd_prescored_indels).map { it -> [[id:'cadd_prescored_indels'], it] }.collect()
                                                         : channel.value([])

    // Input for VEP
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

    // Input for Vcfanno
    ch_vcfanno_extra     = val_vcfanno_extra     ? channel.fromPath(val_vcfanno_extra).collect()
                                                 : []
    ch_vcfanno_lua       = val_vcfanno_lua       ? channel.fromPath(val_vcfanno_lua).collect()
                                                 : channel.value([])
    ch_vcfanno_resources = val_vcfanno_resources ? channel.fromPath(val_vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                 : channel.value([])
    ch_vcfanno_toml      = val_vcfanno_toml      ? channel.fromPath(val_vcfanno_toml).collect()
                                                 : channel.value([])


    // Input for SVDB
    ch_sv_dbs            = val_svdb_query_dbs    ? channel.fromPath(val_svdb_query_dbs)
                                                 : channel.empty()

    // Input for genmod_score
    ch_genmod_score_config = val_genmod_score_config && val_snv_vcf ?
                                                     channel.fromPath(val_genmod_score_config)
                                                     .combine(ch_snv_vcf)
                                                     .map { genmod_score_config, meta, _snv_vcf -> [meta, genmod_score_config] }
                                                     : channel.empty()

    ONCOREFINER (
        samplesheet,
        ch_bam_bai_normal,
        ch_bam_bai_tumor,
        ch_cadd_header,
        ch_cadd_prescored_indels,
        ch_cadd_resources,
        ch_genmod_score_config,
        ch_genome_fasta,
        ch_genome_fai,
        ch_snv_vcf,
        ch_snv_vcf_tbi,
        ch_sv_dbs,
        ch_sv_vcf,
        ch_sv_vcf_tbi,
        ch_vcfanno_extra,
        ch_vcfanno_lua,
        ch_vcfanno_resources,
        ch_vcfanno_toml,
        PREPARE_REFERENCES.out.vep_resources,
        ch_vep_extra_files,
        val_cadd_resources,
        val_genome,
        val_multiqc_config,
        val_multiqc_logo,
        val_multiqc_methods_description,
        val_outdir,
        val_species,
        val_vep_cache_version,
    )

    //
    // Convert BAM to CRAM
    //
    ch_bam_bai = ch_bam_bai_tumor.mix(ch_bam_bai_normal)

    ch_samtools_in = ch_bam_bai.combine(ch_genome_fasta_fai)
            .multiMap { meta_bam_bai, bam, bai, meta_fasta_fai, fasta, fai ->
                bam_bai: tuple(meta_bam_bai, bam, bai)
                fasta_fai: tuple(meta_fasta_fai, fasta, fai)}

    SAMTOOLS_VIEW ( ch_samtools_in.bam_bai, ch_samtools_in.fasta_fai, [[], []], [[],[]], 'crai' )

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
        params.bam_normal,
        params.bai_normal,
        params.bam_tumor,
        params.bai_tumor,
        params.cadd_prescored_indels,
        params.cadd_resources,
        params.genmod_score_config,
        params.genome,
        params.fasta,
        params.fai,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        params.snv_vcf,
        params.species,
        params.sv_vcf,
        params.svdb_query_dbs,
        params.vcfanno_extra,
        params.vcfanno_lua,
        params.vcfanno_resources,
        params.vcfanno_toml,
        params.vep_cache,
        params.vep_cache_version,
        params.vep_plugin_files
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
        CLINICALGENOMICS_ONCOREFINER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
