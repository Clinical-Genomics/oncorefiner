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

include { channelFromMetaAndPath  } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { makeMetadata            } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { ONCOREFINER             } from './workflows/oncorefiner'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references'
include { SAMTOOLS_VIEW           } from './modules/nf-core/samtools/view/main'
include { samplesheetToList       } from 'plugin/nf-schema'

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
    val_bam_normal                  // string:  [optional]  path to BAM file for the normal sample
    val_bai_normal                  // string:  [optional]  path to BAI file for the normal sample
    val_bam_tumor                   // string:  [optional]  path to BAM file for the tumor sample
    val_bai_tumor                   // string:  [optional]  path to BAI file for the tumor sample
    val_cadd_prescored_indels       // string:  [optional]  path to CADD prescored indels file
    val_cadd_resources              // string:  [optional]  path to CADD resources directory
    val_case_id                     // string:  [mandatory] case ID (used in channel metadata)
    val_genmod_score_config         // string:  [optional]  path to Genmod score config file
    val_genome                      // string:  [optional]  genome assembly (e.g. "GRCh38")
    val_genome_fasta                // string:  [optional]  path to genome fasta file
    val_genome_fai                  // string:  [optional]  path to genome fasta index file
    val_linx_breakends_tsv          // string:  [optional]  path to LINX breakends tsv file
    val_linx_fusion_tsv             // string:  [optional]  path to LINX fusion tsv file
    val_linx_sv_tsv                 // string:  [optional]  path to LINX sv tsv file
    val_multiqc_config              // string:  [optional]  path to multiqc config file
    val_multiqc_logo                // string:  [optional]  path to image file to be used as logo in multiqc report
    val_multiqc_methods_description // string:  [optional]  path to text file containing methods description to be included in multiqc report
    val_outdir                      // string:  [mandatory] path to output directory (default: ./results)
    val_sample_id_normal            // string:  [optional]  sample ID for the normal sample (used in channel metadata)
    val_sample_id_tumor             // string:  [optional]  sample ID for the tumor sample (used in channel metadata)
    val_sex                         // string:  [optional]  sex of the patient (used in channel metadata)
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

    // Initialise input channels

    def metadata_case_file = makeMetadata(val_case_id, val_case_id)

    def metadata_normal_sample_file = makeMetadata(
        val_sample_id_normal,
        val_case_id,
        val_sample_id_normal,
        "normal",
        val_sex
    )

    def metadata_tumor_sample_file = makeMetadata(
        val_sample_id_tumor,
        val_case_id,
        val_sample_id_tumor,
        "tumor",
        val_sex
    )

    ch_snv_vcf         = channelFromMetaAndPath(metadata_case_file, val_snv_vcf)
    ch_snv_vcf_tbi     = channelFromMetaAndPath(metadata_case_file, val_snv_vcf + '.tbi')
    ch_sv_vcf          = channelFromMetaAndPath(metadata_case_file, val_sv_vcf)
    ch_sv_vcf_tbi      = channelFromMetaAndPath(metadata_case_file, val_sv_vcf + '.tbi')
    ch_vep_extra_files = channel.empty()

    // Alignment files
    def ch_bam_normal     = channelFromMetaAndPath(metadata_normal_sample_file, val_bam_normal)
    def ch_bai_normal     = channelFromMetaAndPath(metadata_normal_sample_file, val_bai_normal)
    def ch_bam_bai_normal = ch_bam_normal.join(ch_bai_normal, failOnMismatch: true, failOnDuplicate: true)

    def ch_bam_tumor     = channelFromMetaAndPath(metadata_tumor_sample_file, val_bam_tumor)
    def ch_bai_tumor     = channelFromMetaAndPath(metadata_tumor_sample_file, val_bai_tumor)
    def ch_bam_bai_tumor = ch_bam_tumor.join(ch_bai_tumor, failOnMismatch: true, failOnDuplicate: true)

    // Reference files
    def ch_genome_fasta     = channelFromMetaAndPath(metadata_case_file, val_genome_fasta)
    def ch_genome_fai       = channelFromMetaAndPath(metadata_case_file, val_genome_fai)
    def ch_genome_fasta_fai = ch_genome_fasta.join(ch_genome_fai, failOnMismatch: true, failOnDuplicate: true)

    // Input for VCF annotation with LINX
    ch_linx_breakends_tsv = channelFromMetaAndPath(metadata_case_file, val_linx_breakends_tsv)
    ch_linx_fusion_tsv    = channelFromMetaAndPath(metadata_case_file, val_linx_fusion_tsv)
    ch_linx_sv_tsv        = channelFromMetaAndPath(metadata_case_file, val_linx_sv_tsv)
    ch_sv_header          = channelFromMetaAndPath(metadata_case_file, "$projectDir/assets/sv_annotation_header.txt")

    // Input for VCF annotation with LINX
    ch_linx_breakends_tsv = channel.fromPath(val_linx_breakends_tsv).map { tsv -> [[id:tsv.simpleName], tsv] }.collect()
    ch_linx_fusion_tsv    = channel.fromPath(val_linx_fusion_tsv).map { tsv -> [[id:tsv.simpleName], tsv] }.collect()
    ch_linx_sv_tsv        = channel.fromPath(val_linx_sv_tsv).map { tsv -> [[id:tsv.simpleName], tsv] }.collect()

    ch_sv_header = channel.fromPath("$projectDir/assets/sv_annotation_header.txt", checkIfExists: true).collect()

    // CADD input files
    def ch_cadd_header           = channelFromMetaAndPath(metadata_case_file, "$projectDir/assets/cadd_to_vcf_header.txt")
    def ch_cadd_resources        = channelFromMetaAndPath(metadata_case_file, val_cadd_resources)
    def ch_cadd_prescored_indels = channelFromMetaAndPath(metadata_case_file, val_cadd_prescored_indels)

    // Input for VEP
    ch_vep_extra_files = val_vep_plugin_files ? channel.fromList(samplesheetToList(val_vep_plugin_files, 'assets/vep_plugin_files_schema.json')).collect()
                                              : channel.value([])
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
    ch_sv_dbs            = val_svdb_query_dbs    ? channel.fromList(samplesheetToList(val_svdb_query_dbs, 'assets/svdb_query_vcf_schema.json'))
                                                 : channel.empty()

    // Input for genmod_score
    if (val_genmod_score_config) {
        ch_genmod_score_config = channel.fromPath(val_genmod_score_config).map { it -> [[id:it.simpleName], it] }.collect()
        val_run_genmod_score = true
    } else {
        ch_genmod_score_config = channel.empty()
        val_run_genmod_score = false
    }

    ONCOREFINER (
        ch_bam_bai_normal,
        ch_bam_bai_tumor,
        ch_cadd_header,
        ch_cadd_prescored_indels,
        ch_cadd_resources,
        ch_genmod_score_config,
        ch_genome_fasta,
        ch_genome_fai,
        ch_linx_breakends_tsv,
        ch_linx_fusion_tsv,
        ch_linx_sv_tsv,
        ch_snv_vcf,
        ch_snv_vcf_tbi,
        ch_sv_dbs,
        ch_sv_header,
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
        val_run_genmod_score,
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
    cadd_annotated_vcf        = ONCOREFINER.out.cadd_annotated_vcf        // channel: [val(meta), path(vcf)]
    cadd_annotated_tbi        = ONCOREFINER.out.cadd_annotated_tbi        // channel: [val(meta), path(tbi)]
    multiqc_report            = ONCOREFINER.out.multiqc_report            // channel: /path/to/multiqc_report.html
    snv_clinical_filtered_vcf = ONCOREFINER.out.snv_clinical_filtered_vcf // channel: [val(meta), path(vcf)]
    snv_clinical_filtered_tbi = ONCOREFINER.out.snv_clinical_filtered_tbi // channel: [val(meta), path(tbi)]
    snv_vcfanno_vcf           = ONCOREFINER.out.snv_vcfanno_vcf           // channel: [val(meta), path(vcf)]
    snv_vcfanno_tbi           = ONCOREFINER.out.snv_vcfanno_tbi           // channel: [val(meta), path(vcf.tbi)]
    snv_vep_annotated_vcf     = ONCOREFINER.out.snv_vep_annotated_vcf     // channel: [val(meta), path(vcf)]
    snv_vep_annotated_tbi     = ONCOREFINER.out.snv_vep_annotated_tbi     // channel: [val(meta), path(tbi)]
    snv_vep_report            = ONCOREFINER.out.snv_vep_report            // channel: [val(meta), val(process), val(tool), path(html)]
    snv_research_filtered_vcf = ONCOREFINER.out.snv_research_filtered_vcf // channel: [val(meta), path(vcf)]
    snv_research_filtered_tbi = ONCOREFINER.out.snv_research_filtered_tbi // channel: [val(meta), path(vcf.tbi)]
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
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    CLINICALGENOMICS_ONCOREFINER (
        params.bam_normal,
        params.bai_normal,
        params.bam_tumor,
        params.bai_tumor,
        params.cadd_prescored_indels,
        params.cadd_resources,
        params.case_id,
        params.genmod_score_config,
        params.genome,
        params.fasta,
        params.fai,
        params.linx_breakends_tsv,
        params.linx_fusion_tsv,
        params.linx_sv_tsv,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        params.sample_id_normal,
        params.sample_id_tumor,
        params.sex,
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

    //
    // WORKFLOW OUTPUTS: Group files by publish directory
    //
    ch_snv_publish = CLINICALGENOMICS_ONCOREFINER.out.cadd_annotated_vcf
        .mix(CLINICALGENOMICS_ONCOREFINER.out.cadd_annotated_tbi)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_clinical_filtered_vcf)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_clinical_filtered_tbi)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_vcfanno_vcf)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_vcfanno_tbi)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_vep_annotated_vcf)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_vep_annotated_tbi)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_vep_report)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_research_filtered_vcf)
        .mix(CLINICALGENOMICS_ONCOREFINER.out.snv_research_filtered_tbi)

    publish:
    snv = ch_snv_publish
}

output {
    snv {
        path "snv"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
