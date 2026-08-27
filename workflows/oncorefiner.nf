/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PROCESS_SNVS           } from '../subworkflows/local/process_snvs/main.nf'
include { PROCESS_SVS            } from '../subworkflows/local/process_svs/main.nf'
include { PROCESS_CNVS           } from '../subworkflows/local/process_cnvs/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOREFINER {

    take:
    ch_amber_baf_tsv_gz             // channel: [optional]  [val(meta), path(tsv)]
    ch_bam_bai_normal               // channel: [optional]  [val(meta), path(bam), path(bai)]
    ch_bam_bai_tumor                // channel: [mandatory] [val(meta), path(bam), path(bai)]
    ch_cadd_header                  // channel: [mandatory] [path(txt)]
    ch_cadd_prescored_indels        // channel: [optional]  [val(meta), path(dir)]
    ch_cadd_resources               // channel: [optional]  [val(meta), path(dir)]
    ch_cnv_gene_tsv                 // channel: [optional]  [val(meta), path(tsv)]
    ch_cnv_segment_tsv              // channel: [optional]  [val(meta), path(tsv)]
    ch_cobalt_ratio_pcf_normal      // channel: [optional]  [val(meta), path(pcf)]
    ch_cobalt_ratio_pcf_tumor       // channel: [optional]  [val(meta), path(pcf)]
    ch_genmod_score_config_snv      // channel: [optional]  [val(meta), path(ini)]
    ch_genmod_score_config_sv       // channel: [optional]  [val(meta), path(ini)]
    ch_genome_fasta                 // channel: [optional]  [val(meta), path(fasta)]
    ch_genome_fai                   // channel: [optional]  [val(meta), path(fai)]
    ch_linx_breakends_tsv           // channel: [optional]  [val(meta), path(tsv)]
    ch_linx_fusion_tsv              // channel: [optional]  [val(meta), path(tsv)]
    ch_linx_sv_tsv                  // channel: [optional]  [val(meta), path(tsv)]
    ch_snv_vcf                      // channel: [optional]  [val(meta), path(vcf)]
    ch_snv_vcf_tbi                  // channel: [optional]  [val(meta), path(vcf.tbi)]
    ch_sv_dbs                       // channel: [optional]  [path(csv)]
    ch_sv_header                    // channel: [optional]  [path(txt)]
    ch_sv_vcf                       // channel: [optional]  [val(meta), path(vcf)]
    ch_sv_vcf_tbi                   // channel: [optional]  [val(meta), path(vcf.tbi)]
    ch_vcfanno_extra                // channel: [optional]  [path(extra_file1), path(extra_file2), ...]
    ch_vcfanno_lua                  // channel: [optional]  [path(lua_file)]
    ch_vcfanno_resources            // channel: [optional]  [path(resource_file1), path(resource_file2), ...]
    ch_vcfanno_toml                 // channel: [optional]  [path(toml_file)]
    ch_vep_cache                    // channel: [optional]  [vep_cache_files]
    ch_vep_extra_files              // channel: [optional]  [path(plugin_file1), path(plugin_file2), ...]
    val_analysis_type               // string:  [optional]  analysis type, e.g. "tumor_only" or "tumor_normal"
    val_cadd_resources              // string:  [optional]  path to CADD resources directory
    val_genome                      // string:  [optional]  genome assembly (e.g. "GRCh38")
    val_multiqc_config              // string:  [optional]  path to multiqc config file
    val_multiqc_logo                // string:  [optional]  path to image file to be used as logo in multiqc report
    val_multiqc_methods_description // string:  [optional]  path to text file containing methods description to be included in multiqc report
    val_outdir                      // string:  [mandatory] path to output directory (default: ./results)
    val_run_genmod_score_snv        // boolean: [mandatory] whether to skip PROCESS_SNVS:VCF_ANNOTATE_SCORE_GENMOD process for SNVs
    val_run_genmod_score_sv         // boolean: [mandatory] whether to skip PROCESS_SVS:VCF_ANNOTATE_SCORE_GENMOD process for SVs
    val_species                     // string:  [optional]  species (e.g. "homo_sapiens")
    val_vep_cache_version           // string:  [optional]  version of vep cache to use (e.g. "107")

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // Process CNV input files
    PROCESS_CNVS (
        ch_amber_baf_tsv_gz,
        ch_cnv_gene_tsv,
        ch_cnv_segment_tsv,
        ch_cobalt_ratio_pcf_normal,
        ch_cobalt_ratio_pcf_tumor,
        val_analysis_type
    )

    // Process SNV VCF files
    PROCESS_SNVS (
        ch_genome_fasta,
        ch_genome_fai,
        ch_cadd_header,
        ch_cadd_prescored_indels,
        ch_cadd_resources,
        ch_genmod_score_config_snv,
        ch_snv_vcf,
        ch_snv_vcf_tbi,
        ch_vcfanno_extra,
        ch_vcfanno_lua,
        ch_vcfanno_resources,
        ch_vcfanno_toml,
        ch_vep_cache,
        ch_vep_extra_files,
        val_cadd_resources,
        val_genome,
        val_run_genmod_score_snv,
        val_species,
        val_vep_cache_version
    )

    // Process SV VCF files
    PROCESS_SVS(
        ch_bam_bai_normal,
        ch_bam_bai_tumor,
        ch_genmod_score_config_sv,
        ch_linx_breakends_tsv,
        ch_linx_fusion_tsv,
        ch_linx_sv_tsv,
        ch_sv_header,
        ch_sv_vcf,
        ch_sv_vcf_tbi,
        ch_sv_dbs,
        val_genome,
        val_run_genmod_score_sv,
        val_species,
        val_vep_cache_version,
        ch_vep_cache,
        ch_genome_fasta,
        ch_vep_extra_files
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name:  'oncorefiner_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description =     val_multiqc_methods_description
        ? file(val_multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'oncorefiner'],
                files,
                val_multiqc_config
                    ? file(val_multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                val_multiqc_logo ? file(val_multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    cadd_annotated_vcf        = PROCESS_SNVS.out.cadd_annotated_vcf                           // channel: [val(meta), path(vcf)]
    cadd_annotated_tbi        = PROCESS_SNVS.out.cadd_annotated_tbi                           // channel: [val(meta), path(tbi)]
    cnv_baf_normal_bed        = PROCESS_CNVS.out.gens_baf_normal_bed                          // channel: [val(meta), path(bed)]
    cnv_baf_normal_tbi        = PROCESS_CNVS.out.gens_baf_normal_tbi                          // channel: [val(meta), path(tbi)]
    cnv_baf_tumor_bed         = PROCESS_CNVS.out.gens_baf_tumor_bed                           // channel: [val(meta), path(bed)]
    cnv_baf_tumor_tbi         = PROCESS_CNVS.out.gens_baf_tumor_tbi                           // channel: [val(meta), path(tbi)]
    cnv_cov_bed               = PROCESS_CNVS.out.gens_cov_bed                                 // channel: [val(meta), path(bed)]
    cnv_cov_bed_tbi           = PROCESS_CNVS.out.gens_cov_bed_tbi                             // channel: [val(meta), path(tbi)]
    cnv_report_html           = PROCESS_CNVS.out.html_report                                  // channel: [val(meta), path(html)]
    multiqc_data              = MULTIQC.out.data                                              // channel: /path/to/multiqc_data/
    multiqc_plots             = MULTIQC.out.plots                                             // channel: /path/to/multiqc_plots/
    multiqc_report            = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    snv_clinical_filtered_vcf = PROCESS_SNVS.out.clinical_filtered_vcf                        // channel: [val(meta), path(vcf)]
    snv_clinical_filtered_tbi = PROCESS_SNVS.out.clinical_filtered_tbi                        // channel: [val(meta), path(tbi)]
    snv_vcfanno_vcf           = PROCESS_SNVS.out.vcfanno_vcf                                  // channel: [val(meta), path(vcf)]
    snv_vcfanno_tbi           = PROCESS_SNVS.out.vcfanno_tbi                                  // channel: [val(meta), path(vcf.tbi)]
    snv_vep_annotated_vcf     = PROCESS_SNVS.out.vep_annotated_vcf                            // channel: [val(meta), path(vcf)]
    snv_vep_annotated_tbi     = PROCESS_SNVS.out.vep_annotated_tbi                            // channel: [val(meta), path(tbi)]
    snv_vep_report            = PROCESS_SNVS.out.vep_report                                   // channel: [val(meta), val(process), val(tool), path(html)]
    snv_research_filtered_vcf = PROCESS_SNVS.out.research_filtered_vcf                        // channel: [val(meta), path(vcf)]
    snv_research_filtered_tbi = PROCESS_SNVS.out.research_filtered_tbi                        // channel: [val(meta), path(vcf.tbi)]
    sv_clinical_filtered_vcf  = PROCESS_SVS.out.clinical_filtered_vcf                         // channel: [val(meta), path(vcf)]
    sv_clinical_filtered_tbi  = PROCESS_SVS.out.clinical_filtered_tbi                         // channel: [val(meta), path(tbi)]
    sv_research_filtered_vcf  = PROCESS_SVS.out.research_filtered_vcf                         // channel: [val(meta), path(vcf)]
    sv_research_filtered_tbi  = PROCESS_SVS.out.research_filtered_tbi                         // channel: [val(meta), path(vcf.tbi)]
    sv_vcf2cytosure_cgh       = PROCESS_SVS.out.vcf2cytosure_cgh                              // channel: [val(meta), path(cgh)]
    sv_vep_annotated_vcf      = PROCESS_SVS.out.vep_annotated_vcf                             // channel: [val(meta), path(vcf)]
    sv_vep_annotated_tbi      = PROCESS_SVS.out.vep_annotated_tbi                             // channel: [val(meta), path(tbi)]
    sv_vep_report             = PROCESS_SVS.out.vep_report                                    // channel: [val(meta), val(process), val(tool), path(html)]
    versions                  = ch_versions                                                   // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
