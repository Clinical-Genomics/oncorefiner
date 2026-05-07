/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                                  } from '../modules/nf-core/multiqc/main'


//
// SUBWORKFLOWS
//

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oncorefiner_pipeline'

//
// LOCAL SUBWORKFLOWS
//

include { PROCESS_SNVS } from '../subworkflows/local/process_snvs/main.nf'
include { PROCESS_SVS } from '../subworkflows/local/process_svs/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOREFINER {

    take:
        ch_samplesheet           // channel: [mandatory] samplesheet read in from --input
        ch_bam_bai_normal        // channel: [optional]  [val(meta), path(bam), path(bai)]
        ch_bam_bai_tumor         // channel: [mandatory] [val(meta), path(bam), path(bai)]
        ch_cadd_header           // channel: [mandatory] [path(txt)]
        ch_cadd_prescored_indels // channel: [optional]  [val(meta), path(dir)]
        ch_cadd_resources        // channel: [optional]  [val(meta), path(dir)]
        ch_genome_fasta          // channel: [optional]  [val(meta), path(fasta)]
        ch_genome_fai            // channel: [optional]  [val(meta), path(fai)]
        ch_snv_vcf               // channel: [optional]  [val(meta), path(vcf)]
        ch_snv_vcf_tbi           // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_sv_dbs                // channel: [optional]  [path(csv)]
        ch_sv_vcf                // channel: [optional]  [val(meta), path(vcf)]
        ch_sv_vcf_tbi            // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_vcfanno_extra         // channel: [optional]  [path(extra_file1), path(extra_file2), ...]
        ch_vcfanno_lua           // channel: [optional]  [path(lua_file)]
        ch_vcfanno_resources     // channel: [optional]  [path(resource_file1), path(resource_file2), ...]
        ch_vcfanno_toml          // channel: [optional]  [path(toml_file)]
        ch_vep_cache             // channel: [optional]  [vep_cache_files]
        ch_vep_extra_files       // channel: [optional]  [path(plugin_file1), path(plugin_file2), ...]
        val_cadd_resources       // string:  [optional]  path to CADD resources directory
        val_genome               // string:  [optional]  genome assembly (e.g. "GRCh38")
        val_species              // string:  [optional]  species (e.g. "homo_sapiens")
        val_vep_cache_version    // string:  [optional]  version of vep cache to use (e.g. "107")

    main:

        // Initialize input channels
        ch_versions             = channel.empty()
        ch_multiqc_files        = channel.empty()

        // Process SNV VCF files
        PROCESS_SNVS (
            ch_genome_fasta,
            ch_genome_fai,
            ch_cadd_header,
            ch_cadd_prescored_indels,
            ch_cadd_resources,
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
            val_species,
            val_vep_cache_version
        )

        // Process SV VCF files
        PROCESS_SVS(
            ch_bam_bai_normal,
            ch_bam_bai_tumor,
            ch_sv_vcf,
            ch_sv_vcf_tbi,
            ch_sv_dbs,
            val_genome,
            val_species,
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

        softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
            .mix(topic_versions_string)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'oncorefiner_software_'  + 'mqc_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }

        //
        // MODULE: MultiQC
        //
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: ''],
                    files,
                    params.multiqc_config
                    ? file(params.multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )

    emit:
        cadd_annotated_vcf        = PROCESS_SNVS.out.cadd_annotated_vcf    // channel: [val(meta), path(vcf)]
        cadd_annotated_tbi        = PROCESS_SNVS.out.cadd_annotated_tbi    // channel: [val(meta), path(tbi)]
        multiqc_report            = MULTIQC.out.report.toList()            // channel: /path/to/multiqc_report.html
        snv_clinical_filtered_vcf = PROCESS_SNVS.out.clinical_filtered_vcf // channel: [val(meta), path(vcf)]
        snv_clinical_filtered_tbi = PROCESS_SNVS.out.clinical_filtered_tbi // channel: [val(meta), path(tbi)]
        snv_vcfanno_vcf           = PROCESS_SNVS.out.vcfanno_vcf           // channel: [val(meta), path(vcf)]
        snv_vcfanno_tbi           = PROCESS_SNVS.out.vcfanno_tbi           // channel: [val(meta), path(vcf.tbi)]
        snv_vep_annotated_vcf     = PROCESS_SNVS.out.vep_annotated_vcf     // channel: [val(meta), path(vcf)]
        snv_vep_annotated_tbi     = PROCESS_SNVS.out.vep_annotated_tbi     // channel: [val(meta), path(tbi)]
        snv_vep_report            = PROCESS_SNVS.out.vep_report            // channel: [val(meta), val(process), val(tool), path(html)]
        snv_research_filtered_vcf = PROCESS_SNVS.out.research_filtered_vcf // channel: [val(meta), path(vcf)]
        snv_research_filtered_tbi = PROCESS_SNVS.out.research_filtered_tbi // channel: [val(meta), path(vcf.tbi)]
        versions                  = ch_versions                            // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
