/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { ENSEMBLVEP_VEP as VEP_ANNOTATE_SNV } from '../modules/nf-core/ensemblvep/vep/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_postprocessing_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POSTPROCESSING {

    take:
        ch_samplesheet // channel: samplesheet read in from --input

    main:

        // Initialize input channels
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        ch_vep_snv = params.snv_vcf ?
            Channel
                .fromPath(params.snv_vcf, checkIfExists: true)
                .map { vcf ->
                    def meta = [id: vcf.simpleName]
                    def custom_extra_files = params.custom_extra_files ? file(params.custom_extra_files) : []
                    tuple(meta, vcf, custom_extra_files)
                } :
            Channel.empty()

        // Reference files
        ch_genome_fasta              = Channel.fromPath(params.fasta).map { it -> [[id:it.simpleName], it] }.collect()
        ch_vep_cache = params.vep_cache ? Channel.fromPath(params.vep_cache, checkIfExists: true) : Channel.empty()


        //
        // Read and store paths in the vep_plugin_files file
        //
        ch_vep_extra_files_unsplit  = params.vep_plugin_files ? Channel.fromPath(params.vep_plugin_files).collect() : Channel.value([])
        ch_vep_extra_files = Channel.empty()
        if (params.vep_plugin_files) {
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE SNVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
        // Process SNV VCF files
        if (params.snv_vcf) {
            // Always annotate with VEP
            VEP_ANNOTATE_SNV (
                ch_vep_snv,
                params.genome,
                "homo_sapiens",
                params.vep_cache_version,
                ch_vep_cache,
                ch_genome_fasta,
                ch_vep_extra_files
            )
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'nf_core_'  +  'postprocessing_software_'  + 'mqc_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }

        //
        // MODULE: MultiQC
        //
        ch_multiqc_config        = Channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )


        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
