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
include { SVDB_QUERY as SVDB_QUERY_DB              } from '../modules/nf-core/svdb/query/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SV          } from '../modules/nf-core/ensemblvep/vep/main'
include { BCFTOOLS_VIEW as RESEARCH_FILTERING_SV   } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as CLINICAL_FILTERING_SV   } from '../modules/nf-core/bcftools/view/main'

//
// MODULE: Local modules
//

include { GENERATE_CYTOSURE_FILES } from '../subworkflows/local/generate_cytosure_files/main'

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOREFINER {

    take:
        ch_samplesheet        // channel: [mandatory] samplesheet read in from --input
        ch_bam_bai_normal     // channel: [optional]  [val(meta), path(bam), path(bai)]
        ch_bam_bai_tumor      // channel: [mandatory]  [val(meta), path(bam), path(bai)]
        ch_genome_fasta       // channel: [optional]  [val(meta), path(fasta)]
        ch_snv_vcf            // channel: [optional]  [val(meta), path(vcf)]
        ch_snv_vcf_tbi        // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_sv_dbs             // channel: [optional]  [path(csv)]
        ch_sv_vcf             // channel: [optional]  [val(meta), path(vcf)]
        ch_sv_vcf_tbi         // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_vcfanno_extra      // channel: [optional]  [path(extra_file1), path(extra_file2), ...]
        ch_vcfanno_lua        // channel: [optional]  [path(lua_file)]
        ch_vcfanno_resources  // channel: [optional]  [path(resource_file1), path(resource_file2), ...]
        ch_vcfanno_toml       // channel: [optional]  [path(toml_file)]
        ch_vep_cache          // channel: [optional]  [vep_cache_files]
        ch_vep_extra_files    // channel: [optional]  [path(plugin_file1), path(plugin_file2), ...]
        val_genome            // string:  [optional]  genome assembly (e.g. "GRCh38")
        val_species           // string:  [optional]  species (e.g. "homo_sapiens")
        val_vep_cache_version // string:  [optional]  version of vep cache to use (e.g. "107")

    main:

        // Initialize input channels
        ch_versions             = channel.empty()
        ch_multiqc_files        = channel.empty()

        // Process SNV VCF files
        if (ch_snv_vcf) {
            PROCESS_SNVS (
                ch_genome_fasta,
                ch_snv_vcf,
                ch_snv_vcf_tbi,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml,
                ch_vep_cache,
                ch_vep_extra_files,
                val_genome,
                val_species,
                val_vep_cache_version
            )
        }

        // Process SV VCF files
        if (ch_sv_vcf) {

            // VCF2CYTOSURE
            ch_bam_bai = channel.empty().mix(ch_bam_bai_tumor, ch_bam_bai_normal)
            ch_vcf2cytosure_in = ch_bam_bai.combine(
                ch_sv_vcf.join(ch_sv_vcf_tbi, failOnMismatch: true),
                )
                .multiMap { meta_bam_bai, bam, bai, meta_vcf, vcf, tbi ->
                    bam_bai: tuple(meta_bam_bai, bam, bai)
                    vcf: tuple(meta_vcf, vcf)
                    tbi: tuple(meta_vcf, tbi)
                }

            GENERATE_CYTOSURE_FILES (
                ch_vcf2cytosure_in.bam_bai,
                ch_vcf2cytosure_in.tbi,
                ch_vcf2cytosure_in.vcf
            )

            // SVDB QUERY
            ch_sv_dbs
                .splitCsv ( header:true )
                .multiMap { row ->
                    vcf_dbs:  row.filename
                    in_frqs:  row.in_freq_info_key
                    in_occs:  row.in_allele_count_info_key
                    out_frqs: row.out_freq_info_key
                    out_occs: row.out_allele_count_info_key
                }
                .set { ch_svdb_dbs }

            SVDB_QUERY_DB (
                ch_sv_vcf,
                ch_svdb_dbs.in_occs.toList(),
                ch_svdb_dbs.in_frqs.toList(),
                ch_svdb_dbs.out_occs.toList(),
                ch_svdb_dbs.out_frqs.toList(),
                ch_svdb_dbs.vcf_dbs.toList(),
                []
            )


            // Quality Filtering
            SVDB_QUERY_DB.out.vcf
                .map { meta, vcf ->
                    tuple(meta, vcf, []) }
                .set { ch_research_filtering_sv_in }

            RESEARCH_FILTERING_SV(ch_research_filtering_sv_in, [], [], [])

            // VEP
            RESEARCH_FILTERING_SV.out.vcf
                .map { meta, vcf ->
                            tuple(meta, vcf, []) }
                .set { ch_vep_sv }


            ENSEMBLVEP_SV(
                ch_vep_sv,
                val_genome,
                val_species,
                val_vep_cache_version,
                ch_vep_cache,
                ch_genome_fasta,
                ch_vep_extra_files
            )

            // Clinical Filtering
            ENSEMBLVEP_SV.out.vcf
                .join(ENSEMBLVEP_SV.out.tbi)
                .map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi)
                    }
                .set { ch_clinical_filtering_sv_in }
            CLINICAL_FILTERING_SV(ch_clinical_filtering_sv_in, [], [], [])


        }


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
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
