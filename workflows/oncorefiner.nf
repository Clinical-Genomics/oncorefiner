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
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SNV         } from '../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                                  } from '../modules/nf-core/vcfanno/main'
include { BCFTOOLS_VIEW as RESEARCH_FILTERING      } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as CLINICAL_FILTERING      } from '../modules/nf-core/bcftools/view/main'
include { SVDB_QUERY as SVDB_QUERY_DB              } from '../modules/nf-core/svdb/query/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SV          } from '../modules/nf-core/ensemblvep/vep/main'
include { BCFTOOLS_VIEW as RESEARCH_FILTERING_SV   } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as CLINICAL_FILTERING_SV   } from '../modules/nf-core/bcftools/view/main'

//
// MODULE: Local modules
//


//
// SUBWORKFLOWS
//

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oncorefiner_pipeline'
include { PREPARE_REFERENCES     } from '../subworkflows/local/prepare_references'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOREFINER {

    take:
        ch_samplesheet  // channel: [mandatory] samplesheet read in from --input
        ch_snv_vcf      // channel: [optional]  path to input SNV vcf file
        ch_snv_vcf_tbi  // channel: [optional]  path to input SNV vcf tbi file
        ch_sv_vcf       // channel: [optional]  path to input SV vcf file
        ch_sv_vcf_tbi   // channel: [optional]  path to input SV vcf tbi file
        ch_genome_fasta // channel: [optional]  path to genome fasta file
        ch_vep_cache_unprocessed //channel: [optional] path to vep cache tar gzip file

    main:

        // Initialize input channels
        ch_versions             = channel.empty()
        ch_multiqc_files        = channel.empty()




        PREPARE_REFERENCES (
            ch_vep_cache_unprocessed
        )
        .set { ch_references }

        // Gather or get from params
        ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                            : ( params.vep_cache    ? channel.fromPath(params.vep_cache).collect() : channel.value([]) )

        //
        // Read and store paths in the vep_plugin_files file
        //
        ch_vep_extra_files_unsplit  = params.vep_plugin_files ? channel.fromPath(params.vep_plugin_files).collect() : channel.value([])
        ch_vep_extra_files = channel.empty()
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

        // Vcfanno
        ch_vcfanno_resources        = params.vcfanno_resources                  ? channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                                : channel.value([])
        ch_vcfanno_lua              = params.vcfanno_lua                        ? channel.fromPath(params.vcfanno_lua).collect()
                                                                                : channel.value([])
        ch_vcfanno_toml             = params.vcfanno_toml                       ? channel.fromPath(params.vcfanno_toml).collect()
                                                                                : channel.value([])
        ch_vcfanno_extra            = params.vcfanno_extra                      ? channel.fromPath(params.vcfanno_extra).collect()
                                                                                : []

        // SVDB
        ch_sv_dbs                   = params.svdb_query_dbs                  ? channel.fromPath(params.svdb_query_dbs)
                                                                            : channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE SNVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
        // Process SNV VCF files
        if (ch_snv_vcf) {

            // Vcfanno
            ch_snv_vcf
                .join(ch_snv_vcf_tbi)
                .map { meta, vcf, tbi ->
                    def resources = ch_vcfanno_extra
                    tuple(meta, vcf, tbi, resources)
                    }
                .set { ch_vcfanno_in }
            VCFANNO (ch_vcfanno_in, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)


            // Quality Filtering
            VCFANNO.out.vcf
                .join(VCFANNO.out.tbi)
                .map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi)
                    }
                .set { ch_research_filtering_in }
            RESEARCH_FILTERING(ch_research_filtering_in, [], [], [])


            // VEP
            RESEARCH_FILTERING.out.vcf
                    .map { meta, vcf ->
                        tuple(meta, vcf, [])
                    }
                    .set { ch_vep_snv }

            ENSEMBLVEP_SNV (
                ch_vep_snv,
                params.genome,
                params.species,
                params.vep_cache_version,
                ch_vep_cache,
                ch_genome_fasta,
                ch_vep_extra_files
            )

            // Clinical Filtering
            ENSEMBLVEP_SNV.out.vcf
                .join(ENSEMBLVEP_SNV.out.tbi)
                .map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi)
                    }
                .set { ch_clinical_filtering_in }
            CLINICAL_FILTERING(ch_clinical_filtering_in, [], [], [])


        }

        // Process SV VCF files
        if (ch_sv_vcf) {

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
                params.genome,
                params.species,
                params.vep_cache_version,
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
