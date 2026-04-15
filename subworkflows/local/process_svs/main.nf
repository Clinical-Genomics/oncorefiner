/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { SVDB_QUERY as SVDB_QUERY_DB              } from '../../../modules/nf-core/svdb/query/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SV          } from '../../../modules/nf-core/ensemblvep/vep/main'
include { BCFTOOLS_VIEW as RESEARCH_FILTERING_SV   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as CLINICAL_FILTERING_SV   } from '../../../modules/nf-core/bcftools/view/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_SVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_SVS {

    take:
        ch_sv_vcf             // channel: [required]  [val(meta), path(vcf)]
        ch_sv_vcf_tbi         // channel: [required]  [val(meta), path(vcf.tbi)]
        ch_sv_dbs             // channel: [required]  path(svdb_dbs_csv)
        val_genome            // value:   [required]  Genome build (e.g. GRCh38)
        val_species           // value:   [required]  Species
        val_vep_cache_version // value:   [required]  VEP cache
        ch_vep_cache          // channel: [optional]  [val(meta), path(vep_cache)]
        ch_genome_fasta       // channel: [optional]  [val(meta), path(fasta)]
        ch_vep_extra_files    // channel: [optional]  [val(meta), path(vep_extra_files)]

    main:
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
