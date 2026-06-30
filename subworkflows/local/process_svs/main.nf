/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_RESEARCH   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLINICAL   } from '../../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP_VEP                           } from '../../../modules/nf-core/ensemblvep/vep/main'
include { SVDB_QUERY                               } from '../../../modules/nf-core/svdb/query/main'

//
// LOCAL SUBWORKFLOWS
//

include { GENERATE_CYTOSURE_FILES } from '../../../subworkflows/local/generate_cytosure_files/main'
include { ANNOTATE_LINX           } from '../../../subworkflows/local/annotate_linx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_SVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_SVS {

    take:
        ch_bam_bai_normal     // channel: [optional]  [val(meta), path(bam), path(bai)]
        ch_bam_bai_tumor      // channel: [mandatory]  [val(meta), path(bam), path(bai)]
        ch_linx_breakends_tsv // channel: [optional]  [val(meta), path(tsv)]
        ch_linx_fusion_tsv    // channel: [optional]  [val(meta), path(tsv)]
        ch_linx_sv_tsv        // channel: [optional]  [val(meta), path(tsv)]
        ch_sv_header          // channel: [optional]  [path(txt)]
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
    // Annotate VCF with LINX TSVs
    ANNOTATE_LINX(
        ch_linx_breakends_tsv,
        ch_linx_fusion_tsv,
        ch_linx_sv_tsv,
        ch_sv_header,
        ch_sv_vcf,
        ch_sv_vcf_tbi
    )

    ch_sv_linx_vcf = ANNOTATE_LINX.out.vcf
    ch_sv_linx_vcf_tbi = ANNOTATE_LINX.out.tbi


    // SVDB QUERY
    ch_sv_dbs
        .multiMap { filename, in_freq_info_key, in_allele_count_info_key, out_freq_info_key, out_allele_count_info_key ->
            vcf_dbs: filename
            in_frqs: in_freq_info_key
            in_occs: in_allele_count_info_key
            out_frqs: out_freq_info_key
            out_occs: out_allele_count_info_key
        }
        .set { ch_svdb_dbs }

    SVDB_QUERY (
        ch_sv_linx_vcf,
        ch_svdb_dbs.in_occs.toList(),
        ch_svdb_dbs.in_frqs.toList(),
        ch_svdb_dbs.out_occs.toList(),
        ch_svdb_dbs.out_frqs.toList(),
        ch_svdb_dbs.vcf_dbs.toList(),
        []
    )

    // Quality Filtering
    SVDB_QUERY.out.vcf
        .map { meta, vcf ->
            tuple(meta, vcf, []) }
        .set { ch_research_filtering_sv_in }


    BCFTOOLS_VIEW_RESEARCH(ch_research_filtering_sv_in, [], [], [])

    // VEP
    BCFTOOLS_VIEW_RESEARCH.out.vcf
        .map { meta, vcf ->
                    tuple(meta, vcf, []) }
        .set { ch_vep_sv }


    ENSEMBLVEP_VEP(
        ch_vep_sv,
        val_genome,
        val_species,
        val_vep_cache_version,
        ch_vep_cache,
        ch_genome_fasta,
        ch_vep_extra_files
    )

    // Clinical Filtering
    ENSEMBLVEP_VEP.out.vcf
        .join(ENSEMBLVEP_VEP.out.tbi)
        .map { meta, vcf, tbi ->
            tuple(meta, vcf, tbi)
            }
        .set { ch_clinical_filtering_sv_in }
    BCFTOOLS_VIEW_CLINICAL(ch_clinical_filtering_sv_in, [], [], [])

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
}
