/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

//
// MODULE: Local modules
//

include { COMBINE_LINX_TSV   } from '../../../modules/local/combine_linx_tsv/main'
include { ANNOTATE_VCF_BY_ID } from '../../../modules/local/annotate_vcf_by_id/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ANNOTATE_LINX WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_ANNOTATE_LINX_FUSIONS {

    take:
    ch_linx_breakends_tsv // channel: [required]  [val(meta), path(tsv)]
    ch_linx_fusion_tsv    // channel: [required]  [val(meta), path(tsv)]
    ch_linx_sv_tsv        // channel: [required]  [val(meta), path(tsv)]
    ch_sv_header          // channel: [required]  [val(meta), path(txt)]
    ch_sv_vcf             // channel: [required]  [val(meta), path(vcf)]
    ch_sv_vcf_tbi         // channel: [required]  [val(meta), path(vcf.tbi)]

    main:

    def ch_combine_linx_tsv_in = ch_linx_fusion_tsv
        .join(ch_linx_breakends_tsv, failOnMismatch: true)
        .join(ch_linx_sv_tsv, failOnMismatch: true)

    COMBINE_LINX_TSV(ch_combine_linx_tsv_in)

    ch_sv_vcf.view {
    "SV VCF >>> ${it}"
    }

    ch_sv_header.view {
    "SV HEADER >>> ${it}"
    }
    
    COMBINE_LINX_TSV.out.tsv.view {
    "COMBINED LINX >>> ${it}"
    }

    def ch_annotate_vcf_by_id_in = ch_sv_vcf
        .join(COMBINE_LINX_TSV.out.tsv, failOnMismatch: true)
        .join(ch_sv_header, failOnMismatch: true)

    ANNOTATE_VCF_BY_ID(ch_annotate_vcf_by_id_in)

    emit:
    vcf     = ANNOTATE_VCF_BY_ID.out.vcf
    tbi     = ANNOTATE_VCF_BY_ID.out.tbi

}
