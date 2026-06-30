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

workflow ANNOTATE_LINX {

    take:
    ch_linx_breakends_tsv // channel: [optional]  [val(meta), path(tsv)]
    ch_linx_fusion_tsv    // channel: [optional]  [val(meta), path(tsv)]
    ch_linx_sv_tsv        // channel: [optional]  [val(meta), path(tsv)]
    ch_sv_header          // channel: [optional]  [path(txt)]
    ch_sv_vcf             // channel: [required]  [val(meta), path(vcf)]
    ch_sv_vcf_tbi         // channel: [required]  [val(meta), path(vcf.tbi)]

    main:

    ch_linx_input = ch_linx_fusion_tsv
        .join(ch_linx_breakends_tsv, failOnMismatch: true)
        .join(ch_linx_sv_tsv, failOnMismatch: true)

    COMBINE_LINX_TSV(ch_linx_input)

    ch_annotate_vcf_input = ch_sv_vcf
        .join(COMBINE_LINX_TSV.out.tsv, failOnMismatch: true)
        .combine(ch_sv_header)
        .map { meta, vcf_file, tsv_file, header_file ->
            tuple(meta, vcf_file, tsv_file, header_file)
        }


    ANNOTATE_VCF_BY_ID(ch_annotate_vcf_input)

    emit:
    vcf     = ANNOTATE_VCF_BY_ID.out.vcf
    tbi     = ANNOTATE_VCF_BY_ID.out.tbi

}
