//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { BCFTOOLS_VIEW as FILTER_VCF         } from '../../../modules/nf-core/bcftools/view/main'
include { TIDDIT_COV as TIDDIT_COV_TUMOR    } from '../../../modules/nf-core/tiddit/cov/main'
include { TIDDIT_COV as TIDDIT_COV_NORMAL   } from '../../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE as VCF2CYTOSURE_TUMOR  } from '../../../modules/nf-core/vcf2cytosure/main'
include { VCF2CYTOSURE as VCF2CYTOSURE_NORMAL } from '../../../modules/nf-core/vcf2cytosure/main'


workflow GENERATE_CYTOSURE_FILES {
    take:
        ch_vcf        // channel: [mandatory] [val(meta), path(vcf)]
        ch_tbi        // channel: [mandatory] [val(meta), path(tbi)]
        ch_tumor_bam  // channel: [mandatory] [val(meta), path(bam)]
        ch_tumor_bai  // channel: [mandatory] [val(meta), path(bai)]
        ch_normal_bam // channel: [optional]  [val(meta), path(bam)]
        ch_normal_bai // channel: [optional]  [val(meta), path(bai)]

    main:

        // Filter out SNG variants, if present
        ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch: true)
        FILTER_VCF(
            ch_vcf_tbi,
            [],
            [],
            []
        )

        // Generate tumor coverage bed file
        ch_tumor_bam_bai = ch_tumor_bam.join(ch_tumor_bai, failOnMismatch: true)
        TIDDIT_COV_TUMOR (
            ch_tumor_bam_bai,
             [[],[]]
        )

        // Run vcf2cytosure for tumor sample
        VCF2CYTOSURE_TUMOR (
            FILTER_VCF.out.vcf,
            TIDDIT_COV_TUMOR.out.cov,
            [[],[]],
            [[],[]],
            []
        )

        if (ch_normal_bam && ch_normal_bai) {

            // Generate normal coverage bed file
            ch_normal_bam_bai = ch_normal_bam.join(ch_normal_bai, failOnMismatch: true)
            TIDDIT_COV_NORMAL (
                ch_normal_bam_bai,
                [[],[]]
            )

            // Run vcf2cytosure for normal sample
            ch_normal_bam_bai = ch_normal_bam.join(ch_normal_bai, failOnMismatch: true)
            VCF2CYTOSURE_NORMAL (
                FILTER_VCF.out.vcf,
                TIDDIT_COV_NORMAL.out.cov,
                [[],[]],
                [[],[]],
                []
            )
        }

    emit:
        ch_tumor_cgh  = VCF2CYTOSURE_TUMOR.out.cgh  // channel: [val(meta), path(cgh)]
        ch_normal_cgh = VCF2CYTOSURE_NORMAL.out.cgh // channel: [val(meta), path(cgh)]

}
