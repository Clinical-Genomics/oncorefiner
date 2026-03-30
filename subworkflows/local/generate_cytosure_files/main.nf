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
        ch_vcf            // channel: [mandatory] [val(meta), path(vcf)]
        ch_tbi            // channel: [mandatory] [val(meta), path(tbi)]
        ch_bam_bai_tumor  // channel: [mandatory] [val(meta), path(bam), path(bai)]
        ch_bam_bai_normal // channel: [optional]  [val(meta), path(bam), path(bai)]

    main:

        // Filter out SNG variants, if present
        ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch: true)
        FILTER_VCF(
            ch_vcf_tbi,
            [],
            [],
            []
        )

        // Workflow for tumor sample

        // Generate tumor coverage bed file
        TIDDIT_COV_TUMOR (
            ch_bam_bai_tumor,
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

        // Worflow for normal sample, if BAM and BAI files for the normal sample are provided
        ch_bam_bai_normal = channel.empty()

        if (ch_bam_bai_normal) {

            // Generate normal coverage bed file
            TIDDIT_COV_NORMAL (
                ch_bam_bai_normal,
                [[],[]]
            )

            // Run vcf2cytosure for normal sample
            VCF2CYTOSURE_NORMAL (
                FILTER_VCF.out.vcf,
                TIDDIT_COV_NORMAL.out.cov,
                [[],[]],
                [[],[]],
                []
            )
        }

    emit:
        ch_cgh_tumor  = VCF2CYTOSURE_TUMOR.out.cgh  // channel: [val(meta), path(cgh)]
        ch_cgh_normal = VCF2CYTOSURE_NORMAL.out.cgh // channel: [val(meta), path(cgh)]

}
