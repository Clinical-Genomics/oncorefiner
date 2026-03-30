//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { BCFTOOLS_VIEW as FILTER_VCF      } from '../../../modules/nf-core/bcftools/view/main'
include { TIDDIT_COV                       } from '../../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE                     } from '../../../modules/nf-core/vcf2cytosure/main'


workflow GENERATE_CYTOSURE_FILES {
    take:
        ch_vcf // channel: [val(meta), path(vcf)]
        ch_tbi // channel: [val(meta), path(tbi)]
        ch_bam // channel: [val(meta), path(bam)]
        ch_bai // channel: [val(meta), path(bai)]

    main:

        // Filter out SNG variants, if present
        ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch: true)
        FILTER_VCF(
            ch_vcf_tbi,
            [],
            [],
            []
        )

        // Generate coverage bed file
        ch_bam_bai = ch_bam.join(ch_bai, failOnMismatch: true)
        TIDDIT_COV (
            ch_bam_bai,
             [[],[]]
        )

        // vcf2cytosure
        VCF2CYTOSURE (
            FILTER_VCF.out.vcf,
            TIDDIT_COV.out.cov,
            [[],[]],
            [[],[]],
            []
        )

    emit:
        ch_cgh_files = VCF2CYTOSURE.out.cgh // channel: [val(meta), path(cgh)]

}
