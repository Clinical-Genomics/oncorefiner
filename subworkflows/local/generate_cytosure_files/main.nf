//
// Convert VCF with structural variations to the .cgh format used by the CytoSure Interpret Software
//

include { BCFTOOLS_VIEW as FILTER_VCF } from '../../../modules/nf-core/bcftools/view/main'
include { TIDDIT_COV                  } from '../../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE                } from '../../../modules/nf-core/vcf2cytosure/main'


workflow GENERATE_CYTOSURE_FILES {
    take:
    ch_bam_bai // channel: [optional]  [val(meta), path(bam), path(bai)]
    ch_tbi     // channel: [mandatory] [val(meta), path(tbi)]
    ch_vcf     // channel: [mandatory] [val(meta), path(vcf)]

    main:
    // Filter out SGL variants, if present
    ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch: true)
    FILTER_VCF(
        ch_vcf_tbi,
        [],
        [],
        []
    )

    // Generate coverage bed file
    TIDDIT_COV(
        ch_bam_bai,
         [[],[]]
    )

    // Run vcf2cytosure
    VCF2CYTOSURE(
        FILTER_VCF.out.vcf,
        TIDDIT_COV.out.cov,
        [[],[]],
        [[],[]],
        []
    )
}
