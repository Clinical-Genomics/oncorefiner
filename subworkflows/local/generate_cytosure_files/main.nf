//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { BCFTOOLS_VIEW as FILTER_VCF      } from '../../../modules/nf-core/bcftools/view/main'
include { VCF2CYTOSURE                     } from '../../../modules/nf-core/vcf2cytosure/main'


workflow GENERATE_CYTOSURE_FILES {
    take:
        ch_vcf // channel: [val(meta), path(vcf)]
        ch_tbi // channel: [val(meta), path(tbi)]

    main:

        // Filter out SNG variants, if present
        ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch: true)
        FILTER_VCF(
            ch_vcf_tbi,
            [],
            [],
            []
        )

        // vcf2cytosure
        VCF2CYTOSURE (
            FILTER_VCF.out.vcf,
            [[],[]],
            [[],[]],
            [[],[]],
            []
        )

    emit:
        ch_cgh_files = VCF2CYTOSURE.out.cgh // channel: [val(meta), path(cgh)]

}
