//
// Standardise ESVEE structural variant records and create a bgzipped, tabix-indexed VCF
//

include { STANDARDISE_ESVEE_RECORDS } from '../../../modules/local/standardise_esvee_records/main'
include { TABIX_BGZIP              } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX              } from '../../../modules/nf-core/tabix/tabix/main'


workflow STANDARDISE_ESVEE_VCF {
    take:
    ch_sv_vcf // channel: [mandatory] [val(meta), path(vcf)]

    main:
    // Convert paired ESVEE breakend records into symbolic SV records
    STANDARDISE_ESVEE_RECORDS(
        ch_sv_vcf
    )

    // Compress the standardised VCF using bgzip
    TABIX_BGZIP(
        STANDARDISE_ESVEE_RECORDS.out.vcf
    )

    // Create a tabix index for the compressed VCF
    TABIX_TABIX(
        TABIX_BGZIP.out.gz
    )

    emit:
    vcf    = TABIX_BGZIP.out.gz                    // channel: [val(meta), path(vcf.gz)]
    tbi    = TABIX_TABIX.out.index                 // channel: [val(meta), path(vcf.gz.tbi)]
    report = STANDARDISE_ESVEE_RECORDS.out.report  // channel: [val(meta), path(report.tsv)]
}