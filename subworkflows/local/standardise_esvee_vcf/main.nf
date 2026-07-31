//
// Standardise ESVEE structural variant records and create a bgzipped, tabix-indexed VCF
//

include { STANDARDISE_ESVEE_RECORDS } from '../../../modules/local/standardise_esvee_records/main'
include { HTSLIB_BGZIPTABIX         } from '../../../modules/nf-core/htslib/bgziptabix/main'


workflow STANDARDISE_ESVEE_VCF {
    take:
    ch_sv_vcf // channel: [mandatory] [val(meta), path(vcf)]

    main:
    // Convert paired ESVEE breakend records into symbolic SV records
    STANDARDISE_ESVEE_RECORDS(
        ch_sv_vcf
    )

    // Compress and tabix-index the standardised VCF
    HTSLIB_BGZIPTABIX(
        STANDARDISE_ESVEE_RECORDS.out.vcf,
        [],
        [],
        'compress',
        true,
        'vcf'
    )

    emit:
    vcf    = HTSLIB_BGZIPTABIX.out.output          // channel: [val(meta), path(vcf.gz)]
    tbi    = HTSLIB_BGZIPTABIX.out.index           // channel: [val(meta), path(vcf.gz.tbi)]
    report = STANDARDISE_ESVEE_RECORDS.out.report  // channel: [val(meta), path(report.tsv)]
}