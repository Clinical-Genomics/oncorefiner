include { GENMOD_SCORE } from '../../../modules/nf-core/genmod/score/main'

workflow RANK_VARIANTS {

    take:
    ch_vcf           // channel: [mandatory] [ val(meta), [ path(vcf) ] ]
    ch_pedigree_file // channel: [optional]  [ path(ped) ]
    ch_score_config  // channel: [mandatory] [ path(ini) ]


    main:

    ch_genmod_score_in = ch_vcf
        .combine(ch_pedigree_file.collect(), ch_score_config.collect())

    GENMOD_SCORE(ch_genmod_score_in)

    emit:
    ranked_vcf = GENMOD_SCORE.out.vcf // channel: [ val(meta), [ path(vcf) ] ]
}
