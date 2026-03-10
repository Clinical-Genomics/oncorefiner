//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'


workflow PREPARE_REFERENCES {
    take:
        ch_vep_cache                 // channel: [mandatory for annotation] [ path(cache) ]

    main:

        // Untar
        UNTAR_VEP_CACHE (ch_vep_cache)

    emit:
        vep_resources         = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()
}
