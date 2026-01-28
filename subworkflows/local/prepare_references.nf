//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'


workflow PREPARE_REFERENCES {
    take:
        ch_vep_cache                 // channel: [mandatory for annotation] [ path(cache) ]


    main:
        ch_versions      = channel.empty()

        // Untar
        UNTAR_VEP_CACHE (ch_vep_cache)

        // Gather versions
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        vep_resources         = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()
        versions              = ch_versions
}
