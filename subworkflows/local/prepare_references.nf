//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'


workflow PREPARE_REFERENCES {
    take:
        ch_vep_cache // channel: [mandatory] [ path(cache) ]

    main:
        vep_resources = channel.empty()

        //
        // Prepare vep cache files
        //

        if (ch_vep_cache) {
            if (ch_vep_cache.vep_cache.endsWith("tar.gz")) {
                ch_vep_resources = UNTAR_VEP_CACHE (ch_vep_cache.vep_cache).untar.map{ _meta, files -> [files]}.collect()
            } else {
                ch_vep_resources = channel.fromPath(ch_vep_cache.vep_cache).collect()
            }
        }



    emit:
        vep_resources
}
