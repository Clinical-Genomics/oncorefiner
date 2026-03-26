//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'


workflow PREPARE_REFERENCES {
    take:
        val_vep_cache // channel: [mandatory] [ path(cache) ]

    main:
        vep_resources = channel.empty()

        //
        // Prepare vep cache files
        //

        if (val_vep_cache) {
            if (val_vep_cache.endsWith("tar.gz")) {
                ch_vep_resources = UNTAR_VEP_CACHE (channel.fromPath(val_vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()).untar.map{ _meta, files -> [files]}.collect()
            } else {
                ch_vep_resources = channel.fromPath(val_vep_cache).collect()
            }
        }

    emit:
        vep_resources
}
