//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE                           } from '../../../modules/nf-core/untar/main'


workflow PREPARE_REFERENCES {
    take:
        val_bam_normal // string: [optional]  path(bam)
        val_bai_normal // string: [optional]  path(bai)
        val_bam_tumor  // string: [optional]  path(bam)
        val_bai_tumor  // string: [optional]  path(bai)
        val_vep_cache  // string: [mandatory] path(cache)

    main:
        //
        // Prepare bam_bai channels for tumor and normal samples
        //
        ch_bam_bai_tumor = channel.empty()
        if (val_bam_tumor && val_bai_tumor) {
            ch_bam_bai_tumor = channel.fromPath(val_bam_tumor)
                                .join(channel.fromPath(val_bai_tumor))
                                .map { it -> [[id:'tumor'], it] }
        }

        ch_bam_bai_normal = channel.empty()
        if (val_bam_normal && val_bai_normal) {
            ch_bam_bai_normal = channel.fromPath(val_bam_normal)
                                .join(channel.fromPath(val_bai_normal))
                                .map { it -> [[id:'normal'], it] }
        }

        ch_bam_bai_tumor.dump(tag: "ch_bam_bai_tumor")
        ch_bam_bai_normal.dump(tag: "ch_bam_bai_normal")

        //
        // Prepare vep cache files
        //
        ch_vep_resources = channel.value([[]])

        if (val_vep_cache) {
            if (val_vep_cache.endsWith("tar.gz")) {
                ch_vep_resources = UNTAR_VEP_CACHE(
                                    channel.fromPath(val_vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                    )
                                    .untar.map{ _meta, files -> [files]}
                                    .collect()
            } else {
                ch_vep_resources = channel.fromPath(val_vep_cache).collect()
            }
        }

    emit:
        bam_bai_tumor  = ch_bam_bai_tumor  // channel: [val(meta), path(bam), path(bai)]
        bam_bai_normal = ch_bam_bai_normal // channel: [val(meta), path(bam), path(bai)]
        vep_resources  = ch_vep_resources  // channel: [vep_cache_files]
}
