//
// Prepare reference files
//

include { UNTAR as UNTAR_VEP_CACHE           } from '../../modules/nf-core/untar/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_NML } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_TMR } from '../../modules/nf-core/samtools/view/main'


workflow PREPARE_REFERENCES {
    take:
        ch_bam_bai_normal     // channel: [optional]  channel containing BAM and BAI file for normal sample
        ch_bam_bai_tumor      // channel: [optional]  channel containing BAM and BAI file for tumor sample
        ch_genome_fasta_fai   // channel: [optional]  channel containing genome fasta and fai files
        val_vep_cache         // string:  [mandatory] path(cache)

    main:
        ch_vep_resources = channel.value([[]])
        ch_cram_crai_normal = channel.empty()
        ch_cram_crai_tumor = channel.empty()

        //
        // Convert BAM to CRAM
        //
        if (ch_bam_bai_normal) {
            SAMTOOLS_VIEW_NML ( ch_bam_bai_normal, ch_genome_fasta_fai, [[], []], [[],[]], 'crai' )

            SAMTOOLS_VIEW_NML.out.cram
                .join(SAMTOOLS_VIEW_NML.out.crai)
                .set { ch_cram_crai_normal }
        }

        if (ch_bam_bai_tumor) {
            SAMTOOLS_VIEW_TMR ( ch_bam_bai_tumor, ch_genome_fasta_fai, [[], []], [[],[]], 'crai' )

            SAMTOOLS_VIEW_TMR.out.cram
                .join(SAMTOOLS_VIEW_TMR.out.crai)
                .set { ch_cram_crai_tumor }
        }


        //
        // Prepare vep cache files
        //

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
        vep_resources = ch_vep_resources    // channel: [vep_cache_files]
        ch_cram_crai_normal                 // channel: [cram, crai] for normal sample
        ch_cram_crai_tumor                  // channel: [cram, crai] for tumor sample
}
