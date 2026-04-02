


workflow process_snvs {

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE SNVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
        // Process SNV VCF files
        if (ch_snv_vcf) {

            // Vcfanno
            ch_snv_vcf
                .join(ch_snv_vcf_tbi)
                .map { meta, vcf, tbi ->
                    def resources = ch_vcfanno_extra
                    tuple(meta, vcf, tbi, resources)
                    }
                .set { ch_vcfanno_in }
            VCFANNO (ch_vcfanno_in, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)


            // Quality Filtering
            VCFANNO.out.vcf
                .join(VCFANNO.out.tbi)
                .map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi)
                    }
                .set { ch_research_filtering_in }
            RESEARCH_FILTERING(ch_research_filtering_in, [], [], [])


            // VEP
            RESEARCH_FILTERING.out.vcf
                    .map { meta, vcf ->
                        tuple(meta, vcf, [])
                    }
                    .set { ch_vep_snv }

            ENSEMBLVEP_SNV (
                ch_vep_snv,
                val_genome,
                val_species,
                val_vep_cache_version,
                ch_vep_cache,
                ch_genome_fasta,
                ch_vep_extra_files
            )

            // Clinical Filtering
            ENSEMBLVEP_SNV.out.vcf
                .join(ENSEMBLVEP_SNV.out.tbi)
                .map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi)
                    }
                .set { ch_clinical_filtering_in }
            CLINICAL_FILTERING(ch_clinical_filtering_in, [], [], [])


        }
}
