


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { ENSEMBLVEP_VEP as ENSEMBLVEP_SNV         } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                                  } from '../../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_VIEW as RESEARCH_FILTERING      } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as CLINICAL_FILTERING      } from '../../../modules/nf-core/bcftools/view/main'


workflow PROCESS_SNVS {

    take:
        ch_genome_fasta       // channel: [optional]  [val(meta), path(fasta)]
        ch_snv_vcf            // channel: [optional]  [val(meta), path(vcf)]
        ch_snv_vcf_tbi        // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_vcfanno_extra      // channel: [optional]  [path(extra_file1), path(extra_file2), ...]
        ch_vcfanno_lua        // channel: [optional]  [path(lua_file)]
        ch_vcfanno_resources  // channel: [optional]  [path(resource_file1), path(resource_file2), ...]
        ch_vcfanno_toml       // channel: [optional]  [path(toml_file)]
        ch_vep_cache          // channel: [optional]  [vep_cache_files]
        ch_vep_extra_files    // channel: [optional]  [path(plugin_file1), path(plugin_file2), ...]
        val_genome            // string:  [optional]  genome assembly (e.g. "GRCh38")
        val_species           // string:  [optional]  species (e.g. "homo_sapiens")
        val_vep_cache_version // string:  [optional]  version of vep cache to use (e.g. "107")

    main:

        // Annotate with custom databases
        ch_snv_vcf
            .join(ch_snv_vcf_tbi)
            .map { meta, vcf, tbi ->
                def resources = ch_vcfanno_extra
                tuple(meta, vcf, tbi, resources)
                }
            .set { ch_vcfanno_in }

        VCFANNO (ch_vcfanno_in, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)

        // Research filtering
        VCFANNO.out.vcf
            .join(VCFANNO.out.tbi)
            .map { meta, vcf, tbi ->
                tuple(meta, vcf, tbi)
                }
            .set { ch_research_filtering_in }

        RESEARCH_FILTERING(ch_research_filtering_in, [], [], [])

        // Annotate with VEP
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
