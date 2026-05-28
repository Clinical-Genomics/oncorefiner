/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { ENSEMBLVEP_VEP                          } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                                 } from '../../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_RESEARCH } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLINICAL } from '../../../modules/nf-core/bcftools/view/main'
include { ANNOTATE_CADD                           } from '../../../subworkflows/local/annotate_cadd'

//
// SUBWORKFLOW: Installed directly from genomic-medicine-sweden/subworkflows
//

include { VCF_ANNOTATE_SCORE_GENMOD } from '../../../subworkflows/genomic-medicine-sweden/vcf_annotate_score_genmod/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PROCESS_SNVS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_SNVS {

    take:
        ch_genome_fasta          // channel: [optional]  [val(meta), path(fasta)]
        ch_genome_fai            // channel: [optional]  [val(meta), path(fai)]
        ch_cadd_header           // channel: [optional]  [val(meta), path(header_file)]
        ch_cadd_prescored_indels // channel: [optional]  [val(meta), path(dir)]
        ch_cadd_resources        // channel: [optional]  [val(meta), path(dir)]
        ch_genmod_score_config   // channel: [optional]  [val(meta), path(ini)]
        ch_snv_vcf               // channel: [optional]  [val(meta), path(vcf)]
        ch_snv_vcf_tbi           // channel: [optional]  [val(meta), path(vcf.tbi)]
        ch_vcfanno_extra         // channel: [optional]  [path(extra_file1), path(extra_file2), ...]
        ch_vcfanno_lua           // channel: [optional]  [path(lua_file)]
        ch_vcfanno_resources     // channel: [optional]  [path(resource_file1), path(resource_file2), ...]
        ch_vcfanno_toml          // channel: [optional]  [path(toml_file)]
        ch_vep_cache             // channel: [optional]  [path(vep_cache)]
        ch_vep_extra_files       // channel: [optional]  [path(plugin_file1), path(plugin_file2), ...]
        val_cadd_resources       // string:  [optional]  path to CADD resources directory
        val_genome               // string:  [optional]  genome assembly (e.g. "GRCh38")
        val_skip_genmod_score    // boolean: [mandatory] whether to skip VCF_ANNOTATE_SCORE_GENMOD process
        val_species              // string:  [optional]  species (e.g. "homo_sapiens")
        val_vep_cache_version    // string:  [optional]  version of vep cache to use (e.g. "107")

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

        BCFTOOLS_VIEW_RESEARCH(ch_research_filtering_in, [], [], [])

        // Annotate with VEP
        BCFTOOLS_VIEW_RESEARCH.out.vcf
                .map { meta, vcf ->
                    tuple(meta, vcf, [])
                }
                .set { ch_vep_snv }

        // ANNOTATE WITH CADD - currently depends on val_cadd_resources - could be improved?
        if (val_cadd_resources) {

            ch_cadd_in = BCFTOOLS_VIEW_RESEARCH.out.vcf

            ANNOTATE_CADD (
                ch_cadd_in,
                val_genome,
                ch_genome_fai,
                ch_cadd_header,
                ch_cadd_resources,
                ch_cadd_prescored_indels
            )

            ANNOTATE_CADD.out.vcf
                .join(ANNOTATE_CADD.out.tbi)
                .set { ch_vep_snv }
        }

        ENSEMBLVEP_VEP (
            ch_vep_snv,
            val_genome,
            val_species,
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        if (!val_skip_genmod_score) {
            // Rank and add score annotation with genmod score
            ch_annotate_score_genmod_in = ENSEMBLVEP_VEP.out.vcf
                .combine(ch_genmod_score_config)
                .multiMap { meta_vcf, vcf, _meta_genmod_score_config, genmod_score_config ->
                    vcf: tuple(meta_vcf, vcf)
                    score_config: tuple(meta_vcf, genmod_score_config)
                }

            VCF_ANNOTATE_SCORE_GENMOD (
                ch_annotate_score_genmod_in.vcf,
                channel.empty(),
                channel.empty(),
                ch_annotate_score_genmod_in.score_config,
                true,
            )

            ch_clinical_filtering_in = VCF_ANNOTATE_SCORE_GENMOD.out.vcf
                .join(VCF_ANNOTATE_SCORE_GENMOD.out.index)

        } else {
            ch_clinical_filtering_in = ENSEMBLVEP_VEP.out.vcf
                .join(ENSEMBLVEP_VEP.out.tbi)
        }

        // Clinical Filtering
        BCFTOOLS_VIEW_CLINICAL(ch_clinical_filtering_in, [], [], [])
}
