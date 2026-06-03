

process process_cobolt_for_gens {
    tag "$meta.id"

    input:
    tuple val(meta), path(cobolt_pcf), val(sample_type)

    output:
    tuple val(meta), val(sample_type), path("*.cobolt.cov.bed.gz"), path("*.cobolt.cov.bed.gz.tbi"), emit: cobolt_cov

    script:
    def prefix = meta.id ?: "sample"

    """
    cobolt_pcf_for_gens.py \\
        --input-file ${cobolt_pcf} \\
        --output-file ${prefix}.${sample_type}.cobolt.cov.bed

    bgzip ${prefix}.${sample_type}.cobolt.cov.bed

    tabix ${prefix}.${sample_type}.cobolt.cov.bed.gz
    """
}
