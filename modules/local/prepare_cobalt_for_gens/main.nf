process PREPARE_COBALT_FOR_GENS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89a8dce6193522bfbdcab9bc0c4e652a5533ce1c73ad3378033bd0b4e109f6f9/data'
        : 'community.wave.seqera.io/library/click_python_pip_pandas:1ffa4e2ace0e5bf3' }"

    input:
    tuple val(meta), path(cobalt_pcf), val(sample_type)

    output:
    tuple val(meta), val(sample_type), path("*.cobalt.cov.bed.gz"),  emit: gens_cov_bed
    tuple val(meta), val(sample_type), path("*.cobalt.cov.bed.gz.tbi"), emit: gens_cov_tbi
    tuple val("${task.process}"), val('prepare_cobalt_for_gens'), val('1.0'), topic: versions, emit: versions_cobalt_for_gens
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    script:
    def prefix = meta.id ?: "sample"

    """
    prepare_cobalt_for_gens.py \\
        --input-file ${cobalt_pcf} \\
        --output-file ${prefix}.${sample_type}.cobalt.cov.bed

    bgzip ${prefix}.${sample_type}.cobalt.cov.bed

    tabix ${prefix}.${sample_type}.cobalt.cov.bed.gz
    """
}
