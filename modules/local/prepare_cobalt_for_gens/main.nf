process PREPARE_COBALT_FOR_GENS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68a2ae938ce2d285498e0a81233acf6bd423ac910c3fe6fc68f44322559b1187/data' :
    'community.wave.seqera.io/library/tabix_click_pandas_py-bgzip_python:4d1c106d561e2b25' }"

    input:
    tuple val(meta), path(cobalt_pcf)

    output:
    tuple val(meta), path("*${meta.sample_type}.bed.gz"),  emit: bed
    tuple val(meta), path("*${meta.sample_type}.bed.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('prepare_cobalt_for_gens'), val('1.0'), topic: versions, emit: versions_cobalt_for_gens
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    script:
    def prefix = meta.id ?: "sample"

    """
    prepare_cobalt_for_gens.py \\
        --input-file ${cobalt_pcf} \\
        --output-file ${prefix}.${meta.sample_type}.bed

    bgzip ${prefix}.${meta.sample_type}.bed

    tabix ${prefix}.${meta.sample_type}.bed.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip >  ${prefix}.${meta.sample_type}.bed.gz
    touch ${prefix}.${meta.sample_type}.bed.gz.tbi
    """

}
