process PREPARE_COBALT_FOR_GENS {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68a2ae938ce2d285498e0a81233acf6bd423ac910c3fe6fc68f44322559b1187/data' :
    'community.wave.seqera.io/library/tabix_click_pandas_py-bgzip_python:4d1c106d561e2b25' }"

    input:
    tuple val(meta), path(cobalt_pcf)

    output:
    tuple val(meta), path("*.bed.gz"),  emit: bed
    tuple val(meta), path("*.bed.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d; s/^.*bgzip //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('tabix'), eval("tabix --version | sed '1!d; s/^.*tabix //'"), topic: versions, emit: versions_tabix
    tuple val("${task.process}"), val('prepare_cobalt_for_gens'), eval("prepare_cobalt_for_gens.py --version | sed '1!d; s/^.*prepare_cobalt_for_gens //'"), topic: versions, emit: versions_cobalt_for_gens

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    prepare_cobalt_for_gens.py \\
        --input-file ${cobalt_pcf} \\
        --output-file ${prefix}.bed

    bgzip ${prefix}.bed

    tabix ${prefix}.bed.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip >  ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    """

}
