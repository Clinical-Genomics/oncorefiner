process COMBINE_LINX_TSV{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a1a02b55777ce514673dfcc42280441fbbdcae21cbd89d956eb6459a0117f5e4/data' :
        'community.wave.seqera.io/library/pysam_click_python:573fe89e5d35db27' }"

    input:
    tuple val(meta), path(fusion_tsv_file), path(breakends_tsv_file), path(sv_tsv_file), path(output_file)

    output:
    tuple val(meta), path("*tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python -Xgil=0 ./scripts/combine_linx_tsv.py \\
    -f ${fusion_tsv_file} \\
    -b ${breakends_tsv_file} \\
    -sv ${sv_tsv_file} \\
    -o ${output_file}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${output_file}
    """
}
