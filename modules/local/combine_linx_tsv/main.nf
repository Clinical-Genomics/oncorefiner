process COMBINE_LINX_TSV{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89a8dce6193522bfbdcab9bc0c4e652a5533ce1c73ad3378033bd0b4e109f6f9/data' :
        'community.wave.seqera.io/library/click_python_pip_pandas:1ffa4e2ace0e5bf3' }"

    input:
    tuple val(meta), path(fusion_tsv_file), path(breakends_tsv_file), path(sv_tsv_file)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('combine_linx_tsv'), val('1.0'), topic: versions, emit: versions_combine_linx_tsv
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combine_linx_tsv.py \\
    -f ${fusion_tsv_file} \\
    -b ${breakends_tsv_file} \\
    -sv ${sv_tsv_file} \\
    -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
