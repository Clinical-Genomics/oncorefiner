process LINX_TSV{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam_click:b17ee92475905c52' :
        'pysam_click:b17ee92475905c52' }"

    input:
    tuple val(meta), path(fusion_tsv_file), path(breakends_tsv_file), path(sv_tsv_file), path(output_file)

    output:
    tuple val(meta), path("*tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ./scripts/combine_linx_tsv.py \\
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
