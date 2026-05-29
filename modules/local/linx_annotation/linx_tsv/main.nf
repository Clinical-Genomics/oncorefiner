process MERGE_LINX_TSV {
    tag "$meta.id"

    input:
    tuple val(meta), path(fusion_tsv), path(breakend_tsv), path(svs_tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    """
    python3 ${params.combine_linx_tsv_script} \\
    ${fusion_tsv} \\
    ${breakend_tsv} \\
    ${svs_tsv} \\
    ${meta.id}_merged.tsv
    """

    stub:
    """
    touch ${meta.id}_merged.tsv
    """
}
