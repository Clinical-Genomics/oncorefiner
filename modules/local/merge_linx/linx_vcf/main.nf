process MERGE_LINX_VCF {
    tag "$meta.id"

    input:
    tuple val(meta), path(merged_tsv), path(vcfentries_tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    """
    python3 ${params.combine_linx_vcf_script} \\
    ${merged_tsv} \\
    ${vcfentries_tsv} \\
    ${meta.id}_linx_pos.tsv
    """

    stub:
    """
    touch ${meta.id}_linx_pos.tsv
    """
}
