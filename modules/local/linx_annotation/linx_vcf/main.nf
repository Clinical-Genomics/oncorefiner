process LINX_VCF{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'container full path' :
        'short path' }"

    input:
    tuple val(meta), path(vcf_file), path(header_file), path(tsv_file), path(output_file)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python -Xgil=0 ./scripts/annotate_vcf_linx.py \\
    -v ${vcf_file} \\
    -h ${header_file} \\
    -t  ${tsv_file} \\
    -o ${output_file}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${output_file}
    """
}
