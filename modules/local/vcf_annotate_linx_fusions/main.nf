process VCF_ANNOTATE_LINX_FUSIONS{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/62acd8525fffdf87492840e9b6f332cf17ef525057a414601b67ebd82604f9a0/data' :
        'community.wave.seqera.io/library/pysam_click_pandas_python:dd208bd749e55b89' }"
    input:
    tuple val(meta), path(vcf_file), path(tsv_file), path(header_file)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('vcf_annotate_linx_fusions'), eval("vcf_annotate_linx_fusions.py --version | sed '1!d; s/^.*vcf_annotate_linx_fusions //'"), topic: versions, emit: versions_annotate_vcf
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_annotate_linx_fusions.py \\
    -v ${vcf_file} \\
    -h ${header_file} \\
    -t ${tsv_file} \\
    -o ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip >  ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
