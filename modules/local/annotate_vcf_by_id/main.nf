process ANNOTATE_VCF_BY_ID{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ac/ac7bb4d67fee1c5348a574bb9a4fd5b774627a6c8a2671fab65e5f3577b342a7/data' :
        'community.wave.seqera.io/library/pysam_click_python:90b6c5beb8454821' }"
    input:
    tuple val(meta), path(vcf_file), path(tsv_file), path(header_file)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('annotate_vcf_by_id'), val('1.0'), topic: versions, emit: versions_annotate_vcf
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_vcf_by_id.py \\
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
