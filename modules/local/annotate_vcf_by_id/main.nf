process ANNOTATE_VCF_BY_ID{
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a1a02b55777ce514673dfcc42280441fbbdcae21cbd89d956eb6459a0117f5e4/data' :
        'community.wave.seqera.io/library/pysam_click_python:573fe89e5d35db27' }"

    input:
    tuple val(meta), path(vcf_file), path(tsv_file), path(header_file)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('annotate_vcf_by_id'), topic: versions, emit: versions_annotate_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -Xgil=0 ${moduleDir}/scripts/annotate_vcf_by_id.py \\
    -v ${vcf_file} \\
    -h ${header_file} \\
    -t ${tsv_file} \\
    -o ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
