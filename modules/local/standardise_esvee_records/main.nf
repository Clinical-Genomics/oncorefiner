process STANDARDISE_ESVEE_RECORDS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02f777287b58f5abc9a3748d6d33d83d52e83e3afdf424b1b9009a5511789995/data' :
        'community.wave.seqera.io/library/python_click_pysam:5794ef2a03d1ee5e' }"

    input:
    tuple val(meta), path(sv_vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.report.tsv"), emit: report
    tuple val("${task.process}"),
        val('standardise_esvee_records'),
        // keep this in sync with the version in `bin/standardise_esvee_records.py`
        val('1.0.0'),
        topic: versions,
        emit: versions_standardise_esvee_records

    script:
    def prefix = task.ext.prefix ?: meta.id

    """
    standardise_esvee_records.py \\
        ${sv_vcf} \\
        --output ${prefix}.vcf \\
        --report ${prefix}.report.tsv \\
        --verbose-report
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    touch ${prefix}.vcf
    touch ${prefix}.report.tsv
    """
}