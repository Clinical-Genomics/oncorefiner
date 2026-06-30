

process PREPARE_AMBER_FOR_GENS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68a2ae938ce2d285498e0a81233acf6bd423ac910c3fe6fc68f44322559b1187/data' :
    'community.wave.seqera.io/library/tabix_click_pandas_py-bgzip_python:4d1c106d561e2b25' }"


    input:
    tuple val(meta), path(amber_baf_tsv)

    output:
    tuple val(meta), path("*.baf.zoom.tsv.gz"), emit: tsv
    tuple val(meta), path("*.baf.zoom.tsv.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('prepare_amber_for_gens'), val('1.0'), topic: versions, emit: versions_amber_for_gens
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    script:
    def prefix = meta.id ?: "sample"

    """
    amber_baf_for_gens.py \\
        --input-file ${amber_baf_tsv} \\
        --sample_type ${meta.sample_type} \\
        --output-file ${prefix}.${meta.sample_type}.baf.zoom.tsv

    for f in *.baf.zoom.tsv; do
        if [[ -s "\$f" ]]; then
            bgzip -f "\$f"
            tabix -f -s 1 -b 2 -e 3 "\$f.gz"
        fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip >  ${prefix}.${meta.sample_type}.bed.gz
    touch ${prefix}.${meta.sample_type}.bed.gz.tbi
    """
}
