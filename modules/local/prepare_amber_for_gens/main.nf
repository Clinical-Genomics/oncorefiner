process PREPARE_AMBER_FOR_GENS {
    tag "$meta.id"
    tag "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68a2ae938ce2d285498e0a81233acf6bd423ac910c3fe6fc68f44322559b1187/data' :
    'community.wave.seqera.io/library/tabix_click_pandas_py-bgzip_python:4d1c106d561e2b25' }"


    input:
    tuple val(meta), path(amber_baf_tsv)
    val(sampletype)

    output:
    tuple val(meta), path("*tumor.baf.zoom.tsv.gz"), emit: tumor_tsv
    tuple val(meta), path("*normal.baf.zoom.tsv.gz"), emit: normal_tsv, optional: true
    tuple val(meta), path("*tumor.baf.zoom.tsv.gz.tbi"), emit: tumor_tbi
    tuple val(meta), path("*normal.baf.zoom.tsv.gz.tbi"), emit: normal_tbi, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d; s/^.*bgzip //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('tabix'), eval("tabix --version | sed '1!d; s/^.*tabix //'"), topic: versions, emit: versions_tabix
    tuple val("${task.process}"), val('prepare_amber_for_gens'), eval("prepare_amber_for_gens.py --version | sed '1!d; s/^.*prepare_amber_for_gens //'"), topic: versions, emit: versions_amber_for_gens

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    prepare_amber_for_gens.py \\
        --input-file ${amber_baf_tsv} \\
        --sample-type ${sampletype} \\
        --output-file-prefix ${prefix}.${sampletype}

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
    sampletype=${sampletype}
    if [[ $sampletype == *"normal"* ]]; then
        echo | gzip > ${prefix}.${sampletype}.normal.baf.zoom.tsv.gz
        touch ${prefix}.${sampletype}.normal.baf.zoom.tsv.gz.tbi
    fi
    echo | gzip > ${prefix}.${sampletype}.tumor.baf.zoom.tsv.gz
    touch ${prefix}.${sampletype}.tumor.baf.zoom.tsv.gz.tbi
    """
}
