process CNV_REPORT {
    tag "${meta.id}"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-dt_r-readr_r-rmarkdown:8612e5aa5384f180' :
        'community.wave.seqera.io/library/r-base_r-dt_r-readr_r-rmarkdown:5c3b8868d1024810' }"

    input:
    tuple val(meta), path(rmd_template), path(cnv_gene), path(cnv_segment)

    output:
    tuple val(meta), path("*.html")             , emit: report
    path  "versions.yml"                        , emit: versions
    tuple val("${task.process}"), val('rmarkdown'), eval("Rscript -e 'cat(paste(packageVersion('rmarkdown'), collapse='.'))'"), topic: versions, emit: versions_rmarkdown

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_cnv_report"
    """
    Rscript -e 'rmarkdown::render("${rmd_template}", \\
        params = list( \\
            cnv_gene = "${cnv_gene}", \\
            cnv_segment = "${cnv_segment}" \\
        ), \\
        output_file = "${prefix}.html" \\
    )' 
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_cnv_report"
    """
    touch ${prefix}.html

}
