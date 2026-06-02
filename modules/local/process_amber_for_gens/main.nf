

process process_amber_for_gens {
    tag "$meta.id"

    input:
    tuple val(meta), path(amber_baf_tsv)

    output:
    tuple val(meta), path("*.baf.zoom.tsv.gz"), path("*.baf.zoom.tsv.gz.tbi"), emit: baf_zoom

    script:
    def prefix = meta.id ?: "sample"

    """
    amber_baf_for_gens.py \\
        --input-file ${amber_baf_tsv} \\
        --tumor-output ${prefix}.tumor.baf.zoom.tsv \\
        --normal-output ${prefix}.normal.baf.zoom.tsv

    for f in *.baf.zoom.tsv; do
        if [[ -s "\$f" ]]; then
            bgzip -f "\$f"
            tabix -f -s 1 -b 2 -e 3 "\$f.gz"
        fi
    done
    """
}
