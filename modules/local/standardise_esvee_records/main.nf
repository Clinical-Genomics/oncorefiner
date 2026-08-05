process STANDARDISE_ESVEE_RECORDS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02f777287b58f5abc9a3748d6d33d83d52e83e3afdf424b1b9009a5511789995/data' :
        'community.wave.seqera.io/library/python_click_pysam:5794ef2a03d1ee5e' }"

    input:
    tuple val(meta), path(sv_vcf)

    output:
    tuple val(meta), path("*.standardised.vcf"), emit: vcf
    tuple val(meta), path("*.standardised.report.tsv"), emit: report
    tuple val("${task.process}"),
        val('standardise_esvee_records'),
        val('1.0'),
        topic: versions,
        emit: versions_standardise_esvee_records

    script:
    def prefix = task.ext.prefix ?: meta.id

    """
    standardise_esvee_records.py \\
        ${sv_vcf} \\
        --output ${prefix}.standardised.vcf \\
        --report ${prefix}.standardised.report.tsv \\
        --verbose-report
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    cat <<'EOF' > ${prefix}.standardised.vcf
    ##fileformat=VCFv4.2
    ##contig=<ID=chr1,length=1000000>
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">
    ##INFO=<ID=JUNCTIONSEQ,Number=1,Type=String,Description="Non-anchor sequence retained from the canonical ESVEE breakend ALT allele; may include reference-derived assembly sequence and does not necessarily represent a novel insertion">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    EOF

    cat <<'EOF' > ${prefix}.standardised.report.tsv
    section	name	value
    summary	input_records	0
    summary	output_records	0
    summary	records_removed_by_merging	0
    summary	merged_pairs	0
    summary	sgl_converted_to_bnd	0
    summary	symbolic_records_with_junctionseq	0
    summary	unchanged_records	0
    summary	warnings	0
    EOF
    """
}