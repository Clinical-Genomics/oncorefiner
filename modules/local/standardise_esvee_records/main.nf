process STANDARDISE_ESVEE_RECORDS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91b47b3490e3993dbca42c06d76edbb65b887ac5f9a72bcb31f63471fbdcdecf/data' :
        'community.wave.seqera.io/library/python_click:7a177b12e71d4c56' }"

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
    ##INFO=<ID=INSSEQ,Number=1,Type=String,Description="Inserted sequence">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    EOF

    cat <<'EOF' > ${prefix}.standardised.report.tsv
    section	name	value
    summary	input_records	0
    summary	output_records	0
    summary	records_removed_by_merging	0
    summary	merged_pairs	0
    summary	sgl_converted_to_bnd	0
    summary	symbolic_records_with_insseq	0
    summary	unchanged_records	0
    summary	warnings	0
    EOF
    """
}