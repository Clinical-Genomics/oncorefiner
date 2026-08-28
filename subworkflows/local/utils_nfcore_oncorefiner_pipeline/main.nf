//
// Subworkflow with functionality specific to the Clinical-Genomics/oncorefiner pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def citations_list = []
    def vcfanno           = "vcfanno (Pedersen et al. 2016)"
    def bcftools_view     = "bcftools (Danecek et al. 2021)"
    def cadd              = "CADD (Rentzsch et al. 2019)"
    def ensemblvep_vep    = "Ensembl VEP (McLaren et al. 2016)"
    def svdb              = "svdb"
    def genmod            = "Genmod"
    def multiqc           = "MultiQC (Ewels et al. 2016)"
    def rmarkdownnotebook = "R Markdown Notebook"

    if (params.snv_vcf) {
        citations_list =
            citations_list +
            vcfanno        +
            bcftools_view  +
            ensemblvep_vep
            if (params.cadd_resources) {
                citations_list = citations_list + cadd
            }

        if (params.genmod_score_config_snv) {
            citations_list =
                citations_list +
                genmod
        }

        if (params.genmod_score_config_sv) {
            citations_list =
                citations_list +
                genmod
        }

    }

    if (params.sv_vcf) {
        citations_list =
            citations_list +
            svdb           +
            bcftools_view  +
            ensemblvep_vep
    }

    if (params.cnv_gene_tsv || params.cnv_segment_tsv) {
        citations_list =
            citations_list +
            rmarkdownnotebook
    }

    // always run
    citations_list =
        citations_list +
        multiqc

    def citations_text =
            "Tools used in the workflow included: "   +
            citations_list.unique().join(', ').trim() +
            "."

    return citations_text
}

def toolBibliographyText() {

    def bibliography_list   = []
    def vcfanno             = "<li>Pedersen BS, Layer RM, Quinlan AR. Vcfanno: fast, flexible annotation of genetic variants. Genome Biol. 2016 Jun 1;17(1):118. doi: 10.1186/s13059-016-0973-5. PMID: 27250555; PMCID: PMC4888505.</li>"
    def bcftools_view       = "<li>Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590845; PMCID: PMC7898596.</li>"
    def cadd                = "<li>Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. CADD: predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Res. 2019 Jan 8;47(D1):D886-D894. doi: 10.1093/nar/gky1016. PMID: 30371827; PMCID: PMC6323892.</li>"
    def ensemblvep_vep      = "<li>McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. doi: 10.1186/s13059-016-0974-4. PMID: 27268795; PMCID: PMC4893825.</li>"
    def svdb                = "<li>svdb. https://github.com/J35P312/svdb.</li>"
    def multiqc             = "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    def rmarkdownnotebook   = "<li>R Markdown Notebook. https://github.com/rstudio/rmarkdown.</li>"

    if (params.snv_vcf) {
        bibliography_list =
            bibliography_list +
            vcfanno        +
            bcftools_view  +
            ensemblvep_vep
            if (params.cadd_resources) {
                bibliography_list = bibliography_list + cadd
            }
    }

    if (params.sv_vcf) {
        bibliography_list =
            bibliography_list +
            svdb           +
            bcftools_view  +
            ensemblvep_vep
    }

    if (params.cnv_gene_tsv || params.cnv_segment_tsv) {
        bibliography_list =
            bibliography_list +
            rmarkdownnotebook
    }

    // always run
    bibliography_list =
        bibliography_list +
        multiqc

    def bibliography_text = bibliography_list.unique().join(' ').trim()

    return bibliography_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

/**
 * Creates a metadata map with provided values, excluding any null values.
 * @param case_id The case ID (required)
 * @param sample_id The sample ID (can be null)
 * @param sample_type The sample type (can be null)
 * @param sex The sex (can be null)
 * @return Metadata map with provided values
 */
def makeMetadata(id, case_id, sample_id=null, sample_type=null, sex=null) {
    [
        id         : id,
        case_id    : case_id,
        sample_id  : sample_id,
        sample_type: sample_type,
        sex        : sex
    ].findAll { key, value -> value != null }
}

/**
 * Creates a value channel from a metadata tuple and file path if provided and the file exists, otherwise returns an empty channel.
 * @param filePath The path to the file (can be null)
 * @param meta The channel metadata tuple
 * @return Value channel with metadata tuple and collected file path or empty channel
 */
def channelFromMetaAndPath(meta, filePath) {
    if (filePath && meta) {
        return channel.fromPath(filePath, checkIfExists: true).map { file -> [meta, file] }.collect()
    }
    if (filePath && !meta) {
        error "Metadata must be provided when a file path is given. Please provide metadata for the file: ${filePath}"
    }
    return channel.empty()
}
