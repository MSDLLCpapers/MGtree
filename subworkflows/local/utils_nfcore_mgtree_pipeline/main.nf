//
// Subworkflow with functionality specific to the MSDLLCPapers/MGtree pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
include { BCL2FASTQ               } from '../../../modules/nf-core/bcl2fastq/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    bcl_samplesheet
    fastqdir
    params_map        //     map: copy of the params

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
        params_map,
    )


    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )


    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(params_map)

    //
    // Create channel from input file provided through params.input
    //

    if (input) {
        Channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [meta.id, meta + [single_end: true], [fastq_1]]
                }
                else {
                    return [meta.id, meta + [single_end: false], [fastq_1, fastq_2]]
                }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map { meta, fastqs ->
                fastqs = fastqs
                    .flatten()
                    .collect { fq ->
                        def resolved = fastqdir ? file(fastqdir).resolve(fq) : file(fq)
                        if (!resolved.exists()) {
                            error("unable to find fastq file for sample ${meta.id}: ${resolved.toUriString()}")
                        }
                        resolved
                    }

                return [meta, fastqs]
            }
            .set { ch_samplesheet }
    }
    else if (fastqdir && bcl_samplesheet) {
        rows_with_header = file(bcl_samplesheet)
            .splitCsv()
            .dropWhile { row -> row[0] != '[Data]' }
            .drop(1)
            .takeWhile { row -> row[0] != '' && !row[0].startsWith('[') }
        rows = rows_with_header[1..-1].collect { row -> [rows_with_header[0], row].transpose().collectEntries() }
        ch_samplesheet = Channel.fromList(
            rows.collect { row ->
                def fastq = file("${fastqdir}/${row.Sample_Name}_S[1-9]*_L[0-9]*_R[12]_001.fastq.gz")
                if (!fastq) {
                    error("Fastq for sample ${row.Sample_Name} not found in ${fastqdir}, or this is not bcl2fastq output")
                }
                def meta = [id: row.Sample_Name, single_end: fastq.find { fq -> fq.name =~ /_R2_001\.fastq\.gz$/ } == null]
                return [meta, fastq.sort { it -> it.baseName }]
            }
        )
    }
    else if (fastqdir) {
        ch_samplesheet = Channel
            .fromPath("${fastqdir}/*[._]f{ast,}q.gz")
            .filter { fq -> !fq.baseName.startsWith('Undetermined_S0') }
            .map { fq ->
                log.debug("${fq.name}")
                def patterns = [
                    /^(\S+)_S\d+_L\d+_R[12](_001)?\.fastq\.gz$/,
                    /^(\S+)[._][12][._]f(ast)?q\.gz$/,
                ]
                def id = patterns.findResult { pt -> fq.name =~ pt }?.getAt(0)?.getAt(1) ?: fq.name
                return [id, fq]
            }
            .groupTuple()
            .map { id, fastq ->
                def meta = [id: id, single_end: fastq.find { fq -> fq.name =~ /(_R2_001\.|[_.]2[_.])f(ast)?q\.gz$/ } == null]
                return [meta, fastq.flatten().sort { it -> it.baseName }]
            }
    }
    else {
        ch_samplesheet = Channel.empty()
    }

    ch_samplesheet = ch_samplesheet.map { meta, fastq ->
        [meta + [negative_control: meta.negative_control != null ? meta.negative_control : meta.id.startsWith('NTC')], fastq]
    }

    ch_samplesheet.ifEmpty {
        error("No samples found (probable cause: improper specification of --input, --samplesheet, and/or --fastqdir)")
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email                  //  string: email address
    email_on_fail          //  string: email address sent on pipeline failure
    plaintext_email        // boolean: Send plain-text email instead of HTML
    max_multiqc_email_size // MemoryValue
    outdir                 //    path: Path to output directory where results will be published
    monochrome_logs        // boolean: Disable ANSI colour codes in log output
    hook_url               //  string: hook URL for notifications
    multiqc_report         //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

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
                multiqc_report.toList(),
                max_multiqc_email_size,
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters(params) {
    genomeExistsError(params)
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(params, attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Get attribute from genome config file e.g. fasta
//
def getControlGenomeAttribute(params, attribute) {
    if (params.genomes && params.control_genome && params.genomes.containsKey(params.control_genome)) {
        if (params.genomes[params.control_genome].containsKey(attribute)) {
            return params.genomes[params.control_genome][attribute]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError(params) {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "FastP (Chen 2023)",
        "NGmerge (Gaspar 2018)",
        "Bowtie2 (Langmead 2012)",
        "Samtools (Li 2009)",
        "Bedtools (Quinlan 2010)",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Chen S, (2023) Ultrafast one-pass FASTQ data preprocessing, quality contorl, and deduplication using fastqp. iMeta 2: e107. https://doi.org/10.1002/imt2.107</li>",
        "<li>Gaspar JM. NGmerge: merging paired-end reads via novel empirically-derived models of sequencing errors. BMC Bioinformatics. 2018 Dec 20;19(1):536.</li>",
        "<li>Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.</li>",
        "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9</li>",
        "<li>Quinlan AR, Hall IM. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6) 841–842. doi: 10.1093/bioinformatics/btq033</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        def temp_doi_ref = manifest_doi
            .collect { doi_ref ->
                "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
            }
            .join('')
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
