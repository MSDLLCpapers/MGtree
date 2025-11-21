/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { GENOTYPESAMPLE         } from '../subworkflows/local/genotypesample/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mgtree_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGTREE {
    take:
    ch_samplesheet                  // channel: samplesheet read in from --input
    val_fasta
    val_bowtie2
    val_newick
    val_run_fastp
    val_skip_fastqc
    val_skip_ngmerge
    val_no_leaf_genotypes
    val_outdir
    val_multiqc_config
    val_multiqc_logo
    val_multiqc_methods_description

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Module: Run bowtie2-build
    //
    ch_fasta = Channel.value([[id: 'ref'], file(val_fasta)])
    ch_bt2_build_in = Channel.empty()
    if (val_bowtie2) {
        ch_bowtie2 = Channel.value([[id: 'ref'], file(val_bowtie2)])
    }
    else {
        ch_bowtie2 = Channel.empty()
        ch_bt2_build_in = ch_bt2_build_in.mix(ch_fasta)
    }

    ch_control_fasta = Channel.empty()
    ch_control_bowtie2 = Channel.empty()
    BOWTIE2_BUILD(
        ch_bt2_build_in
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    BOWTIE2_BUILD.out.index
        .branch { meta, _index ->
            ref: meta.id == 'ref'
            unknown: true
        }
        .set { ch_bt2_build_out }
    ch_bowtie2 = ch_bowtie2.mix(ch_bt2_build_out.ref).collect()
    ch_newick = Channel.value(file(val_newick))

    //
    // SUBWORKFLOW: Genotype sample
    //
    GENOTYPESAMPLE(
        ch_samplesheet,
        ch_fasta,
        ch_bowtie2,
        ch_newick,
        val_run_fastp,
        val_skip_fastqc,
        val_skip_ngmerge,
        val_no_leaf_genotypes,
        val_outdir,
    )
    ch_versions = ch_versions.mix(GENOTYPESAMPLE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        GENOTYPESAMPLE.out.fastqc_raw_zip.collect { it -> it[1] },
        GENOTYPESAMPLE.out.fastp_trim_log.collect { it -> it[1] },
        GENOTYPESAMPLE.out.fastqc_trim_zip.collect { it -> it[1] },
        GENOTYPESAMPLE.out.bowtie2_log.collect { it -> it[1] },
        GENOTYPESAMPLE.out.mgtree_count_json_mqc.collect(),
        GENOTYPESAMPLE.out.genotype_failed_meta_mqc,
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name: '' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = val_multiqc_config
        ? Channel.fromPath(val_multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = val_multiqc_logo
        ? Channel.fromPath(val_multiqc_logo, checkIfExists: true)
        : Channel.empty()


    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    ch_multiqc_custom_methods_description = val_multiqc_methods_description
        ? file(val_multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
