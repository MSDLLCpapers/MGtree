#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MSDLLCPapers/MGtree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/MSDLLCPapers/MGtree
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MGTREE                    } from './workflows/mgtree'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_mgtree_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_mgtree_pipeline'
include {
    getGenomeAttribute ;
    getControlGenomeAttribute
} from './subworkflows/local/utils_nfcore_mgtree_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.samplesheet,
        params.fastqdir,
        params,
    )

    //
    // WORKFLOW: Run main workflow
    //
    MGTREE(
        PIPELINE_INITIALISATION.out.samplesheet,
        params.fasta ?: getGenomeAttribute(params, 'fasta'),
        params.bowtie2 ?: getGenomeAttribute(params, 'bowtie2'),
        params.newick ?: getGenomeAttribute(params, 'newick'),
        params.run_fastp,
        params.skip_fastqc,
        params.skip_ngmerge,
        params.no_leaf_genotypes,
        params.outdir,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.getOrDefault('max_multiqc_email_size', 0) as MemoryUnit,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MGTREE.out.multiqc_report,
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow MGTREE {
    take:
    samplesheet                     // channel: samplesheet read in from --input
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

    //
    // WORKFLOW: Run pipeline
    //
    MGTREE(
        samplesheet,
        val_fasta,
        val_bowtie2,
        val_newick,
        val_run_fastp,
        val_skip_fastqc,
        val_skip_ngmerge,
        val_no_leaf_genotypes,
        val_outdir,
        val_multiqc_config,
        val_multiqc_logo,
        val_multiqc_methods_description,
    )

    emit:
    multiqc_report = MGTREE.out.multiqc_report // channel: /path/to/multiqc_report.html
}
