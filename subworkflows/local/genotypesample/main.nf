include { FASTQ_TRIM_FASTP_FASTQC    } from '../../../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { NGMERGE                    } from '../../../modules/nf-core/ngmerge/main'
include { BOWTIE2_ALIGN              } from '../../../modules/nf-core/bowtie2/align/main'
include { PARSESAM                   } from '../../../modules/local/parsesam/main'
include { MGTREE_COUNT               } from '../../../modules/local/mgtree/count/main'

workflow GENOTYPESAMPLE {
    take:
    ch_samplesheet
    ch_fasta
    ch_bowtie2
    ch_newick
    val_run_fastp
    val_skip_fastqc
    val_skip_ngmerge
    val_no_leaf_genotypes
    val_outdir

    main:

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Run FastP and FastQC
    //
    FASTQ_TRIM_FASTP_FASTQC(
        ch_samplesheet,
        [],
        false,
        true,
        false,
        !val_run_fastp,
        val_skip_fastqc,
    )
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)
    ch_reads = FASTQ_TRIM_FASTP_FASTQC.out.reads

    //
    // MODULE: Run NGMerge
    //
    if (!val_skip_ngmerge) {
        // NGmerge requires paired-end
        ch_reads
            .branch { meta, _fastq ->
                single_end: meta.single_end
                paired_end: true
            }
            .set { ch_reads_by_arity }
        NGMERGE(
            ch_reads_by_arity.paired_end
        )
        ch_versions = ch_versions.mix(NGMERGE.out.versions.first())
        ch_reads = NGMERGE.out.unstitched_read1
            .join(NGMERGE.out.unstitched_read2)
            .map { meta, r1, r2 ->
                [meta, [r1, r2]]
            }
            .mix(ch_reads_by_arity.single_end)
    }

    //
    // Module: Run bowtie2
    //
    BOWTIE2_ALIGN(
        ch_reads,
        ch_bowtie2,
        ch_fasta,
        false,
        false,
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Module: SAM PE to TSV
    //
    PARSESAM(
        BOWTIE2_ALIGN.out.sam
    )
    ch_versions = ch_versions.mix(PARSESAM.out.versions.first())
    ch_parse_tsv = PARSESAM.out.parsed
    ch_sam = PARSESAM.out.sam

    //
    // Module: Count reads on tree
    //
    MGTREE_COUNT(
        ch_parse_tsv,
        ch_newick,
        val_no_leaf_genotypes,
    )
    ch_versions = ch_versions.mix(MGTREE_COUNT.out.versions.first())
    MGTREE_COUNT.out.console
        .filter { _meta, message ->
            !message.isEmpty()
        }
        .subscribe { meta, message ->
            log.info("MGTREE_COUNT: ${meta.id}: ${message}")
        }
    MGTREE_COUNT.out.genotype
        .ifEmpty {
            error("No samples mapped to the reference")
        }
        .map { meta, table ->
            def lines = table.readLines()
            def header = lines[0].split('\t')
            return [
                meta,
                lines[1..-1].collect { line ->
                    [
                        header,
                        line.split('\t'),
                    ].transpose().collectEntries()
                },
            ]
        }
        .set { ch_genotype_read_in }
    ch_samplesheet
        .join(ch_genotype_read_in, remainder: true)
        .map { it -> [it[0], it[-1]] }
        .branch { meta, geno ->
            negative_control_called: meta.negative_control && geno != null
            negative_control_ok: meta.negative_control
            positive_not_called: geno == null
            positive_ok: true
        }
        .set { ch_genotype_branch }

    ch_genotype_branch.positive_ok
        .toList()
        .dump(tag: 'genotype_counts')
        .collectFile(name: "mgtree_genotype_mqc.json") { maps ->
            def json_map = [
                id: "mgtree_genotypes",
                section_name: "Genotype calls",
                description: "These bar graphs show the percent composition of each sample",
                plot_type: "bargraph",
                pconfig: [
                    id: "mgtree_genotype_barplots",
                    title: "Genotype predictions by MGtree",
                    ylab: "Genotype composition",
                ],
                data: maps.collectEntries { meta, rows ->
                    [
                        "${meta.id}",
                        rows?.collectEntries { row ->
                            return [row.Genotype, row.Count]
                        } ?: [undetermined: 0],
                    ]
                },
            ]
            return groovy.json.JsonOutput.toJson(json_map)
        }
        .dump(tag: 'mgtree_count_json_mqc')
        .set { mgtree_count_json_mqc }

    Channel
        .empty()
        .mix(
            ch_genotype_branch.negative_control_called.map { meta, _geno -> meta + [reason: "Negative control fail"] },
            ch_genotype_branch.positive_not_called.map { meta, _geno -> meta + [reason: "No genotypes"] },
        )
        .dump(tag: 'genotype_failed_meta')
        .tap { ch_genotype_failed_meta }
        .collectFile(
            name: 'no_genotype_samples_list.txt',
            storeDir: "${val_outdir}/mgtree",
            newLine: true,
            seed: 'sample\treason',
        ) { meta ->
            "${meta.id}\t${meta.reason}"
        }
        .collectFile(
            name: 'no_genotype_samples_mqc.tsv',
            seed: '''\
            # id: 'failed_genotypes'
            # section_name: 'Samples with no genotype called'
            # description: 'This is a list of samples for which genotyping failed, and the reason for the failure. This usually happens when the sample does not map to the reference library, usually when the sample is negative for the pathogen in question. Negative control samples that are nevertheless mapped to a genotype will also be reported here.'
            # plot_type: "table"
            '''.stripIndent(),
        ) { fl -> fl.text }
        .set { ch_genotype_failed_meta_mqc }

    emit:
    fastqc_raw_zip           = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip
    fastqc_trim_zip          = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip
    fastp_trim_log           = FASTQ_TRIM_FASTP_FASTQC.out.trim_log
    trimmed_reads            = ch_reads
    sam                      = ch_sam
    genotype                 = MGTREE_COUNT.out.genotype
    tally                    = MGTREE_COUNT.out.tally
    bowtie2_log              = BOWTIE2_ALIGN.out.log
    genotype_failed_meta     = ch_genotype_failed_meta
    genotype_failed_meta_mqc = ch_genotype_failed_meta_mqc
    mgtree_count_json_mqc
    versions                 = ch_versions // channel: [ versions.yml ]
}
