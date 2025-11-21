## Table of Contents
* [Introduction](#intro)
* [Quick usage](#quick_usage)
* [Reference sequences](#reference)
  * [Creating a phylogenetic tree](#creating)
  * [`updateNewick.py`](#updatenewick)
  * [Reference indexing](#index)
* [Sample processing](#sample)
  * [Adapter removal](#adapter)
  * [Alignment with `bowtie2`](#align)
  * [`parseSAM.py`](#parsesam)
  * [`MGtree.py`](#mgtree)
* [Usage with Nextflow](#usage)
  * [Input specification](#input)
  * [Nextflow Samplesheet Input](#samplesheet)
  * [Running the pipeline](#running)
  * [Output files Nextflow](#output)
  * [Other Nextflow information and parameters](#nextflow)
*  [Accessory scripts](#accessory)
  * [`treeMaker.py`](#treemaker)
  * [`utils.py`](#utils)
  * [`test.sh`](#test)
* [Miscellaneous](#misc)
<br><br>

<a name="intro"></a>
## Introduction

This pipeline performs metagenomics classification of short reads against a set of reference sequences via a phylogenetic tree.  It is alignment-based, making it particularly suitable for querying against references that are very similar to each other.

The pipeline produces a list of the reference genotypes and relative abundances that are found in a given sample.  Note that the pipeline, and this manual, always refer to the output classifications as `genotypes`, but the output could be at any of a number of levels (e.g. serotype, species, strain, etc.), depending on the reference sequences that are being queried against.

This pipeline can be run from the main branch, which contaims the core scripts and can be used to process a single sample with the alinger of your choice, and a nextflow branch, which allows for batch processing of samples from QC and alignment with `bowtie2`, through classification. 
* For the main branch, use `git checkout MGtree_main`
* For the nextflow branch, use `git checkout MGtree_nextflow`

![MGtree-](https://github.com/user-attachments/assets/d22b715d-cd51-456e-b732-1aded4d9e1c0)
**MGtree workflow.**

The pipeline utilizes several open-source programs:
* [`MEGA-X`](https://www.megasoftware.net/): creates a phylogenetic tree
* [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) v2.3.5.1: aligns the short reads to a set of reference sequences
* [`NGmerge`](https://github.com/jsh58/NGmerge) v0.3: adapter trimming of fastq sequencing files
* [`fastQC`](https://github.com/s-andrews/FastQC) v0.12.1: quality report of sequencing fastq files
* [`fastp`](https://github.com/OpenGene/fastp) v0.23.4: optional QC and trimming utility
* [`MultiQC`](https://multiqc.info) v1.25.1: quality report aggregation

The rest of the pipeline, `MGtree.py` and the other scripts, are written in python with minimal dependencies.  They should work with either python2 or python3, but [everyone should be using python3 now](https://pythonclock.org/).  The scripts are described below, the their usage with Nextflow for batch processing is described. To process samples in batch with Nextflow, the user will also need Nextflow (v24.04.2+) and java (v8+). 
<br><br>


<a name="quick_usage"></a>
## Quick usage

The typical command for running the pipeline is as follows:

With an input samplesheet (additional details below):
The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
```

```bash
nextflow run MGtree --input ./samplesheet.csv --outdir ./results --fasta reference.fasta --newick newick_labeled.nwk -profile docker
```

With a fastq directory:
With an input samplesheet (format described below):
```bash
nextflow run MGtree --fastqdir ./fastq_directory --outdir ./results --fasta reference.fasta --newick newick_labeled.nwk -profile docker
```

The reference and labeled newick file requirements are described below. 
This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

The pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

<a name="reference"></a>
## Reference sequences

Reference sequences should be collected in a single fasta-formatted file.  They will be used to create a phylogenetic tree, and a reference index against which to align the reads.  These procedures need to be done only once for a given set of reference sequences, even if multiple samples are to be processed against them.
<br>

<a name="creating"></a>
#### Creating a phylogenetic tree

Before running `MGtree.py`, a phylogenetic tree must be created from the reference sequences to be queried. The tree must be exported as a [Newick string](https://en.wikipedia.org/wiki/Newick_format) for use by the rest of the pipeline.
<br>

<a name="updatenewick"></a>
#### `updateNewick.py`

This script reads in a [Newick string](https://en.wikipedia.org/wiki/Newick_format) from an input file and outputs the same.  It gives every unnamed node a unique name (an integer).  It also can (re)assign genotypes to nodes that are listed in an optional input csv file.  The genotype assignment is done based on node names, so this script may need to be run twice for a given tree: once (without `-n <file>`) to name every node, and then (with `-n <file>`) to assign genotypes.

```
usage: updateNewick.py -i <file> -o <file> [-n <file>] [-v] [-h]

Required arguments:
  -i <file>      Input Newick file
  -o <file>      Output Newick file

Optional arguments:
  -n <file>      Input csv file to genotype nodes
  -v, --verbose  Run in verbose mode
  -h, --help     Show help message and exit
```
<br>

```
  -i <file>      Input Newick file
```
* The nodes in the input Newick string may contain names.  The names must be unique.  Names are not allowed to contain the following characters: `;(),:`
* Nodes with names that contain `_` will have genotypes interpreted as the first token before `_`.  Two or more `_` characters are not allowed in node names.
<br><br>

```
  -n <file>      Input csv file to genotype nodes
```
* This file should list, on each line, a node name and the genotype to be assigned to it, comma-separated.  Genotypes are not allowed to contain the following characters: `;(),:_`
<br>

<a name="index"></a>
#### Reference indexing

In order to align reads to a set of references, the reference sequences must be indexed.  This pipeline uses the short read aligner `bowtie2`, whose indexes are created by [`bowtie2-build`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).  
<br>


<a name="sample"></a>
## Sample processing

The rest of the pipeline analyzes samples in batch (starting from a directory of fastq files or an input samplesheet) and produces a list of the reference genotypes and relative abundances found therein for each sample.  It requires the phylogenetic tree and reference index previously created.  The pipeline is designed to analyze samples that have been sequenced via some type of massively-parallel sequencing, such as Illumina.  The sequencing can be targeted or shotgun, and single-end or paired-end. The scripts involved are described indiviudually and the usage with nextflow is described below. 

<a name="adapter"></a>
#### Adapter removal

For sequence reads derived from short DNA fragments, the 3' ends may contain portions of a sequencing adapter. This adapter contamination may prevent the reads from aligning to the reference sequences and adversely affect the downstream analysis.  Please consult [this remarkably well-written reference](https://github.com/jsh58/ATAC-seq2#adapter-removal) that discusses two programs for removing adapters: [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) and [NGmerge](https://github.com/jsh58/NGmerge). The user can perform their own adatper removal prior to processing (`--skip_ngmerge`), or `NGmerge` can be implemented in the pipeline. 

<a name="align"></a>
#### Alignment with `bowtie2`

As stated above, this pipeline is particularly suited for querying against references that are very similar.  To this end, it works best when *all* valid alignments of a read/fragment are available to be analyzed.  This can be accomplished with our preferred aligner, [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) in [`-a` mode](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#a-mode-search-for-and-report-all-alignments), which is currently implemented via nextflow in our pipeline. 

<a name="parsesam"></a>
#### `parseSAM.py`

This script analyzes a SAM file, considering primary alignments of reads/fragments, as well as secondary alignments that are equivalent (or within a specified threshold).  For paired-end sequencing, if the primary alignment is in a properly paired configuration, only secondary alignments that are also properly paired are evaluated.  The script produces a file that lists, for each read/fragment, the alignment type and the reference sequence(s) to which it aligns.

```
usage: parseSAM.py -i <file> -o <file> [-m <int>] [-s <float>] [-S] [-v] [-h]

Required arguments:
  -i <file>      Input SAM file (name-sorted, with header; use '-' for stdin)
  -o <file>      Output tsv file of reads, aln types, and refs

Optional arguments:
  -m <int>       Minimum MAPQ (def. 0)
  -s <float>     Keep sec alns with AS >= bestAS - <float> (def. 0)
  -S             Skip name-sorting requirement
  -v, --verbose  Run in verbose mode
  -h, --help     Show help message and exit
```
<br>

```
  -i <file>      Input SAM file (name-sorted, with header; use '-' for stdin)
```
* As stated, the input SAM file must have a header and be sorted by queryname (so the script can match up primary and secondary alignments for reads/fragments).
* The requirement for strict queryname sorting can be waived with `-S`.
<br><br>

```
  -o <file>      Output tsv file of reads, aln types, and refs
```
* The output is a headerless tsv file that lists, for each read/fragment, the name (QNAME), the alignment type, and the reference sequence(s) to which the read/fragment had a valid alignment
* The alignment type can be any of the following:
  * `PE`: for paired-end sequencing, where both reads aligned in a properly paired configuration; only one record is given for the pair of reads
  * `R1`: for the R1 read of a pair, where the reads did not align in a properly paired configuration
  * `R2`: for the R2 read of a pair, where the reads did not align in a properly paired configuration
  * `SE`: for single-end sequencing
* The third column lists the reference sequences, comma-separated.  Two alternative entries are:
  * `unmapped`: the read/fragment did not have a valid alignment to any reference sequence
  * `lowMAPQ`: the read/fragment had a MAPQ value that was too low (see `-m <int>` below)
* Here is a small section of the output file for three read pairs:
```
NDX551337_RUO:11:HLG2WBGXV:1:11102:16719:4459    PE   GII.17_PP423019.1
NDX551337_RUO:11:HLG2WBGXV:1:11102:14718:14360   R1   unmapped
NDX551337_RUO:11:HLG2WBGXV:1:11102:14718:14360   R2   GII.14_OR536456.1,GII.14_OR536453.1,GII.14_OR536457.1
NDX551337_RUO:11:HLG2WBGXV:1:11103:10045:17536   PE   GII.14_OR536456.1,GII.14_OR536453.1,GII.14_OR536457.1
```
<br>

```
  -m <int>       Minimum MAPQ (def. 0)
```
* Alignments with a MAPQ below the `-m` threshold are reported as `lowMAPQ` in the output.  The default value of 0 means that no alignments will be reported as such.
* A good explanation of mapping quality scores can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mapping-quality-higher-more-unique).  Note that the MAPQs reported by `bowtie2` are not meaningful when it is run with `-k` or `-a`.
<br><br>

```
  -s <float>     Keep sec alns with AS >= bestAS - <float> (def. 0)
```
* `parseSAM.py` considers all secondary alignments of multimapping reads/fragments, but, by default, it reports only the alignments whose scores are equal to the best score for the read/fragment.  Setting a value such as `-s 20` causes `parseSAM.py` also to report secondary alignments whose scores are within 20 of the best.  This is equivalent to the [`-s` parameter of Genrich](https://github.com/jsh58/Genrich#sparam).
* The SAM records should have alignment scores under the extra field `AS`.  If not, all alignments are considered equivalent.
* For paired-end sequencing, the alignment score for a fragment's alignment is equal to the sum of the two alignments' individual scores.  Properly paired alignments take precedence over unpaired alignments, regardless of the alignment scores.
<br><br>

<a name="mgtree"></a>
#### `MGtree.py`

This script analyzes a tsv file of alignments and a phylogenetic tree (via a Newick string).  For each read/fragment, it adds counts to the tree at the lowest common ancestor (LCA) node of the set of alignments.  Properly paired alignments are given a weight of two, and unpaired alignments are each counted as one, as are alignments for single-end reads.

Once the alignments are added to the tree, the script interprets the LCA counts for nodes with assigned genotypes.  There are optional output files for counts at leaf nodes (corresponding to the original reference sequences used to build the tree), and for the LCA nodes assigned to each read.

```
usage: MGtree.py -i <file> -n <file> -o <file> [-t <file>] [-q <file>] [-g] [-v] [-h]

Required arguments:
  -i <file>      Input tsv file listing read names, aln types, and references
  -n <file>      Input Newick file
  -o <file>      Output tsv file of genotypes

Optional arguments:
  -t <file>      Output tsv file of leaf node counts
  -q <file>      Output JSON file of reads and assigned LCA nodes
  -g             Do not interpret leaf node names with '_' as having genotypes
  -v, --verbose  Run in verbose mode
  -h, --help     Show help message and exit
```
<br>

```
  -i <file>      Input tsv file listing read names, aln types, and references
```
* This file must be of the format produced by [`parseSAM.py`](#parsesam).
<br><br>

```
  -n <file>      Input Newick file
```
* The nodes in the input Newick string may contain names.  The names must be unique.  Names are not allowed to contain the following characters: `;(),:`
* Nodes with names that contain `_` will have genotypes interpreted as the first token before `_`.  Two or more `_` characters are not allowed in node names.  If `-g` is set, then leaf nodes will *not* have genotypes interpreted.
<br><br>

```
  -o <file>      Output tsv file of genotypes
```
* The primary output is a sorted list of the genotypes observed in the sample, along with the read counts and percentages.  For example:
```
Genotype   Count   Percentage
GII.14     794     70.39%
GII.17     222     19.68%
GII.6      109     9.66%
GI.1       3       0.27%
```
* Only nodes with genotypes are considered.  The reported count is the sum of the LCA counts at the genotype node and every node in its subtree.
* The percentages may add up to more than 100%, if a genotype node has an ancestor node with a genotype.  A warning is printed to `stderr` in such cases.
* LCA counts at a node that is an ancestor of multiple genotype nodes will be reported as ambiguous, such as `ambig[GII.17,GII.14]`.  If no descendant genotypes can be found (including leaf nodes, when run with `-g`), the result will be `ambig[?]`.
* Genotypes may be listed multiple times in the output file, if multiple nodes have the same genotype designation.  If this occurs with leaf nodes, consider using `-g`.
<br>

```
  -t <file>      Output tsv file of leaf node counts
```
* This headerless file has a sorted list of every leaf node (reference sequence) and the count of reads aligning to it.  Reads/fragments with multiple alignments have their counts split evenly among the alignments, so the totals may be fractional.  For example:
```
GII.14_OR536456.1    265.83
GII.14_OR536453.1    265.83
GII.14_OR536457.1    262.33
GII.17_PP423019.1    140.00
GII.17_PP423020.1    82.00
GII.6_OR536442.1     74.50
GII.6_OR536441.1     34.50
GI.1_MK956174.1      2.33
GI.1_MK956175.1      0.33
GI.1_MK956173.1      0.33
```
<br>

```
  -q <file>      Output JSON file of reads and assigned LCA nodes
```
* This file lists the nodes and the reads assigned to them based on the LCA of their alignments, in JSON format.
<br><br>

<a name="usage"></a>
## Usage with Nextflow

<a name="input"></a>
## Input specification

There are multiple ways to specify the input samples data:

- Using `--input`
  - This is the nf-core standard way to specify input files.
  - See section ** Nextflow Samplesheet Input**.
- Using `--input` and `--fastqdir`
  - This is similar to just `--input` alone but the fastq filenames inside the samplesheet are specified relative to `--fastqdir`.
- Using `--fastqdir` without `--input`
  - This will discover all the fastq files in that directory and group them by filename prefix.
  - For best results, filenames should be as generated by `bcl2fastq`, i.e. the expected filename pattern ([regular expression](https://en.wikipedia.org/wiki/Regular_expression)) is `\S+_S\d+_L\d{3}_R[12]_001.fastq.gz`.
- Using `--fastqdir` and `--samplesheet`
  - This is similar to the `--input` and `--fastqdir` specification but the `--samplesheet` argument should be an Illumina experiment `SampleSheet.csv`.
  - The directory should contain fastq generated by `bcl2fastq` without renaming.

<a name="samplesheet"></a>
## Nextflow Samplesheet Input

You will need to create a samplesheet (for nextflow input versus an Illumina SampleSheet) with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column                | Description                                                                                                                                                                            |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`             | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`             | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `is_negative_control` | Boolean, optional. If `true`, sample is treated as a negative control, and genotype calls will be flagged as failures. If missing, samples starting with "NTC" are negative controls.  |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

<a name="running"></a>
## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run MGtree --input ./samplesheet.csv --outdir ./results --fasta reference.fasta --newick newick_labeled.nwk -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run MGtree -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
fasta: CDC_NoV_genotypes.fasta
newick: CDC_NoV_genotypes.nwk
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Pipeline parameters

```
Typical pipeline command:

  nextflow run MGtree -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

--help                                                [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.
--help_full                                           [boolean]         Show the help message for all non-hidden parameters.
--show_hidden                                         [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--help_full`.

Input/output options
  --input                                             [string] Path to comma-separated file containing information about the samples in the experiment.
  --samplesheet                                       [string] Path to Illumina NGS samplesheet to be used with --fastqdir input.
  --fastqdir                                          [string] Path to discover fastq files.
  --outdir                                            [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
  --email                                             [string] Email address for completion summary.
  --multiqc_title                                     [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Reference files options
  --fasta                                             [string] Fasta file containing reference sequences
  --bowtie2                                           [string] Path to prebuilt bowtie2 indices
  --newick                                            [string] Path to file containing a single Newick tree

Pipeline parameters
  --run_fastp                                         [boolean] Run fastp as a first trimming and QC pass
  --skip_fastqc                                       [boolean] Generate FastQC reports
  --skip_ngmerge                                      [boolean] Don't use NGmerge to trim adapters (using dovetail alignments)
  --no_leaf_genotypes                                 [boolean] Ignore genotype information in leaf nodes of the Newick-encoded tree

Generic options
  --multiqc_methods_description                       [string]  Custom MultiQC yaml file containing HTML including a methods description.
  --save_temps                                        [boolean] Publish intermediate files
  --save_index                                        [boolean] Save bowtie2 index when BOWTIE2_BUILD is run

```

<a name="output"></a>
## Output files Nextflow

This section describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

### Pipeline output overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MGtree

<details markdown=1>
<summary>Output files</summary>

- `mgtree/`
  - `*_genotype.tsv`: A TSV file reporting the genotype composition of the sample.
  - `*_tally.tsv`: A headerless TSV file reporting the number of reads mapping to each reference sequence.

#### *_genotype.tsv

A TSV file reporting the genotype composition of the sample. It has three columns:
  - `Genotype` - the genotype call
  - `Count` - number of reads supporting this call
  - `Percentage` - percent of genotyped reads supporting this call

#### *_tally.tsv

A headerless TSV file reporting the number of reads mapping to each reference sequence. It has two columns:
  - `Reference` - the reference accession
  - `Count` - number of reads mapping to this reference

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

Note: Any samples that fail during the pipeline will be reported in `multiqc_report.html`, as well as the reason for the failure.
[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline output information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<a name="nextflow"></a>
## Other Nextflow information and parameters

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull MSDLLCPapers/MGtree
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [MSDLLCPapers/MGtree releases page](https://github.com/MSDLLCPapers/MGtree/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some  setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<a name="accessory"></a>
## Accessory script information

The remaining scripts in the pipeline are for support and do not need to be invoked directly in the normal course of operating the pipeline.

<a name="treemaker"></a>
#### `treeMaker.py`

This script contains constructors and methods for Node and Tree classes.  A Tree is instantiated from a Newick string, and it instantiates all of the Nodes and links.  A Node contains a unique name, a genotype (optionally), pointers to parent and children Nodes, and counters of reads and LCAs.

Invoking this script directly will run a series of unit tests on a toy Tree.

<a name="utils"></a>
#### `utils.py`

This script contains utility functions for file I/O.  Specifying `-` allows the use of `stdin` or `stdout`, and gzip-compressed files are identified by `.gz` suffixes.

<a name="test"></a>
#### `test.sh`

This script systematically tests the pipeline scripts using test files located in the `test-data/` folder:
* `updateNewick.py` without `-n <file>`
* `updateNewick.py` with `-n <file>`
* `parseSAM.py`
* `MGtree.py`

For each command, the output files are checked against expected results.  For MGtree.py, all three output files are checked.  The script prints the diffs of the output files, and a summary of the test results.

<a name="misc"></a>
## Miscellaneous

* All scripts can be run in verbose mode, printing various information and stats to `stderr`.

* Warning messages are ignored at the user's peril.
