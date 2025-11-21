#!/bin/bash
#$ -pe threads 4
#$ -l mem_reserve=8G,h_vmem=200G,fourperhost=1
#$ -N Run_dir_process
#$ -o $JOB_NAME.$JOB_ID.out -e $JOB_NAME.$JOB_ID.err
#$ -P vlong

# Run genotyping analysis on a run directory containing fastq files
    # Runs genotypeSample.sh and runSummarize.py
# Author: Samantha Sholes
# Version 1; 2024-10-21

# Nextflow wrapper. Requires nextflow 24.04.2 or later
# Author: Scott Norton
# Version 2; 2024-11-12

module purge
module load singularity/3.10.5

# Set defaults
nf_profile=""
fastq_folder=""
output_folder=""
fastp_flag=""
ngmerge_flag=""
g_flag=""
scripts_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Display help information
usage="Usage: $0 {req_args} [opt_args]
  Required arguments:
    -i <fol>   Input fastq folder
    -r <file>  Reference sequence file
    -n <file>  Newick string file
    -o <fol>   Output folder
  Optional arguments:
    -s <fol>   Folder containing MGtree and accessory scripts
    -f         Run fastp preprocessing step
    -N         Skip NGmerge step
    -g         For MGtree, do not interpret leaf node names with '_'
                 as having genotypes
    -p <str>   Comma-separated list of nextflow configuration profiles
                 to use (def. $nf_profile)"

# Parse options
while getopts ":hi:r:n:o:s:fNt:c:d:gO:p:" opt; do
  case $opt in
    h )
      echo "$usage" 1>&2
      exit 0 ;;
    i )
      fastq_folder=$OPTARG ;;
    r )
      reference=$OPTARG ;;
    n )
      newick=$OPTARG ;;
    o )
      output_folder=$OPTARG ;;
    s )
      scripts_path=$OPTARG ;;
    f )
      fastp_flag="--run_fastp" ;;
    N )
      ngmerge_flag="--skip_ngmerge" ;;
    g )
      g_flag="--no_leaf_genotypes" ;;
    p )
        nf_profile=$OPTARG ;;
    : )
      echo -e "Invalid option: -$OPTARG requires an argument\n$usage" 1>&2
      exit 255 ;;
    \? )
      echo -e "Error! Invalid option\n$usage" 1>&2
      exit 255 ;;
  esac
done

# Shift away the processed options
shift $((OPTIND - 1))
if [ -n "$1" ]; then
  echo -e "Invalid args: $*\n$usage" 1>&2
  exit 255
fi

# Check if required arguments are provided
if [[ -z $fastq_folder ]]; then
    echo "Input fastq folder does not exist. Please provide a valid directory path."
    echo "$usage" 1>&2
    exit 255
fi

if [[ -z $output_folder ]]; then
    echo "Output folder not specified. Please provide a valid directory path."
    echo "$usage" 1>&2
    exit 255
fi

if [[ -z $reference ]]; then
    echo "Reference file not specified. Please provide a valid file path."
    echo "$usage" 1>&2
    exit 255
fi

if [[ ! -f $newick ]]; then
    echo "Newick file \"$newick\" does not exist. Please provide a valid file path."
    echo "$usage" 1>&2
    exit 255
fi

profile_arg=''
if [[ -n "$nf_profile" ]]; then
    profile_arg="-profile $nf_profile"
fi

# Loop over R1 reads in the folder
echo 'sample,fastq_1,fastq_2' > "$output_folder/samplesheet.csv"
for x in "$fastq_folder"/*_R1.fastq.gz; do
    # Get the R2 file path based on the current R1 file
    t2="${x/_R1/_R2}"

    # Print the file paths to confirm
    echo "Checking for paired files: $x and $t2" 1>&2

    # Check if the R2 file exists
    if [ -e "$t2" ]; then
        # Extract the base name without extension
        base_name=$(basename "${x}" _R1.fastq.gz)
        echo "$base_name,$x,$t2" >> "$output_folder/samplesheet.csv"

    else
        echo "Error: Paired file not found for $x"
    fi
done

# shellcheck disable=SC2086
nextflow run "$scripts_path/main.nf" \
    $profile_arg \
    --input "$output_folder/samplesheet.csv" \
    --fasta "$reference.fasta" \
    --bowtie2 "$(dirname $reference)" \
    --newick "$newick" \
    --outdir "$output_folder" \
    $fastp_flag $ngmerge_flag $g_flag
