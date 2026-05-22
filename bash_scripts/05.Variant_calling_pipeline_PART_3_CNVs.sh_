#!/bin/bash

# Creation date: 12.12.25
# Based on Lopez Diaz et al (2025), Nat. comm.
# Author: Ana Rodríguez López
# Reference genome used: https://doi.org/10.1128/mbio.00951-25

# Use bedtools genomecov (with -d option) for calculate read depth at a singlenucleotide resolution for final populations and WT. 
# The resulting file can be used to determine Copy Number Variations (CNVs)

# Activate conda environment
source /home/ana/anaconda3/etc/profile.d/conda.sh
conda activate bedtools_env

# Create directory variables
input_dir="" # Directory containing files
output_dir=""
logfile_dir=""

#Create file variables
sample_list="" # List with sample names
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name="$(echo $line | awk '{print $1}')"

# Variables for code
input_BAMfile="${input_dir}/${file_name}.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
output_coverage_file="${output_dir}/${file_name}_genomecov.bedgraph"
log_fileBEDtools="${logfile_dir}/bedtoolsGenomeCov_${file_name}.log"

# Software code
{
    echo "BEDtools genomecov running for ${file_name}..."
    bedtools genomecov -d -ibam "$input_BAMfile" > "$output_coverage_file"
    if [[ $? -eq 0 ]]; then
        echo "BEDtools genomecov completed successfully for ${file_name}"
    else
        echo "BEDtools genomecov failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_fileBEDtools" 2>&1

# Deactivate conda
conda deactivate
