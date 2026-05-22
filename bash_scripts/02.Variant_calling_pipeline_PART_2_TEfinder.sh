#!/bin/bash

# Creation date: 01.12.25
# Based on Lopez Diaz et al (2025), Nat. comm.
# Author: Ana Rodríguez López
# Reference genome used: https://doi.org/10.1128/mbio.00951-25

# With the BAM files obtained after running Picard MarkDuplicates, run TEfinder (Sohrab et al, 2021) to analyze potential new TE insertions in the evolved populations.

# Activate conda environment
source /home/ana/anaconda3/etc/profile.d/conda.sh
conda activate picard_env

# Create file variables
sample_list="" # Text file with sample names
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name="$(echo $line | awk '{print $1}')"
sample_name="$(echo $file_name | awk -F "_" '{print $1}')"

# TE finder pipeline variables
TEfinder_dir="" 
TEfinder_script_path="${TEfinder_dir}/TEfinder" # TEfinder pipeline path
BAM_file="/path/to/directory/${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
reference=""
gtfTE="${TEfinder_dir}/TEs_Fuxo3_assembly.gtf"
teList="${TEfinder_dir}/TE_list_foxy3.txt"
picardPath="${TEfinder_dir}/picard.jar"

# Code
$TEfinder_script_path -alignment $BAM_file -fa $reference -gtf $gtfTE -te $teList -threads 4 -picard $picardPath -maxHeapMem 25000 -workingdir TEfinder_Foxy3_${sample_name} -outname $sample_name 

# Deactivate conda
conda deactivate

# When the script is complete, run bedtools intersect for each population and the WT (-a $Evolved_bed -b $WT_bed -v), to remove the TEs present in the ancestor --> Obtaining the real Transposon Insertion Variations (TIVs)
