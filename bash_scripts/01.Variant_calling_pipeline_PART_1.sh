#!/bin/bash

# Creation date: 26.11.25
# Based on Lopez Diaz et al (2025), Nat. comm.
# Author: Ana Rodríguez López
# Reference genome used: https://doi.org/10.1128/mbio.00951-25

# First part of the pipeline for variant calling. Run this part first and check the statistics at the end before continuing with the next parts

# After a first round of FastQC using the raw data, clean reads following the pipeline: samplesDNA_pretreatments.sh

# Run FastQC again before mapping. If everything is correct, continue with the pipeline.

# Create file variables
sample_list="path" # A list containing sample names. Example can be found in 
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name="$(echo $line | awk '{print $1}')"
sample_name="$(echo $file_name | awk -F "_" '{print $1}')"

# Create directory variables
clean_samples_dir="/data/ana/Sequencing_Novogene/YP_YC_DNAseq_analysis/input"
output_dir="/data/ana/Sequencing_Novogene/YP_YC_DNAseq_analysis/output"
logfile_dir="/data/ana/Sequencing_Novogene/YP_YC_DNAseq_analysis/logfiles"

    ## Mapping to the reference genome using BWAmem

# Activate bwa_env
source /home/ana/anaconda3/etc/profile.d/conda.sh
conda activate bwa_env

# BWA-mem variables
refseq="/data/ana/Sequencing_Novogene/YP_YC_DNAseq_analysis/input/foxy3_bwa_index"
forward_bwa="${clean_samples_dir}/${file_name}_1_cleanP.fq.gz"
reverse_bwa="${clean_samples_dir}/${file_name}_2_cleanP.fq.gz"
sam_file="${output_dir}/${file_name}.bwa.sam"

# Create Read Groups
header=$(zcat $forward_bwa | head -n 1)
ID="${sample_name}_L$(echo "$header" | awk -F':' '{print $4}')"
SM="$sample_name"
PL="ILLUMINA"
read_group="@RG\\tID:$ID\\tSM:$SM\\tPL:$PL"

# Log file BWA 
log_bwa="${logfile_dir}/bwamem_${file_name}.log"

# Verify that the input file exists
if [[ ! -f "$forward_bwa" ]]; then
    echo "Error: Input file ${forward_bwa} does not exist." >&2
    exit 1
fi

if [[ ! -f "$reverse_bwa" ]]; then
    echo "Error: Input file ${reverse_bwa} does not exist." >&2
    exit 1
fi

#Software code
{
    echo "Aligning ${file_name}_1_cleanP.fq.gz and ${file_name}_2_cleanP.fq.gz to the reference genome (Ayhan et al, 2025)"
    bwa mem -t 8 -M -R "$read_group" -a -v 1 "$refseq" "$forward_bwa" "$reverse_bwa" > "$sam_file" # I have added -R option to add Read Groups information (neccesary for later steps with GATK)
    if [[ $? -eq 0 ]]; then
        echo "BWA mem completed successfully for ${file_name}"
    else
        echo "BWA mem failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_bwa" 2>&1

#Deactivate conda
conda deactivate

    ## Picard CleanSAM (Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads), SAM to BAM file conversion and SAM removal

    ## PICARD CLEANSAM

# Activate picard_env
conda activate picard_env

# Create log file
log_CS_StB=${logfile_dir}/PicardCleanSAM_SAMtoBAM_${file_name}.log

# Verify that the input file exists
if [[ ! -f "$sam_file" ]]; then
    echo "Error: Input file ${sam_file} does not exist." >> "$log_CS_StB" 2>&1
    exit 1
fi

# Software code
cleanSAM_file=${output_dir}/${file_name}.bwa.cleanSAM.sam
{
    echo "Cleaning SAM file ${file_name}.bwa.sam"
    picard CleanSam I=$sam_file O=$cleanSAM_file

    if [[ $? -eq 0 ]]; then
        echo "CleanSam completed successfully for ${file_name}.bwa.sam"
    else
        echo "CleanSam failed for ${file_name}.bwa.sam"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_CS_StB" 2>&1

# Deactivate conda
conda deactivate

# Remove raw SAM files
rm $sam_file
echo "$sam_file has been removed succesfully" >> "$log_CS_StB" 2>&1

    ## SAMTools VIEW AND SORT (SAM to BAM) and SAM removal

# Activate conda environment
conda activate samtools_env

# Verify input files exist
if [[ ! -f "$cleanSAM_file" ]]; then
    echo "Error: Cleaned SAM file ${cleanSAM_file} does not exist." >&2
    exit 1
fi

# Software code
bam_file=${output_dir}/${file_name}.bwa.cleanSAM.bam
{
    echo "Coverting SAM to BAM for ${file_name}.bwa.cleanSAM.sam"
    samtools view -bS $cleanSAM_file | samtools sort -o $bam_file
    if [[ $? -eq 0 ]]; then
        echo "Samtools view and sort completed successfully for ${file_name}.bwa.cleanSAM.sam"
    else
        echo "Samtools view and sort failed for ${file_name}.bwa.cleanSAM.sam"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_CS_StB" 2>&1

# Deactivate conda
conda deactivate

# Remove SAM files
rm $cleanSAM_file
echo "$cleanSAM_file has been removed succesfully" >> "$log_CS_StB" 2>&1

    ## Picard FixMateInformation (Verify mate-pair information between mates and fix if needed) and Picard MarkDuplicates (Locate and tag duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA)

# Activate picard conda environment
conda activate picard_env

    ## PICARD FIXMATEINFORMATION

# Log file  
log_picard_fixMate="${logfile_dir}/picard_FixMateInfo_${file_name}.log"

# Verify that the input file exists
if [[ ! -f "$bam_file" ]]; then
    echo "Error: Input file ${bam_file} does not exist." >> "$log_picard_fixMate" >&2
    exit 1
fi
#Software codes
bam_file_FMI="${output_dir}/${file_name}.bwa.cleanSAM.fixMate.bam"
{
    echo "Processing ${file_name}.bwa.cleanSAM.bam for Fix Mate Information..."
    picard FixMateInformation I="$bam_file" O="$bam_file_FMI"
    if [[ $? -eq 0 ]]; then
        echo "Picard Fix Mate Information completed successfully for ${file_name}"
    else
        echo "Picard Fix Mate Information failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_picard_fixMate" 2>&1 

    ## PICARD MARKDUPLICATES and BAM indexing

# Log file  
log_picard_markDup_index="${logfile_dir}/picard_markDup_index_${file_name}.log"

# Software codes
bam_file_MD="${output_dir}/${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
metrics_file_MD="${output_dir}/${file_name}_markDupMetrics.txt"
{
    echo "Executing Picard MarkDuplicates for ${file_name}.bwa.cleanSAM.fixMate.bam"
    picard MarkDuplicates I="$bam_file_FMI" O="$bam_file_MD" M="$metrics_file_MD"
    if [[ $? -eq 0 ]]; then
        echo "Picard Mark Duplicates completed successfully for ${file_name}"
    else
        echo "Picard Mark Duplicates failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_picard_markDup_index" 2>&1

# Deactivate conda
conda deactivate

# Activate conda environment
conda activate samtools_env

# Software code
{
    echo "Indexing ${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
    samtools index "$bam_file_MD"
    if [[ $? -eq 0 ]]; then
        echo "Indexing completed successfully for ${file_name}"
    else
        echo "Indexing failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_picard_markDup_index" 2>&1

# Deactivate conda
conda deactivate

    ## Collecting mapping stats with SAMtools Flagstats and Picard CollectRawWgsMetrics. If you have samples sequenced in multiple lanes, merged them before running the following codes!!! Use this code:

# sbatch --wrap="samtools merge -o YC2_MKDN240008445-1A_22CHMMLT4_L26.bwa.cleanSAM.fixMate.markDup.bam YC2_MKDN240008445-1A_227Y3HLT4_L2.bwa.cleanSAM.fixMate.markDup.bam YC2_MKDN240008445-1A_22CHMMLT4_L6.bwa.cleanSAM.fixMate.markDup.bam"
# Check that it was done correctly with:
# samtools view -H YC2_MKDN240008445-1A_22CHMMLT4_L26.bwa.cleanSAM.fixMate.markDup.bam | grep '^@RG' # It should appear read groups from the different lanes
# samtools view -c YC2_MKDN240008445-1A_227Y3HLT4_L2.bwa.cleanSAM.fixMate.markDup.bam # Check nº of mapped reads 
# samtools view -c YC2_MKDN240008445-1A_22CHMMLT4_L6.bwa.cleanSAM.fixMate.markDup.bam # Check nº of mapped reads 
# samtools view -c YC2_MKDN240008445-1A_22CHMMLT4_L26.bwa.cleanSAM.fixMate.markDup.bam # Check nº of mapped reads 
# Calculate if the number of mapped lectures in the fusion file correspond with the sum of the separated files


    ## SAMtools Flagstats (collects statistics from BAM files and outputs in a text format)

# Activate conda environment
conda activate samtools_env

# Samtools flagstats code
samtools flagstat -@ 2 "$bam_file_MD" > "${output_dir}/${file_name}_flagstat.txt"

# Deactivate conda
conda deactivate

    ## Picard CollectRawWgsMetrics (Collect whole genome sequencing-related metrics.  This tool computes metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments. These metrics include the percentages of reads that pass minimal base- and mapping- quality filters as well as coverage (read-depth) levels)

# Reference genome variable
reference_genome="/data/ana/Sequencing_Novogene/YP_YC_DNAseq_analysis/input/fusox_3_1_AssemblyScaffolds.fasta"

# Activate conda environment
conda activate picard_env

# Picard code
picard CollectRawWgsMetrics -I "$bam_file_MD" -O "${output_dir}/${file_name}_picardRawWgsMetrics.txt" -R "$reference_genome" --INCLUDE_BQ_HISTOGRAM true

# Deactivate conda
conda deactivate
