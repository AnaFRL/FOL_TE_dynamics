#!/bin/bash

# Creation date: 01.12.25
# Based on Lopez Diaz et al (2025), Nat. comm.
# Author: Ana Rodríguez López
# Reference genome used: Fol4287 2025 version

# Use the BAM files obtained after running Picard MarkDuplicates to analyze Single Nucleotide Variants (SNVs) and small Insertion/Deletions (INDELs) with GATK software

# Create file variables
sample_list="" # Lisy with sample names
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name="$(echo $line | awk '{print $1}')"

# Create directory variables
markDup_BAM="/path/to/directory/${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
output_dir=""
logfile_dir=""

    ## Run FreeBayes for an initial variant calling. --> FreeBayes performs variant calling on .bam files that contain aligned DNA sequence data. FreeBayes analyzes the alignment of sequences against a reference genome to detect variants, such as SNPs (single nucleotide polymorphisms), indels (insertions and deletions), and other genetic variations in the DNA.

# Activate conda environment
source /home/ana/anaconda3/etc/profile.d/conda.sh
conda activate freebayes_env

# Log file FreeBayes
logfile_FreeBayesVC="${logfile_dir}/freebayesVC_${file_name}.log"

# Software code
reference=""
reference_fai="" # You can obtain this file running samtools index $BAMfile
freebayes_file="${output_dir}/${file_name}.bwa.cleanSAM.fixMate.markDup.bam.fbayes.vcf"
{
    echo "Running FreeBayes for ${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
    freebayes-parallel  <(fasta_generate_regions.py "$reference_fai" 100000) 16 -f "$reference" "$markDup_BAM" >"$freebayes_file" # 16 CPUs
    if [[ $? -eq 0 ]]; then
        echo "Free Bayes Variant Calling completed successfully for ${file_name}"
    else
        echo "Free Bayes Variant Calling failed for ${file_name}"
        exit 1  # Exit the script if it fails
    fi
} >> "$logfile_FreeBayesVC" 2>&1

    # NOTES--> fasta_generate_regions.py "$reference_fai" 100000: This script splits the reference file (indexed with .fai) into intervals of 100,000 base pairs. This allows FreeBayes to parallelize the process, dividing the analysis into "regions" to make variant calling faster and more efficient.
    # freebayes-parallel: Runs FreeBayes in parallel across the specified regions, using 16 threads.

# Zip the output file and index it
bgzip -c $freebayes_file > ${freebayes_file}.gz; tabix -p vcf ${freebayes_file}.gz

# Deactivate conda
conda deactivate

    ## GATK BaseRecalibrator --> Generates recalibration table for Base Quality Score Recalibration (BQSR)

# Before running GATK, a reference sequence dictionary is neccesary (save it in the same directory where the reference sequence is saved). Run: picard CreateSequenceDictionary R="/reference/genome/path"  O="/output/directory/path/reference.dict"

# Activate conda environment
conda activate gatk4_env

# Analyze patterns of covariation in the sequence dataset
# Variables
input_knownSites="${output_dir}/${file_name}.bwa.cleanSAM.fixMate.markDup.bam.fbayes.vcf.gz"
recalData_File="${output_dir}/${file_name}_recalData.table"
logfile_GATKBaseRecal="${logfile_dir}/gatk4BaseRecal_${file_name}.log"

# Check if input files exist
if [[ ! -f "$markDup_BAM" || ! -f "$input_knownSites" || ! -f "$reference" ]]; then
    echo "Error: One or more input files are missing. Check $markDup_BAM, $input_knownSites, or $reference."
    exit 1
fi

# Software code
{
    echo "Creating GATKReport covariation data for ${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
    gatk BaseRecalibrator -R "$reference" -I "$markDup_BAM" --known-sites "$input_knownSites" -O "$recalData_File" #--filter_bases_not_stored #include this last option to avoid errors due to "empty" sequences in the BAMfile (in GATK3)
    if [[ $? -eq 0 ]]; then
        echo "GATK4 recalibration data recovery completed successfully for ${file_name}"
    else
        echo "GATK4 recalibration data recovery failed for ${file_name}"
        exit 1  
    fi 
} >> "$logfile_GATKBaseRecal" 2>&1
    # This creates a GATKReport file called “$recalData_File” containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data. 

    # Apply the recalibration to your sequence data.
#Variables
recalibrated_BAM="${output_dir}/${file_name}.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"

# Software code
{
    echo "Recalibrating bases for ${file_name}.bwa.cleanSAM.fixMate.markDup.bam"
    gatk ApplyBQSR -R "$reference" -I "$markDup_BAM" --bqsr-recal-file "$recalData_File" -O "$recalibrated_BAM" #--filter_bases_not_stored #include this last option to avoid errors due to "empty" sequences in the BAMfile (in GATK3)
    if [[ $? -eq 0 ]]; then
        echo "GATK4 recalibration completed successfully for ${file_name}"
    else
        echo "GATK4 recalibration failed for ${file_name}"
        exit 1  
    fi 
} >> "$logfile_GATKBaseRecal" 2>&1

# Deactivate conda
conda deactivate

    #This creates a file called “$recalibrated_BAM” containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down.
