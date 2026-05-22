#!/bin/bash

# Creation date: 01.12.25
# Based on Lopez Diaz et al (2025), Nat. comm.
# Author: Ana Rodríguez López
# Reference genome used: https://doi.org/10.1128/mbio.00951-25

    ## GATK Mutect2: Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations. 
    
# Activate conda
source /home/ana/anaconda3/etc/profile.d/conda.sh
conda activate gatk4_env

# Create directory variables
files_dir="" # DIrectory containing files
output_dir=""
logfile_dir=""

# BAM files
input_BAM_WT="${files_dir}/SRR7690004_WT.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"

YC_samples=(
    "YC1_MKDN240008444-1A_22CHMMLT4_L6.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YC2_MKDN240008445-1A_22CHMMLT4_L26.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YC3_MKDN240008446-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YC4_MKDN240008447-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YC5_MKDN240008448-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
)

YP_samples=(
    "YP1_MKDN240008449-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YP2_MKDN240008450-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YP3_MKDN240008451-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YP4_MKDN240008452-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
    "YP5_MKDN240008453-1A_22CHMMLT4_L5.bwa.cleanSAM.fixMate.markDup.fbayes.vcf.gatkBQSR.bam"
)

# Output files
Mutect2_YC="${output_dir}/YCpopulation_somatic_GATKMutect2.vcf.gz"
Mutect2_YP="${output_dir}/YPpopulation_somatic_GATKMutect2.vcf.gz"

# Log files
log_file_GATKMutect2_YC="${logfile_dir}/gatkMutect2YC.log"
log_file_GATKMutect2_YP="${logfile_dir}/gatkMutect2YP.log"

# Check input files
for bam in "$input_BAM_WT" "${YC_samples[@]/#/$files_dir/}" "${YP_samples[@]/#/$files_dir/}"
do
    if [[ ! -f "$bam" ]]; then
        echo "ERROR: Missing BAM file: $bam"
        exit 1
    fi
done

# Software code
reference="" # Reference genome
{
    echo "Running GATK Mutect2 for YC population samples..."
    gatk --java-options "-Xmx12g" Mutect2 -R "$reference" -I "$input_BAM_WT" -normal SRR7690004 \
    $(for bam in "${YC_samples[@]}"; do echo -I "${files_dir}/${bam}"; done) \
    -O "$Mutect2_YC"
    if [[ $? -eq 0 ]]; then
        echo "GATK Mutect2 SNP and INDEL detection completed successfully for YC population"
    else
        echo "GATK Mutect2 failed for YC population"
        exit 1  
    fi 
} >> "$log_file_GATKMutect2_YC" 2>&1

{
    echo "Running GATK Mutect2 for YP population samples..."
    gatk --java-options "-Xmx12g" Mutect2 -R "$reference" -I "$input_BAM_WT" -normal SRR7690004 \
    $(for bam in "${YP_samples[@]}"; do echo -I "${files_dir}/${bam}"; done) \
    -O "$Mutect2_YP"
    if [[ $? -eq 0 ]]; then
        echo "GATK Mutect2 SNP and INDEL detection completed successfully for YP population"
    else
        echo "GATK Mutect2 failed for YP population"
        exit 1  
    fi 
} >> "$log_file_GATKMutect2_YP" 2>&1

    ## GATK FilterMutectCalls: Filter variants in a Mutect2 VCF callset. It applies filters to the raw output of Mutect2. Parameters are contained in M2FiltersArgumentCollection and described in https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf

# Log files
log_file_GATK_FMC_YC="${logfile_dir}/gatkFilterMutectCalls_YC.log"
log_file_GATK_FMC_YP="${logfile_dir}/gatkFilterMutectCalls_YP.log"

# Output files
Mutect2_YC_filtered="${output_dir}/YCpopulation_somatic_GATKMutect2_filtered.vcf.gz"
Mutect2_YP_filtered="${output_dir}/YPpopulation_somatic_GATKMutect2_filtered.vcf.gz"

# Code
{
    echo "Running GATK FilterMutectCalls for YC population..."
    gatk FilterMutectCalls -R "$reference" -V "$Mutect2_YC" -O "$Mutect2_YC_filtered"
    if [[ $? -eq 0 ]]; then
        echo "GATK FilterMutectCalls completed successfully for YC population"
    else
        echo "GATK FilterMutectCalls failed for YC population"
        exit 1
    fi 
} >> "$log_file_GATK_FMC_YC" 2>&1

# Running GATK Mutect2 and FilterMutectCalls for YP population
{
    echo "Running GATK FilterMutectCalls for YP population..."
    gatk FilterMutectCalls -R "$reference" -V "$Mutect2_YP" -O "$Mutect2_YP_filtered"
    if [[ $? -eq 0 ]]; then
        echo "GATK FilterMutectCalls completed successfully for YP population"
    else
        echo "GATK FilterMutectCalls failed for YP population"
        exit 1
    fi
} >> "$log_file_GATK_FMC_YP" 2>&1

# Deactivate conda
conda deactivate

## After running this script and checking the log files, extract only the PASS variants to a new file (use grep '^#\|PASS' in the VCF file obtained after FilterMutectCalls). Then, run GATK VariantsToTable for extracting specified fields for each variant in a VCF file to a tab-delimited table. 
# zcat YCpopulation_somatic_GATKMutect2_filtered.vcf.gz | grep '^#\|PASS' > YCpopulation_somatic_GATKMutect2_filtered.PASS.vcf
# zcat YPpopulation_somatic_GATKMutect2_filtered.vcf.gz | grep '^#\|PASS' > YPpopulation_somatic_GATKMutect2_filtered.PASS.vcf
# gatk VariantsToTable -V "${workingdir}/YCpopulation_somatic_GATKMutect2_filtered.PASS.vcf" -F CHROM -F POS -F TYPE -F REF -F ALT -GF GT -GF AF -GF AD -GF DP -O "${workingdir}/YCpopulation_somatic.PASS.table"
# gatk VariantsToTable -V "${workingdir}/YPpopulation_somatic_GATKMutect2_filtered.PASS.vcf" -F CHROM -F POS -F TYPE -F REF -F ALT -GF GT -GF AF -GF AD -GF DP -O "${workingdir}/YPpopulation_somatic.PASS.table"
