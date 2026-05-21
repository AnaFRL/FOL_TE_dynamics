## Author = Ana Rodríguez López
## Date.of.creation = 11.10.24

# bbtools instalation with mamba (conda)
# conda create -n bbtools_env
# conda activate bbtools_env
# mamba install bioconda::bbmap
# All the info about this tool can be reached at:
# /path/directory/instalation  # You can see where it is installed in the MAMBA installation output


#SCRIPT

##sbatch --array=numb_samples_to_processed%X -c 4 --mem=32G myscript_name.sh

#Initialize conda 
source /home/ana/anaconda3/etc/profile.d/conda.sh

##conda activate bbtools_env (here I have installed all the programs I need from bbmap)
conda activate bbtools_env

# REFORMAT.SH FROM BBMAP. 
# Preprocessing: merge forward and reverse reads in one file (interleaved files).

#!/bin/bash
#Create variables with directories
sample_list="path" #A text list with 3 columns: 1- Sample name, 2- Forward reads path, 3- Reverse reads path
output_interleaved="output/dir" #Where to allocate the output files
mkdir -p $output_interleaved
#Create variables for iterate
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name=$(echo $line | awk '{print $1}')
forward=$(echo $line | awk '{print $2}')
reverse=$(echo $line | awk '{print $3}')
interleaved=${output_interleaved}/${file_name}_interleaved.fq.gz
#Create log file variable
log_file=${output_interleaved}/${file_name}_interleaved.log
#Software codes
{
echo "Processing ${file_name}..."  
reformat.sh in1="$forward" in2="$reverse" out="$interleaved"
if [[ $? -eq 0 ]]; then
    echo "Interleaving completed successfully for ${file_name}."
else
    echo "Interleaving failed for ${file_name}."
    exit 1
fi
} >> "$log_file" 2>&1  # Redirects sstdout and stserr to a common file

# TESTFORMAT2 FROM BBMAP. Reads the entire file to find extended information about the format and contents.
#Create variables 
interleaved_dir="out/dir/with/interleaved/reads"
input_interleaved=${interleaved_dir}/${file_name}_interleaved.fq.gz
log_file_TF=${interleaved_dir}/interleaved_testformat2_${file_name}.log
# Verify that the input file exists
if [[ ! -f "$input_interleaved" ]]; then
    echo "Error: Input file ${input_interleaved} does not exist." >&2
    exit 1
fi
#Sotware execution
{
    echo "Processing ${file_name}..."  
testformat2.sh $input_interleaved qhist=${interleaved_dir}/${file_name}_interleaved_qhist.txt ihist=${interleaved_dir}/${file_name}_interleaved_ihist.txt bhist=${interleaved_dir}/${file_name}_interleaved_bhist.txt lhist=${interleaved_dir}/${file_name}_interleaved_lhist.txt gchist=${interleaved_dir}/${file_name}_interleaved_gchist.txt zmwhist=null khist=null merge=t trim=t speed=t
if [[ $? -eq 0 ]]; then
        echo "testformat2 completed successfully for ${file_name}_interleaved.fq.gz."
    else
        echo "testformat2 failed for ${file_name}_interleaved.fq.gz."
        exit 1
    fi
} >> "$log_file_TF" 2>&1


# BBDUK FROM BBMAP. It combines most common data-quality-related trimming, filtering, and masking operations
    ##ADAPTER TRIMMING first, for removing adapter sequences and do a quality filtering
input_reads=${interleaved_dir}/${file_name}_interleaved.fq.gz
output_reads=${interleaved_dir}/${file_name}_interleaved_ADAPTERTRIM.fq.gz
adapters_file="/bbmap-39.10-0/resources/adapters.fa" #you should find where is this file in the computer
log_file_BBD=${interleaved_dir}/interleaved_bbduck_${file_name}.log
# Verify that the input file exists
if [[ ! -f "$input_reads" ]]; then
    echo "Error: Input file for bbduk adapter and quality trimming ${input_reads} does not exist." >&2
    exit 1
fi
#Sotware execution
{
    echo "Processing ${file_name}_interleaved..."  
    bbduk.sh -Xmx6g t=4 in=$input_reads out=$output_reads ref=$adapters_file interleaved=true ktrim=r k=23 mink=11 hdist=1 tpe tbo entropy=0.1 filterpolyg=10 qtrim=r trimq=10 maq=10 ordered=t ftm=5
    if [[ $? -eq 0 ]]; then
        echo "bbduk adapter and quality trimming completed successfully for ${file_name}_interleaved.fq.gz."
    else
        echo "bbduk adapter and quality trimming failed for ${file_name}_interleaved.fq.gz."
        exit 1  
    fi
} >> "$log_file_BBD" 2>&1

   ##CONTAMINANT FILTERING, for removing phiX sequences and other possible contaminants
# Variables
input_trim_files=$output_reads
unmatched_output=${interleaved_dir}/${file_name}_interleaved_DECONTAMINATED.fq.gz
matched_output=${interleaved_dir}/${file_name}_interleaved_PHIXMATCHED.fq.gz
contaminants_file="/home/ana/anaconda3/envs/bbtools_env/opt/bbmap-39.10-0/resources/phix174_ill.ref.fa.gz"
stats_file=${interleaved_dir}/${file_name}_interleaved_decontamination_stats.txt
# Verify that the input file exists
if [[ ! -f "$input_trim_files" ]]; then
    echo "Error: Input file for bbduk decontamination ${input_trim_files} does not exist." >&2
    exit 1
fi
#Sotware execution
{
    echo "Processing ${file_name}_interleaved_ADAPTERTRIM.fq.gz..."  
    bbduk.sh -Xmx6g t=4 in=$input_trim_files out=$unmatched_output interleaved=true outm=$matched_output ref=$contaminants_file k=31 hdist=1 stats=$stats_file ordered=t
    if [[ $? -eq 0 ]]; then
        echo "bbduk decontamination completed successfully for ${file_name}_interleaved_ADAPTERTRIM.fq.gz."
    else
        echo "bbduk decontamination failed for ${file_name}_interleaved_ADAPTERTRIM.fq.gz."
        exit 1  
    fi
} >> "$log_file_BBD" 2>&1


#ERROR CORRECTION with TADPOLE + TESTFORMAT2 (BBMAP).
#Tadpole: Uses kmer counts to assemble contigs, extend sequences, or error-correct reads --> It uses a lot of memory!!

#Log file variable
log_file_TPTF=${interleaved_dir}/interleaved_TadpoleTestformat2_${file_name}.log

    ##ERROR CORRECTION WITH TADPOLE
input_tadpole=$unmatched_output #output from bbduk decontamination step
output_tadpole=${interleaved_dir}/${file_name}_interleaved_CORRECTED.fq.gz
# Verify that the input file exists
if [[ ! -f "$input_tadpole" ]]; then
    echo "Error: Input file for tadpole error correction ${input_tadpole} does not exist." >&2
    exit 1
fi
#Sotware execution
{
    echo "Processing ${file_name}_interleaved_DECONTAMINATED.fq.gz..."  
    tadpole.sh in=$input_tadpole out=$output_tadpole mode=correct k=50 prefilter=2 t=4 
    if [[ $? -eq 0 ]]; then
        echo "tadpole error correction completed successfully for ${file_name}_interleaved_DECONTAMINATED.fq.gz"
    else
        echo "tadpole error correction failed for ${file_name}_interleaved_DECONTAMINATED.fq.gz"
        exit 1  
    fi
} >> "$log_file_TPTF" 2>&1

    ##TESTFORMAT2
# Verify that the input file exists
if [[ ! -f "$output_tadpole" ]]; then
    echo "Error: Input file ${output_tadpole} does not exist." >&2
    exit 1
fi
#Sotware execution
{
    echo "Processing ${file_name}_interleaved_CORRECTED.fq.gz..."  
testformat2.sh $output_tadpole qhist=${interleaved_dir}/${file_name}_interleaved_corrected_qhist.txt ihist=${interleaved_dir}/${file_name}_interleaved_corrected_ihist.txt bhist=${interleaved_dir}/${file_name}_interleaved_corrected_bhist.txt lhist=${interleaved_dir}/${file_name}_interleaved_corrected_lhist.txt gchist=${interleaved_dir}/${file_name}_interleaved_corrected_gchist.txt zmwhist=null khist=null merge=t trim=t speed=t
if [[ $? -eq 0 ]]; then
        echo "testformat2 completed successfully for ${file_name}_interleaved_CORRECTED.fq.gz."
    else
        echo "testformat2 failed for ${file_name}_interleaved_CORRECTED.fq.gz."
    fi
} >> "$log_file_TPTF" 2>&1


# DEINTERLEAVE READS with REFORMAT (bbmap) and FASTQC

# Log file variable
log_file_RefFQ=${interleaved_dir}/reformat_deinterleaveFastQC_${file_name}.log

    ##REFORMAT: Deinterleave reads
# Create variables for iterate
interleaved_corrected_files=$output_tadpole
output_read1=${interleaved_dir}/${file_name}_1_corrected.fq.gz
output_read2=${interleaved_dir}/${file_name}_2_corrected.fq.gz
# Verify that the input file exists
if [[ ! -f "$interleaved_corrected_files" ]]; then
    echo "Error: Input file ${interleaved_corrected_files} does not exist." >&2
    exit 1
fi
# Software execution
{
    echo "Processing ${file_name}_interleaved_CORRECTED.fq.gz..."
    reformat.sh in=$interleaved_corrected_files out1=$output_read1 out2=$output_read2
    if [[ $? -eq 0 ]]; then
        echo "Reformat deinterleave completed successfully for ${file_name}_interleaved_CORRECTED.fq.gz"
    else
        echo "Reformat deinterleave failed for ${file_name}_interleaved_CORRECTED.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_RefFQ" 2>&1

#Deactivate bbtools_env
conda deactivate

    ## FASTQC ANALYSIS OF THE CLEAN READS
# Activate conda env
conda activate fastqc_env
# Verify that the input file exists
if [[ ! -f "$output_read1" ]]; then
    echo "Error: Input file ${output_read1} does not exist." >&2
    exit 1
fi
if [[ ! -f "$output_read2" ]]; then
    echo "Error: Input file ${output_read2} does not exist." >&2
    exit 1
fi
# Variables
fastqc_corrected_dir=${interleaved_dir}/fastQC_bbmapCorrectedReads
mkdir -p $fastqc_corrected_dir
# Command
{
    echo "FastQC processing of ${file_name}_1_corrected.fq.gz..."
    fastqc $output_read1 -o $fastqc_corrected_dir
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for ${file_name}_1_corrected.fq.gz"
    else
        echo "FastQC failed for ${file_name}_1_corrected.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_RefFQ" 2>&1
{
    echo "FastQC processing of ${file_name}_2_corrected.fq.gz..."
    fastqc $output_read2 -o $fastqc_corrected_dir
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for ${file_name}_2_corrected.fq.gz"
    else
        echo "FastQC failed for ${file_name}_2_corrected.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_RefFQ" 2>&1

#Deactivate conda environment
conda deactivate


# TRIMMOMATIC + FASTQC

#Create log file variable 
log_file_trimmomatic_FQ=${interleaved_dir}/trimmomatic_fastqc_${file_name}.log

    ##TRIMMOMATIC
#Activate conda environment
conda activate trimmomatic_env
#Input reads
Read1_input=${interleaved_dir}/${file_name}_1_corrected.fq.gz
Read2_input=${interleaved_dir}/${file_name}_2_corrected.fq.gz
#Output paired reads
Read1_output_paired=${interleaved_dir}/${file_name}_1_cleanP.fq.gz
Read2_output_paired=${interleaved_dir}/${file_name}_2_cleanP.fq.gz
#Output unpaired reads
Read1_output_unpaired=${interleaved_dir}/${file_name}_1_clean_UP.fq.gz
Read2_output_unpaired=${interleaved_dir}/${file_name}_2_clean_UP.fq.gz
#Adapter file
adapters_file_trimm="/home/ana/anaconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/adapters/TruSeq-PE_adapters.fa" # Merged of all Trimmomatic defautl PE-adapter files
#Summary stats file
summary_file=${interleaved_dir}/${file_name}_trimmommatic_sumfile.txt
#Software code
# Verify that the input file exists
if [[ ! -f "$Read1_input" ]]; then
    echo "Error: Input file ${Read1_input} does not exist." >&2
    exit 1
fi
if [[ ! -f "$Read2_input" ]]; then
    echo "Error: Input file ${Read2_input} does not exist." >&2
    exit 1
fi
{
    echo "Trimmomatic running for ${file_name}_1_corrected.fq.gz and ${file_name}_2_corrected.fq.gz"
    trimmomatic PE -threads 4 -phred33 -summary $summary_file $Read1_input $Read2_input $Read1_output_paired $Read1_output_unpaired $Read2_output_paired $Read2_output_unpaired ILLUMINACLIP:$adapters_file_trimm:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 
    if [[ $? -eq 0 ]]; then
        echo "Trimmomatic completed successfully for ${file_name}_1_corrected.fq.gz and ${file_name}_2_corrected.fq.gz"
    else
        echo "Trimmomatic failed for ${file_name}_1_corrected.fq.gz and ${file_name}_2_corrected.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_trimmomatic_FQ" 2>&1
    
#Deactivate conda environment for trimmomatic
conda deactivate
    
    ##FASTQC
#Activate conda environment
conda activate fastqc_env
#Create directory
fastqc_clean_dir=${interleaved_dir}/fastQC_trimmomaticCleanReads
mkdir -p $fastqc_clean_dir
# Verify that the input file exists
if [[ ! -f "$Read1_output_paired" ]]; then
    echo "Error: Input file ${Read1_output_paired} does not exist." >&2
    exit 1
fi

if [[ ! -f "$Read2_output_paired" ]]; then
    echo "Error: Input file ${Read2_output_paired} does not exist." >&2
    exit 1
fi
#Command line
{
    echo "FastQC processing of ${file_name}_1_cleanP.fq.gz"
    fastqc $Read1_output_paired -o $fastqc_clean_dir
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for ${file_name}_1_cleanP.fq.gz"
    else
        echo "FastQC failed for ${file_name}_1_cleanP.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_trimmomatic_FQ" 2>&1
{
    echo "FastQC processing of ${file_name}_2_cleanP.fq.gz"
    fastqc $Read2_output_paired -o $fastqc_clean_dir
    if [[ $? -eq 0 ]]; then
        echo "FastQC completed successfully for ${file_name}_2_cleanP.fq.gz"
    else
        echo "FastQC failed for ${file_name}_2_cleanP.fq.gz"
        exit 1  # Exit the script if it fails
    fi
} >> "$log_file_trimmomatic_FQ" 2>&1
#Deactivate conda
conda deactivate
