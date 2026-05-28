#!/bin/bash

# TEcount allows to count reads that correspond to TEs and genes, for correlating with the expression levels. TEcount requirements: Python: 2.7.x or >= 3.2.x (tested on Python 2.7.11 and 3.7.7), pysam: 0.9.x or greater, R: 2.15.x or greater and DESeq2: 1.10.x or greater.

# Export PATH and PYTHONPATH
# Define path to .egg files
EGG_PYSAM="/home/ana/local/lib/python3.12/dist-packages/pysam-0.23.0-py3.12-linux-x86_64.egg"
EGG_TE="/home/ana/local/lib/python3.12/dist-packages/TEtranscripts-2.2.3-py3.12.egg"

export PYTHONPATH="${EGG_PYSAM}:${EGG_TE}:/home/ana/local/lib/python3.12/dist-packages/:$PYTHONPATH"

# Verify Python versions
python_version=$(python3 --version)
echo "Python version: $python_version" >> "$logfile_TEcount" 2>&1

# Verify pysam version
python3 -c "import pysam; print('¡Éxito total! Pysam encontrado por SLURM. Versión:', pysam.__version__)" >> "$logfile_TEcount" 2>&1

# Verify R version
R_version=$(R --version | head -n 1)
echo "R version: $R_version" >> "$logfile_TEcount" 2>&1

# Verify if DESeq2 is installed 
Rscript -e 'if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 no está instalado"); print(packageVersion("DESeq2"))' >> "$logfile_TEcount" 2>&1

# Create file variables
sample_list="/legserv/Archive/Ana/Sequencing_data/pdb_12M16E60E_Ana/sample_list_BAMRNAseqTC.txt"
line=$(head -n $SLURM_ARRAY_TASK_ID $sample_list| tail -1)
file_name="$(echo $line | awk '{print $1}' | tr -d '\r')"
sample_name="$(echo $file_name | awk -F "_" '{print $5}')"

# Create directory and file variables for the code
input_dir="/data/ana/Sequencing_Novogene/TEcount_RNAseqTC"
output_dir="/data/ana/Sequencing_Novogene/TEcount_RNAseqTC/output"
BAM_file="${input_dir}/${file_name}_Fusox3_hisat.bam"
gene_GTF="${input_dir}/Fusox3_1_GeneCatalog_20230106.gtf"
te_GTF="${input_dir}/TEs_Foxy3_newTElist_NoComments_Customized.gtf"
output_name="${sample_name}_TEcount"
logfile_TEcount="${output_dir}/TEcount_${sample_name}.log"

# Check if input files exist
{
if [[ ! -f "$BAM_file" || ! -f "$gene_GTF" || ! -f "$te_GTF" ]]; then
    echo "Error: One or more files are missing. Verify $BAM_file, $gene_GTF, o $te_GTF."
    exit 1
fi
} >> "$logfile_TEcount" 2>&1

# TEcount code
{
    echo "Running TEcount for ${file_name}_Fusox3_hisat.bam"
    python3 /home/ana/local/bin/TEcount --sortByPos --format BAM --mode multi -b "$BAM_file" --GTF "$gene_GTF" --TE "$te_GTF" --project "$output_name" --outdir "$output_dir"
    if [[ $? -eq 0 ]]; then
        echo "TEcount succesfully completed for ${file_name}"
    else
        echo "TEcount failed for ${file_name}"
        exit 1  
    fi 
} >> "$logfile_TEcount" 2>&1
