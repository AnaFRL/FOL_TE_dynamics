#!/bin/bash

## Customize repeatmasker output, adding TE information, for downstream use with TEcount / TEtranscripts (https://github.com/mhammell-laboratory/TEtranscripts)

# Files
directory="/data/ana/Sequencing_Novogene/TEcount_RNAseqTC"
gtf_file="${directory}/TEs_Foxy3_newTElist_NoComments.gtf" # TE GTF without comments '##'. Use grep -v "^##" file.gtf > file_NoComments.gtf
te_info="${directory}/TE_family_class_Foxy3.txt" #If you have created this in Excel or windows apps, first run this: sed -i 's/\r$//' /data/ana/Sequencing_Novogene/TEcount_RNAseqTC/TE_family_class_Foxy3.txt
# Output file
output_file="${directory}/TEs_Foxy3_newTElist_NoComments_Customized.gtf"

# Initialize counter for TE copies
declare -A count

# Process GTF file line by line
while IFS=$'\t' read -r chr source type start end score strand phase attributes; do

    # Extract TE name from the Target field
    TE_name=$(echo "$attributes" | grep -oP 'Target "Motif:\K[^"]+')
    
    # Assign "Unknown" if TE name is missing
    if [[ -z "$TE_name" ]]; then
        TE_name="Unknown"
    fi

    # Retrieve TE family and class information
    read -r family class <<< $(grep -w "$TE_name" "$te_info" | awk '{print $2, $3}')

    # Assign "Unknown" if family or class are missing
    if [[ -z "$family" ]]; then
        family="Unknown"
    fi
    if [[ -z "$class" ]]; then
        class="Unknown"
    fi

    # Generate unique identifier for each TE copy
    count["$TE_name"]=$((count["$TE_name"] + 1))
    TE_id="TE_${TE_name}_dup_${count["$TE_name"]}"

    # Write modified entry to output file
    echo -e "$chr\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\tgene_id \"$TE_name\"; transcript_id \"$TE_id\"; family_id \"$family\"; class_id \"$class\";" >> "$output_file"
done < "$gtf_file"  # Procesar todo el archivo

echo "Output file generated: $output_file"
