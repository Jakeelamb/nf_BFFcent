#!/bin/bash

# Script to create a GATK map file for GenomicsDBImport from a directory of GVCF files
# Usage: ./create_gatk_map_file.sh <gvcf_directory> <interval> [output_file]

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <gvcf_directory> <interval> [output_file]"
    echo "Example: $0 ./results/gatk_haplotype 1 interval_1_map.txt"
    echo "If output_file is not specified, it will default to <interval>_sample_map.txt"
    exit 1
fi

# Directory containing the GVCF files
GVCF_DIR="$1"
INTERVAL="$2"
OUTPUT_FILE="${3:-${INTERVAL}_sample_map.txt}"
CONCAT_DIR="${GVCF_DIR}/concatenated"

# Check if directory exists
if [ ! -d "$GVCF_DIR" ]; then
    echo "Error: Directory $GVCF_DIR does not exist"
    exit 1
fi

# Create directory for concatenated files
mkdir -p "$CONCAT_DIR"

# Get unique sample names from all GVCF files
echo "Identifying unique samples from GVCF files..."
declare -A samples
for gvcf_file in "$GVCF_DIR"/*.g.vcf.gz; do
    if [ -f "$gvcf_file" ]; then
        base_name=$(basename "$gvcf_file")
        # Extract sample name (everything before the last dot before .g.vcf.gz)
        sample_id=${base_name%.*.g.vcf.gz}
        samples["$sample_id"]=1
    fi
done

# Function to extract chromosome number from filename
get_chrom_number() {
    local filename=$1
    # Extract the number between the last dot and .g.vcf.gz
    local chrom_num=$(echo "$filename" | sed -E 's/.*\.([0-9]+)\.g\.vcf\.gz/\1/')
    # Convert to integer for proper numeric sorting
    echo $((10#$chrom_num))
}

# Concatenate GVCF files for each sample
echo "Concatenating GVCF files for each sample..."
for sample_id in "${!samples[@]}"; do
    echo "Processing sample: $sample_id"
    
    # Find all GVCF files for this sample and sort them by chromosome number
    declare -a sorted_files
    while IFS= read -r file; do
        sorted_files+=("$file")
    done < <(find "$GVCF_DIR" -name "${sample_id}.*.g.vcf.gz" | sort -t. -k2,2n)
    
    if [ ${#sorted_files[@]} -gt 0 ]; then
        echo "Found ${#sorted_files[@]} chromosome files for sample $sample_id"
        
        # Create a file list for bcftools concat
        file_list="${CONCAT_DIR}/${sample_id}_file_list.txt"
        printf "%s\n" "${sorted_files[@]}" > "$file_list"
        
        # Verify the order of files
        echo "Files to be concatenated in order:"
        while IFS= read -r file; do
            chrom_num=$(get_chrom_number "$file")
            echo "Chromosome $chrom_num: $file"
        done < "$file_list"
        
        # Concatenate the files
        echo "Concatenating files for sample $sample_id..."
        bcftools concat -f "$file_list" -O z -o "${CONCAT_DIR}/${sample_id}.g.vcf.gz"
        bcftools index "${CONCAT_DIR}/${sample_id}.g.vcf.gz"
        
        # Clean up the file list
        rm "$file_list"
    else
        echo "Warning: No GVCF files found for sample $sample_id"
    fi
done

# Initialize counter and map file
count=0
> "$OUTPUT_FILE"

# Create map file from concatenated GVCFs
echo "Creating map file from concatenated GVCFs..."
for concat_file in "$CONCAT_DIR"/*.g.vcf.gz; do
    if [ -f "$concat_file" ]; then
        base_name=$(basename "$concat_file")
        sample_id=${base_name%%.g.vcf.gz}
        
        # Write sample ID and concatenated GVCF path to map file (tab-delimited)
        echo -e "$sample_id\t$concat_file" >> "$OUTPUT_FILE"
        count=$((count + 1))
    fi
done

if [ $count -eq 0 ]; then
    echo "Warning: No concatenated GVCF files found in $CONCAT_DIR"
    exit 1
else
    echo "Created $OUTPUT_FILE with $count samples for interval $INTERVAL"
    echo "Map file location: $(readlink -f $OUTPUT_FILE)"
fi

# Make the script executable
chmod +x "$0"

