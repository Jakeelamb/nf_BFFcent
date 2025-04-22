#!/bin/bash

# Check if directory argument was provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory_with_fastq_files>"
    echo "Example: $0 /scratch/alpine/soenksen@colostate.edu/variant_calling"
    exit 1
fi

# Directory containing the fastq files - remove trailing slash if present
FASTQ_DIR="${1%/}"
SAMPLE_LIST="samples.txt"

# Check if directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Directory $FASTQ_DIR does not exist"
    exit 1
fi

# Create header for samples.txt
echo -e "sample_id\tread1\tread2" > "$SAMPLE_LIST"

# Create a temporary file to store sample IDs and their corresponding R1 and R2 files
declare -A R1_FILES
declare -A R2_FILES

# Find all R1 and R2 files and store them by sample ID
for fastq_file in "$FASTQ_DIR"/*_R1.fastq.gz; do
    if [ ! -f "$fastq_file" ]; then
        echo "No R1 FASTQ files found in $FASTQ_DIR"
        exit 1
    fi
    
    # Extract sample ID from filename (remove _R1.fastq.gz)
    base_name=$(basename "$fastq_file")
    sample_id=${base_name%_R1.fastq.gz}
    
    # Store R1 file path
    R1_FILES["$sample_id"]="$fastq_file"
    
    # Construct and check for R2 file
    r2_file="${FASTQ_DIR}/$(basename "${fastq_file/_R1.fastq.gz/_R2.fastq.gz}")"
    if [ -f "$r2_file" ]; then
        R2_FILES["$sample_id"]="$r2_file"
    else
        echo "Warning: R2 file not found for sample $sample_id"
    fi
done

# Write samples with both R1 and R2 files to the sample list
for sample_id in "${!R1_FILES[@]}"; do
    if [ -n "${R2_FILES[$sample_id]}" ]; then
        echo -e "$sample_id\t${R1_FILES[$sample_id]}\t${R2_FILES[$sample_id]}" >> "$SAMPLE_LIST"
    fi
done

# Count samples and report
sample_count=$(grep -v "^sample_id" "$SAMPLE_LIST" | wc -l)
echo "Created $SAMPLE_LIST with $sample_count samples"
echo "Sample list location: $(readlink -f $SAMPLE_LIST)"

# Make the script executable
chmod +x "$0"


