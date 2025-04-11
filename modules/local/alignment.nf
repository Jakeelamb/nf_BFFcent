process alignment {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode
    cpus 24
    memory '32 GB'
    time '24h'

    input:
    tuple val(sample_id), path(trimmed_reads), path(ref_genome)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_files
    path "${sample_id}_readgroup_info.txt", emit: readgroup_info

    script:
    """
    # Function to extract read group information
    extract_read_groups() {
        local fastq_file=\$1
        local sample_name=\$2
        
        # Extract the first header line
        local header=\$(zcat \${fastq_file} | head -n 1)
        
        # Extract components
        local flowcell=\$(echo \${header} | awk -F: '{print \$3}')
        local lane=\$(echo \${header} | awk -F: '{print \$4}')
        local barcode=\$(echo \${header} | awk -F: '{print \$10}')
        
        # Hardcoded values
        local LB="Lib2"
        local PL="Illumina"
        
        # Construct ID and PU
        local ID="\${sample_name}.\${LB}.\${flowcell}.\${lane}"
        local PU="\${flowcell}.\${lane}.\${barcode}"
        
        # Save info to file
        printf "%-15s\\t%-30s\\t%-15s\\t%-10s\\t%-30s\\t%-10s\\t%-50s\\n" \\
            "\${sample_name}" "\${ID}" "\${sample_name}" "\${LB}" "\${PU}" "\${PL}" "\${header}" > ${sample_id}_readgroup_info.txt
        
        # Also return the formatted RG string for BWA-MEM2
        echo "@RG\\tID:\${ID}\\tSM:\${sample_name}\\tLB:\${LB}\\tPU:\${PU}\\tPL:\${PL}"
    }
    
    # Print header for the readgroup table
    printf "%-15s\\t%-30s\\t%-15s\\t%-10s\\t%-30s\\t%-10s\\t%-50s\\n" \\
        "Sample" "ID" "SM" "LB" "PU" "PL" "Raw Header" > ${sample_id}_readgroup_info.txt
    printf "%s\\n" "\$(printf '=%.0s' {1..160})" >> ${sample_id}_readgroup_info.txt
    
    # Extract read group information
    RG_STRING=\$(extract_read_groups "${trimmed_reads[0]}" "${sample_id}")
    
    # Validate read group information
    if grep -P "\\t\\t" ${sample_id}_readgroup_info.txt > /dev/null; then
        echo "WARNING: Read group information has missing values. Please review the output."
        exit 1
    else
        echo "Read group information validated successfully."
    fi
    
    # Run BWA-MEM2 alignment
    echo "Running BWA-MEM2 alignment..."
    bwa-mem2 mem -t ${task.cpus} \\
        -R "\${RG_STRING}" \\
        "${ref_genome}" \\
        "${trimmed_reads[0]}" \\
        "${trimmed_reads[1]}" \\
        > "${sample_id}.sam"
    
    # Convert SAM to BAM
    samtools view -@ ${task.cpus} -b -o "${sample_id}.bam" "${sample_id}.sam"
    
    # Remove SAM file to save space
    rm "${sample_id}.sam"
    """
}