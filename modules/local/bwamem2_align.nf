process BWAMEM2_ALIGN {
    tag "$meta.id"
    label 'process_high_memory' // Adjust label if needed based on your setup

    // Use conda environment defined in config
    conda params.conda_env_path

    input:
    tuple val(meta), path(reads)
    // Reference genome parameter implicitly available via params.reference_genome

    output:
    // Corrected: Emit the output channel as 'bam_files' as expected by main.nf
    tuple val(meta), path("*.aligned.bam"), emit: bam_files
    path "workdir_contents_before.txt"      , emit: log_before // Optional log output
    path "workdir_contents_after.txt"       , emit: log_after  // Optional log output

    script:
    // Convert the string path parameter into a Path object
    def ref_path_obj = file(params.reference_genome)
    def ref_dir = ref_path_obj.parent
    def ref_name = ref_path_obj.name
    def ref_path_str = params.reference_genome // Keep original string for shell script use

    """
    #!/bin/bash
    set -e -o pipefail

    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt

    echo "Using reference path: ${ref_path_str}"
    echo "Reference Dir: ${ref_dir}"
    echo "Reference Name: ${ref_name}"

    echo "Creating symlinks for BWA-MEM2 index files..."
    # UPDATED: Added '0123' to the list of extensions to be linked
    for ext in amb ann bwt.2bit.64 pac sa 0123; do
        src_file="${ref_dir}/${ref_name}.\$ext"
        tgt_link="genome.fa.\$ext"
        if [ -e "\$src_file" ]; then
            ln -s "\$src_file" "\$tgt_link"
            if [ ! -e "\$tgt_link" ]; then
                echo "Error: Failed to create symlink \$tgt_link from \$src_file"
                exit 1
            fi
        else
            # Ensure all required index files exist in the source directory
            echo "Error: Source index file \$src_file does not exist. Required extension: \$ext"
            exit 1
        fi
    done

    # Check for additional optional BWA index file (.bwt)
    # This might be redundant if .bwt.2bit.64 is always present and used by bwa-mem2,
    # but keep it for now unless it causes issues.
    bwt_src="${ref_dir}/${ref_name}.bwt"
    bwt_tgt="genome.fa.bwt"
    if [ -e "\$bwt_src" ]; then
        # Avoid creating genome.fa.bwt if genome.fa.bwt.2bit.64 exists, as bwa-mem2 uses the latter
        if [ ! -e "genome.fa.bwt.2bit.64" ]; then
            ln -s "\$bwt_src" "\$bwt_tgt"
            echo "Created symlink for \$bwt_tgt"
        else
            echo "Skipping redundant symlink for optional .bwt file as .bwt.2bit.64 exists."
        fi
    fi

    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt

    # Function to extract read group information
    get_read_group() {
        local fastq=\$1
        local sample_name=\$2

        # Extract header information using the read group info from the previous stage
        local rg_file="${meta.id}_rg_info.txt"
        if [ -f "\$rg_file" ]; then
            cat "\$rg_file"
        else
            echo "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:lib1\\tPL:ILLUMINA"
        fi
    }

    # Get read group information
    RG_STRING=\$(get_read_group "${reads[0]}" "${meta.id}")
    echo "Constructed RG string: \$RG_STRING"

    echo "Running BWA-MEM2 alignment..."
    # Pass the shell variable RG_STRING to -R, ensure it's quoted
    bwa-mem2 mem -t ${task.cpus} -M -R "\$RG_STRING" genome.fa ${reads[0]} ${reads[1]} | \\
    samtools view -@ ${task.cpus-1} -bS - > ${meta.id}.bam

    echo "Moving output BAM..."
    # Ensure the output matches the defined output name pattern
    mv ${meta.id}.bam ${meta.id}.aligned.bam

    echo "Alignment finished."
    """
}