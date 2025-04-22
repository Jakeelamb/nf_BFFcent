process REFERENCE_INDEX {
    publishDir "$params.outdir/ref_index", mode: params.publish_dir_mode
    cpus 4
    memory '128 GB'
    time '24h'

    input:
    tuple path(ref_genome), path(ref_genome_dict), path(ref_genome_fasta_index)

    output:
    tuple path(ref_genome), path(ref_genome_dict), path(ref_genome_fasta_index), emit: indexed_ref

    script:
    def ref_base = ref_genome.toString().replaceFirst(/\.gz$/, '')
    """
    # Check if all required index files exist
    missing_files=0

    # Get uncompressed reference name for BWA-MEM2 index checking
    ref_for_bwa="\${PWD}/${ref_base}"
    orig_ref="\${PWD}/${ref_genome}"
    echo "Reference genome: \$orig_ref"
    echo "Reference for BWA indexes: \$ref_for_bwa"

    # Check BWA-MEM2 indexes
    for ext in 'amb' 'ann' 'bwt.2bit.64' 'pac' 'sa' '0123'; do
        # Index files won't have .gz extension even if reference is compressed
        if [ ! -f "\${ref_for_bwa}.\${ext}" ] && [ ! -f "${ref_genome}.\${ext}" ]; then
            missing_files=\$((missing_files+1))
            echo "Missing BWA-MEM2 index file: \${ref_for_bwa}.\${ext} or ${ref_genome}.\${ext}"
        else
            echo "Found BWA-MEM2 index: \${ref_for_bwa}.\${ext} or ${ref_genome}.\${ext}"
        fi
    done

    # Check dictionary
    if [ ! -s "${ref_genome_dict}" ]; then
        missing_files=\$((missing_files+1))
        echo "Missing reference dictionary: ${ref_genome_dict}"
    else
        echo "Found dictionary: ${ref_genome_dict}"
    fi

    # Check FASTA index
    if [ ! -s "${ref_genome_fasta_index}" ]; then
        missing_files=\$((missing_files+1))
        echo "Missing FASTA index: ${ref_genome_fasta_index}"
    else
        echo "Found FASTA index: ${ref_genome_fasta_index}"
    fi

    # If all indexes exist, do nothing
    if [ \$missing_files -eq 0 ]; then
        echo "All reference index files already exist. Skipping indexing."
        exit 0
    fi

    # Otherwise, create the missing indexes
    echo "Need to create \$missing_files index files."

    # Decompress reference if needed
    if [[ "${ref_genome}" == *.gz ]]; then
        echo "Decompressing reference genome..."
        gunzip -c "${ref_genome}" > "${ref_base}"
        ref_to_use="${ref_base}"
    else
        ref_to_use="${ref_genome}"
    fi

    # Index the reference genome if it is not already indexed
    if [ ! -f "\${ref_to_use}.bwt" ] || [ ! -f "\${ref_to_use}.bwt.2bit.64" ]; then
        echo "Indexing reference genome with BWA-MEM2..."
        bwa-mem2 index "\${ref_to_use}"
    else
        echo "BWA-MEM2 index already exists, skipping..."
    fi

    # Create GATK sequence dictionary if it is not already created
    if [ ! -s "${ref_genome_dict}" ]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R "\${ref_to_use}" -O "${ref_genome_dict}"
    else
        echo "Sequence dictionary already exists, skipping..."
    fi

    # Create reference fasta index if it is not already created
    if [ ! -s "${ref_genome_fasta_index}" ]; then
        echo "Creating FASTA index..."
        samtools faidx "\${ref_to_use}"
        mv "\${ref_to_use}.fai" "${ref_genome_fasta_index}"
    else
        echo "FASTA index already exists, skipping..."
    fi
    """
}