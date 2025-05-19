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
    # Decompress reference if needed
    if [[ "${ref_genome}" == *.gz ]]; then
        gunzip -c "${ref_genome}" > "${ref_base}"
        ref_to_use="${ref_base}"
    else
        ref_to_use="${ref_genome}"
    fi

    # Index the reference genome if it is not already indexed
    if [ ! -f "\${ref_to_use}.bwt" ]; then
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