process samtools {
    tag "$sample_id"
    publishDir "${params.outdir}/samtools", mode: params.publish_dir_mode
    cpus 24
    memory '32 GB'
    time '24h'

    input: 
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_processed.bam"), emit: processed_bam

    script:
    """
    # Create temporary directory for stats
    mkdir -p stats

    # Filter for mapped reads
    samtools view -@ ${task.cpus} -u \\
        -o "${sample_id}_justmap.bam" \\
        "${bam_file}"

    # Run samtools pipeline (collate, fixmate, sort, markdup)
    samtools collate -@ ${task.cpus} -O -u \\
        "${sample_id}_justmap.bam" - | \\
    samtools fixmate -@ ${task.cpus} -m -u - - | \\
    samtools sort -@ ${task.cpus} -u -O bam - | \\
    samtools markdup -@ ${task.cpus} \\
        -f "stats/${sample_id}_stats_file.txt" \\
        -S -d 2500 \\
        --mode t \\
        --include-fails \\
        - "${sample_id}_processed.bam"

    # Index the final BAM file
    samtools index -@ ${task.cpus} "${sample_id}_processed.bam"

    # Clean up intermediate files
    rm "${sample_id}_justmap.bam"
    """
}