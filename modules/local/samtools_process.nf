process SAMTOOLS_PROCESS {
    tag "$meta.id"
    publishDir "${params.outdir}/samtools", mode: params.publish_dir_mode
    cpus 24
    memory '32 GB'
    time '24h'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}.sorted.markdup.bam"), emit: processed_bam
    tuple val(meta), path("${meta.id}.sorted.markdup.bam.bai"), emit: bam_index
    path "${meta.id}_samtools_stats.txt", emit: stats

    script:
    """
    # First sort by queryname for fixmate
    samtools sort -n -@ ${task.cpus} -o ${meta.id}.namesorted.bam ${bam_file}
    
    # Run fixmate to add mate score tag
    samtools fixmate -m ${meta.id}.namesorted.bam ${meta.id}.fixmate.bam
    
    # Sort by coordinate for markdup
    samtools sort -@ ${task.cpus} -o ${meta.id}.positionsorted.bam ${meta.id}.fixmate.bam
    
    # Mark duplicates
    samtools markdup -@ ${task.cpus} ${meta.id}.positionsorted.bam ${meta.id}.sorted.markdup.bam
    
    # Index the final BAM file
    samtools index -@ ${task.cpus} ${meta.id}.sorted.markdup.bam
    
    # Generate statistics
    samtools stats ${meta.id}.sorted.markdup.bam > ${meta.id}_samtools_stats.txt
    
    # Clean up intermediate files
    rm ${meta.id}.namesorted.bam ${meta.id}.fixmate.bam ${meta.id}.positionsorted.bam
    """
}