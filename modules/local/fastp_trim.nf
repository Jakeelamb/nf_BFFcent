process FASTP_TRIM {
    tag "$meta.id"
    publishDir "${params.outdir}/trimming", mode: params.publish_dir_mode
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "${meta.id}_fastp.{html,json}", emit: report

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${meta.id}_trimmed_R1.fastq.gz \\
        -O ${meta.id}_trimmed_R2.fastq.gz \\
        --html ${meta.id}_fastp.html \\
        --json ${meta.id}_fastp.json \\
        --detect_adapter_for_pe \\
        -w ${task.cpus} \\
        -q 20 \\
        --compression 6 \\
	--correction
    """
} 