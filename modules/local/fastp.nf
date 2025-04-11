process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimming", mode: params.publish_dir_mode
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.{html,json}", emit: report

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}_trimmed_R1.fastq.gz \\
        -O ${sample_id}_trimmed_R2.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --detect_adapter_for_pe \\
        --dedup \\
        -w ${task.cpus} \\
        -q 20 \\
        --compression 6 \\
	--correction
    """
} 