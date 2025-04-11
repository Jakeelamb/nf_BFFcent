process GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/gatk_haplotype", mode: params.publish_dir_mode
    cpus 28
    memory '64 GB'
    time '24h'

    input:
    tuple val(sample_id), path(bam_file), path(ref_genome)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcfs

    script:
    """
    # Create temporary directory
    mkdir -p temp
    
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=temp" HaplotypeCaller \\
        -R "${ref_genome}" \\
        -I "${bam_file}" \\
        -O "${sample_id}.g.vcf.gz" \\
        --native-pair-hmm-threads ${task.cpus} \\
        -ERC GVCF

    # Index the GVCF file
    gatk IndexFeatureFile -I "${sample_id}.g.vcf.gz"
    """
}