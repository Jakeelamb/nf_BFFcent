process GATK_HAPLOTYPECALLER {
    tag "$meta.id - Interval $interval"
    publishDir "${params.outdir}/gatk_haplotype", mode: params.publish_dir_mode, pattern: "*.g.vcf.gz*"
    cpus 28
    memory '64 GB'
    time '24h'

    input:
    tuple val(meta), path(bam_file), val(interval), path(reference, stageAs: "genome.fa")

    output:
    tuple val(meta), val(interval), path("${meta.id}.${interval}.g.vcf.gz"), path("${meta.id}.${interval}.g.vcf.gz.tbi"), emit: gvcfs_interval

    script:
    """
    # List current directory for debugging
    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    # Define reference path variables
    REF_PATH="${params.reference_genome}"
    REF_DIR=\$(dirname "\$REF_PATH")
    REF_NAME=\$(basename "\$REF_PATH")
    REF_BASE=\$(echo "\$REF_NAME" | sed 's/\\.[^.]*\$//')
    
    # Create symlink for .fai index
    ln -s "\${REF_PATH}.fai" "genome.fa.fai"
    if [ ! -e "genome.fa.fai" ]; then
        echo "Error: genome.fa.fai not found. Check if \${REF_PATH}.fai exists."
        exit 1
    fi
    
    # Create symlink for .dict file
    ln -s "\${REF_DIR}/\${REF_BASE}.dict" "genome.dict"
    if [ ! -e "genome.dict" ]; then
        echo "Error: genome.dict not found. Check if \${REF_DIR}/\${REF_BASE}.dict exists."
        exit 1
    fi
    
    # List current directory after linking for debugging
    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    # Create temporary directory for GATK within task directory
    mkdir -p gatk_tmp
    
    # Run HaplotypeCaller
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./gatk_tmp" HaplotypeCaller \\
        -R genome.fa \\
        -I ${bam_file} \\
        -O ${meta.id}.${interval}.g.vcf.gz \\
        -L ${interval} \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        --smith-waterman JAVA

    # Index the GVCF file
    gatk IndexFeatureFile -I "${meta.id}.${interval}.g.vcf.gz"
    """
}