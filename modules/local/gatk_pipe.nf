process gatk_pipe {
    tag "joint_genotyping"
    publishDir "${params.outdir}/gatk_pipe", mode: params.publish_dir_mode
    cpus 32
    memory '64 GB'
    time '24h'

    input:
    path gvcfs
    tuple path(ref_genome), path(ref_genome_dict), path(ref_genome_fasta_index)

    output:
    path "filtered_snps.vcf.gz", emit: filtered_vcf
    path "filtered_snps.vcf.gz.tbi", emit: filtered_vcf_index

    script:
    // Create GVCF parameters string
    def gvcf_params = gvcfs.collect { "-V $it" }.join(' ')
    """
    # Create temporary directory
    mkdir -p temp
    mkdir -p gendb

    # Step 1: Import GVCFs to GenomicsDB
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=temp" GenomicsDBImport \\
        ${gvcf_params} \\
        --genomicsdb-workspace-path "gendb" \\
        --batch-size 50 \\
        --interval-padding 100 \\
        --intervals "${workflow.projectDir}/intervals_chromo_only_nopos.list" \\
        --genomicsdb-shared-posixfs-optimizations \\
        --reader-threads 2

    # Step 2: GenotypeGVCFs
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=temp" GenotypeGVCFs \\
        -R "${ref_genome}" \\
        -V "gendb://gendb" \\
        -O "raw_variants.vcf.gz" \\
        --reader-threads 2

    # Step 3: Select and filter SNPs
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=temp" SelectVariants \\
        -R "${ref_genome}" \\
        -V "raw_variants.vcf.gz" \\
        --select-type-to-include SNP \\
        -O "raw_snps.vcf.gz"

    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=temp" VariantFiltration \\
        -R "${ref_genome}" \\
        -V "raw_snps.vcf.gz" \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "FS > 60.0" --filter-name "FS60" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -filter "AF < 0.05" --filter-name "minMAF" \\
        -filter "DP < 5 || DP > 100" --filter-name "Depth" \\
        -O "temp_filtered_snps.vcf.gz"

    gatk SelectVariants \\
        -V "temp_filtered_snps.vcf.gz" \\
        --select-type-to-include SNP \\
        --restrict-alleles-to BIALLELIC \\
        --exclude-filtered true \\
        -O "filtered_snps.vcf.gz"
    """
}