process plink2 {
    tag "population_analysis"
    publishDir "${params.outdir}/plink2", mode: params.publish_dir_mode
    cpus 16
    memory '64 GB'
    time '24h'

    input:
    path vcf_file

    output:
    path "pca_results.*", emit: pca_results
    path "binary_data.*", emit: binary_data
    path "allele_freq.afreq", emit: allele_freq

    script:
    """
    # Step 1: LD pruning
    plink2 \\
        --vcf "${vcf_file}" \\
        --double-id \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.1 \\
        --out pruned_snps

    # Step 2: Create binary files
    plink2 \\
        --vcf "${vcf_file}" \\
        --double-id \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --extract pruned_snps.prune.in \\
        --make-pgen \\
        --out binary_data

    # Step 3: Calculate frequencies
    plink2 \\
        --pfile binary_data \\
        --freq \\
        --allow-extra-chr \\
        --out allele_freq

    # Step 4: Run PCA
    plink2 \\
        --pfile binary_data \\
        --read-freq allele_freq.afreq \\
        --pca \\
        --allow-extra-chr \\
        --out pca_results
    """
}