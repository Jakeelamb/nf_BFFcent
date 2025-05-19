process CONCAT_GVCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/concatenated_gvcfs", mode: params.publish_dir_mode
    cpus 2
    memory '16 GB'
    time '2h'

    input:
    tuple val(sample_id), val(interval_gvcf_tuples) // [sample_id, [[interval1, gvcf1, tbi1], [interval2, gvcf2, tbi2], ...]]

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: concatenated_gvcfs

    script:
    // Create a list of GVCF files to concatenate, sorted by interval
    def gvcf_files_str = interval_gvcf_tuples
        .sort { a, b -> a[0].toString() <=> b[0].toString() } // Sort by interval
        .collect { interval, gvcf, tbi -> "${gvcf}" }
        .join(" ")
    """
    # Create a file list for bcftools concat
    echo "${gvcf_files_str}" | tr ' ' '\\n' > ${sample_id}_file_list.txt
    
    # Verify the order of files
    echo "Files to be concatenated in order:"
    cat ${sample_id}_file_list.txt
    
    # Concatenate the files
    echo "Concatenating GVCFs for sample ${sample_id}..."
    bcftools concat -f ${sample_id}_file_list.txt -O z -o ${sample_id}.g.vcf.gz
    
    # Index the concatenated GVCF
    bcftools index -t ${sample_id}.g.vcf.gz
    
    echo "Concatenation complete for sample ${sample_id}"
    """
} 