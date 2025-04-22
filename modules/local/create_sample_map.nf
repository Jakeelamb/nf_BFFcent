process CREATE_SAMPLE_MAP {
    tag "Interval $interval"
    publishDir "${params.outdir}/sample_maps", mode: params.publish_dir_mode
    cpus 1
    memory '2 GB'
    time '1h'

    input:
    tuple val(interval), val(sample_gvcf_tuples)

    output:
    tuple val(interval), path("${interval}_sample_map.txt"), emit: sample_map

    script:
    """
    # Create sample map for GenomicsDBImport
    touch ${interval}_sample_map.txt
    
    # Add each sample and its GVCF to the map
    ${sample_gvcf_tuples.collect { sample_id, gvcf_path -> "echo \"${sample_id}\t${gvcf_path}\" >> ${interval}_sample_map.txt" }.join('\n')}
    
    echo "Created sample map for interval ${interval} with ${sample_gvcf_tuples.size()} samples"
    """
} 