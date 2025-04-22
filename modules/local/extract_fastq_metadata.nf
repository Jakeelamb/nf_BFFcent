process EXTRACT_FASTQ_METADATA {
    tag "$meta.id"
    cpus 1
    memory '2 GB'
    time '1h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_rg_info.txt"), emit: rg_info

    script:
    """
    # Extract the first header line from read1
    zcat "${reads[0]}" | head -n 1 > header.txt

    # Parse Illumina header format
    # Expected: @A01180:202:H5JHKDSX7:3:1101:7563:1000 1:N:0:CTGTGGTGAC+CCTGTCTGTC
    HEADER=\$(cat header.txt)
    
    # Extract fields
    INSTRUMENT=\$(echo "\$HEADER" | awk -F':' '{print \$1}' | sed 's/@//')
    FLOWCELL=\$(echo "\$HEADER" | awk -F':' '{print \$3}')
    LANE=\$(echo "\$HEADER" | awk -F':' '{print \$4}')
    BARCODE=\$(echo "\$HEADER" | awk '{print \$2}' | awk -F':' '{print \$4}' | sed 's/+//g')
    
    # Set defaults if fields are missing
    [ -z "\$INSTRUMENT" ] && INSTRUMENT="unknown"
    [ -z "\$FLOWCELL" ] && FLOWCELL="unknown"
    [ -z "\$LANE" ] && LANE="1"
    [ -z "\$BARCODE" ] && BARCODE="none"
    
    # Construct read group fields
    RG_ID="${meta.id}.\${FLOWCELL}.\${LANE}"
    RG_SM="${meta.id}"
    RG_LB="${meta.id}.lib"
    RG_PU="\${FLOWCELL}.\${LANE}.\${BARCODE}"
    RG_PL="ILLUMINA"
    
    # Write to output file
    echo -e "ID:\${RG_ID}\tSM:\${RG_SM}\tLB:\${RG_LB}\tPU:\${RG_PU}\tPL:\${RG_PL}" > ${meta.id}_rg_info.txt
    """
} 