#!/usr/bin/env nextflow

// Import modules
include { REFERENCE_INDEX } from './modules/local/reference_index'
include { FASTP_TRIM } from './modules/local/fastp_trim'
include { BWAMEM2_ALIGN } from './modules/local/bwamem2_align'
include { SAMTOOLS_PROCESS } from './modules/local/samtools_process'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk_haplotypecaller'
include { GATK_JOINT_GENOTYPE } from './modules/local/gatk_joint_genotype'
include { PLINK2_ANALYSIS } from './modules/local/plink2_analysis'

// Define a process to collect readgroup info
process READGROUP_COLLECT {
    publishDir "${params.outdir}/readgroups", mode: params.publish_dir_mode
    
    input:
    path(readgroup_files)
    
    output:
    path "combined_readgroup_info.txt"
    
    script:
    """
    # Get the header from the first file
    head -n 2 ${readgroup_files[0]} > combined_readgroup_info.txt
    
    # Concatenate the data from all files, skipping headers
    for file in ${readgroup_files}; do
        if [ "\$file" != "${readgroup_files[0]}" ]; then
            tail -n +3 \$file >> combined_readgroup_info.txt
        else
            tail -n +3 \$file >> combined_readgroup_info.txt
        fi
    done
    
    # Add summary information
    echo -e "\\nSummary:" >> combined_readgroup_info.txt
    echo "--------" >> combined_readgroup_info.txt
    total_samples=\$(grep -v "Sample\\|===" combined_readgroup_info.txt | wc -l)
    echo "Total samples processed: \$total_samples" >> combined_readgroup_info.txt
    
    # Perform validation
    echo -e "\\nValidation Checks:" >> combined_readgroup_info.txt
    echo "----------------" >> combined_readgroup_info.txt
    if grep -P "\\t\\t" combined_readgroup_info.txt > /dev/null; then
        echo "WARNING: Some entries have missing values. Please review the output." >> combined_readgroup_info.txt
    else
        echo "All entries have complete information." >> combined_readgroup_info.txt
    fi
    """
}

// Define workflow
workflow {
    // Parse input samples
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, [file(row.read1), file(row.read2)]) }
        .set { samples_ch }
    
    // Reference genome processing
    ref_genome = file(params.reference_genome)
    ref_genome_dict = file("${params.reference_genome.replaceFirst(/\.fa(sta)?$/, '.dict')}")
    ref_genome_fasta_index = file("${params.reference_genome}.fai")
    
    // Run reference indexing if needed
    REFERENCE_INDEX(tuple(ref_genome, ref_genome_dict, ref_genome_fasta_index))
    
    // Trimming can start immediately
    FASTP_TRIM(samples_ch)
    
    // Alignment needs indexed reference genome
    // Combine the trimmed reads with indexed reference
    FASTP_TRIM.out.trimmed_reads
        .combine(REFERENCE_INDEX.out.indexed_ref.first())
        .map { it -> [it[0], it[1], it[2]] }
        .set { reads_and_ref }
    
    // Run alignment with reference dependency
    BWAMEM2_ALIGN(reads_and_ref)
    
    // Collect and combine all readgroup info files
    BWAMEM2_ALIGN.out.readgroup_info
        .collect()
        .set { all_readgroup_info }
    
    // Generate combined readgroup report
    READGROUP_COLLECT(all_readgroup_info)
    
    // SAM/BAM processing
    SAMTOOLS_PROCESS(BWAMEM2_ALIGN.out.bam_files)
    
    // Wait for all samples to complete processing before running variant calling
    // Combine with reference genome for GATK
    SAMTOOLS_PROCESS.out.processed_bam
        .combine(REFERENCE_INDEX.out.indexed_ref.first())
        .map { sample_id, bam, ref, dict, fai -> [sample_id, bam, ref] }
        .set { all_processed_bams_with_ref }
    
    // Variant calling with reference genome
    GATK_HAPLOTYPECALLER(all_processed_bams_with_ref)
    
    // Collect all GVCF files for joint genotyping
    GATK_HAPLOTYPECALLER.out.gvcfs
        .map { sample_id, gvcf, idx -> gvcf }
        .collect()
        .set { all_gvcfs }
    
    // Get reference for joint genotyping
    REFERENCE_INDEX.out.indexed_ref
        .first()
        .map { ref, dict, fai -> [ref, dict, fai] }
        .set { indexed_ref_for_gatk }
    
    // Joint genotyping and filtering
    GATK_JOINT_GENOTYPE(all_gvcfs, indexed_ref_for_gatk)
    
    // Population genetics analysis
    PLINK2_ANALYSIS(GATK_JOINT_GENOTYPE.out.filtered_vcf)
}

// Log workflow completion
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}
