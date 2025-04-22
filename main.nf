#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define Parameters
params.reads = "samples.txt" // Input samplesheet (sample,fastq_1,fastq_2)
params.reference_genome = "$projectDir/reference/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa" // Reference genome in reference directory
params.intervals_file = "intervals_chromo_only_nopos.list" // File listing intervals (1-10)
params.outdir = "./results"
params.publish_dir_mode = 'copy' // Or 'link', 'rellink', etc.
params.skip_index_check = false // Option to skip index checking (speeds up initialization)

// Define Log Level
log.info """\
         V A R I A N T - C A L L I N G   P I P E L I N E
         ================================================
         reads              : ${params.reads}
         reference_genome   : ${params.reference_genome}
         intervals          : ${params.intervals_file}
         outdir             : ${params.outdir}
         skip_index_check   : ${params.skip_index_check}
         ------------------------------------------------
         """
         .stripIndent()

// Include Modules
include { EXTRACT_FASTQ_METADATA } from './modules/local/extract_fastq_metadata'
include { FASTP_TRIM } from './modules/local/fastp_trim'
include { BWAMEM2_ALIGN } from './modules/local/bwamem2_align'
include { SAMTOOLS_PROCESS } from './modules/local/samtools_process'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk_haplotypecaller'
include { REFERENCE_INDEX } from './modules/local/reference_index'
include { CREATE_SAMPLE_MAP } from './modules/local/create_sample_map'
include { GATK_JOINT_GENOTYPE } from './modules/local/gatk_joint_genotype'
include { PLINK2_ANALYSIS } from './modules/local/plink2_analysis'
include { CONCAT_GVCFS } from './modules/local/concat_gvcfs'

// Workflow Definition
workflow {

    // --- Input Channels ---
    ch_reads = Channel.fromPath(params.reads)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [id: row.sample_id]
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            
            [ meta, [read1, read2] ]
        }
    
    ch_intervals = Channel.fromPath(params.intervals_file).splitText().map { it.trim() }.filter { it } // Emits interval names
    
    // Reference genome file
    ch_ref_input = Channel.value(file(params.reference_genome, checkIfExists: true))

    // --- Step 0: Index Reference Genome (if needed) ---
    // Create paths for expected index/dict files based on reference genome path
    def ref_file = file(params.reference_genome)
    def ref_dir = ref_file.getParent()
    def ref_name = ref_file.getName()
    def ref_base_name = ref_name.replaceFirst(/\.gz$/, '').replaceFirst(/\.fasta$/, '').replaceFirst(/\.fa$/, '') // Handle .fa, .fasta, .gz extensions

    def ref_dict_path = file("${ref_dir}/${ref_base_name}.dict")
    def ref_fai_path = file("${ref_dir}/${ref_name}.fai") // FAI should match original extension
    
    // Get uncompressed name for BWA indexes (index files won't have .gz even if the reference is compressed)
    def ref_for_bwa_indexes = params.reference_genome.replaceFirst(/\.gz$/, '')

    // Flag to check if indexing is needed
    def indexing_needed = true

    // Skip index checking if requested (useful for repeated runs where indexes are known to exist)
    if (!params.skip_index_check) {
        // Check for required index files
        def required_extensions = ['.amb', '.ann', '.bwt.2bit.64', '.pac', '.sa', '.0123']
        indexing_needed = required_extensions.any { ext ->
            !file("${ref_for_bwa_indexes}${ext}").exists()
        } || !ref_dict_path.exists() || !ref_fai_path.exists()
        
        // Log the index paths we're checking
        log.debug "Checking for reference index files:"
        required_extensions.each { ext -> log.debug " - ${ref_for_bwa_indexes}${ext}" }
        log.debug " - ${ref_dict_path}"
        log.debug " - ${ref_fai_path}"
    }

    // Run indexing only if needed
    if (indexing_needed) {
        log.info "Reference genome index files not found or check requested. Will run indexing."
        REFERENCE_INDEX(
            ch_ref_input
                .map { ref -> tuple(ref, ref_dict_path, ref_fai_path) } // Pass ref path and target index/dict paths
        )
        ch_indexed_reference = REFERENCE_INDEX.out.indexed_ref // [ref_genome, ref_genome_dict, ref_genome_fasta_index]
    } else {
        log.info "All reference genome index files found. Skipping indexing step."
        // Create the reference channel directly if indexes exist
        ch_indexed_reference = Channel.value(tuple(ref_file, ref_dict_path, ref_fai_path))
    }

    // --- Step 1: Extract FASTQ Metadata ---
    EXTRACT_FASTQ_METADATA(ch_reads)

    // --- Step 2: Read Trimming ---
    FASTP_TRIM(ch_reads)

    // --- Step 3: Alignment ---
    // BWAMEM2_ALIGN expects [meta, reads, ref, rg_info]
    // ch_indexed_reference provides [ref_genome, ref_genome_dict, ref_genome_fasta_index]
    ch_align_input = FASTP_TRIM.out.trimmed_reads
        .join(EXTRACT_FASTQ_METADATA.out.rg_info) // [meta, reads, rg_info]
        .combine(ch_indexed_reference.map{ ref, dict, fai -> ref }) // Combine with just the ref path from indexed_ref
        .map { meta, reads, rg_info, ref -> [meta, reads, ref, rg_info] } // Reorder for BWAMEM2_ALIGN

    BWAMEM2_ALIGN(ch_align_input)

    // --- Step 4: SAMtools Processing ---
    SAMTOOLS_PROCESS(BWAMEM2_ALIGN.out.bam_files)

    // --- Step 5: Haplotype Calling (per sample, per interval) ---
    // GATK_HAPLOTYPECALLER expects [meta, bam, interval, ref]
    ch_haplotype_input = SAMTOOLS_PROCESS.out.processed_bam
                             .combine(ch_intervals)   // [meta, bam, interval]
                             .combine(ch_indexed_reference.map{ ref, dict, fai -> ref })   // Combine with just the ref path
    // Map to ensure correct order: [meta, bam_file, interval, reference]
    ch_haplotype_input = ch_haplotype_input.map { meta, bam, interval, ref -> tuple(meta, bam, interval, ref) }

    GATK_HAPLOTYPECALLER(ch_haplotype_input)

    // --- Step 6: Prepare for Joint Genotyping (Create Sample Map per Interval) ---
    ch_gvcf_info = GATK_HAPLOTYPECALLER.out.gvcfs_interval
        .map { meta, interval, gvcf_path, tbi_path ->
            // Check if the TBI file actually exists
            if (file(tbi_path).exists()) {
                return tuple(interval, tuple(meta.id, gvcf_path)) // Format for grouping: [interval, [sample_id, gvcf_path]]
            } else {
                log.warn "Missing index file ${tbi_path} for ${gvcf_path}. Skipping this GVCF for interval ${interval}."
                return null // Indicate failure
            }
        }
        .filter { it != null } // Remove entries where index was missing

    // Group by interval -> [interval, [[sample1, path1], [sample2, path2], ...]]
    ch_grouped_gvcf_info = ch_gvcf_info.groupTuple()

    // --- Step 7: Concatenate GVCFs per sample (across intervals) ---
    // Group GVCFs by sample ID instead of by interval
    ch_sample_gvcfs = GATK_HAPLOTYPECALLER.out.gvcfs_interval
        .map { meta, interval, gvcf_path, tbi_path ->
            return tuple(meta.id, tuple(interval, gvcf_path, tbi_path))
        }
        .groupTuple() // Group by sample ID: [sample_id, [[interval1, gvcf1, tbi1], [interval2, gvcf2, tbi2], ...]]

    // Concatenate GVCFs for each sample
    CONCAT_GVCFS(ch_sample_gvcfs)

    // --- Step 8: Create sample map for joint genotyping ---
    CREATE_SAMPLE_MAP(ch_grouped_gvcf_info) // Emits: [interval, map_file_path]

    // --- Step 9: Joint Genotyping (per interval) ---
    ch_joint_genotype_input = CREATE_SAMPLE_MAP.out.sample_map
                                 .combine(ch_indexed_reference.map{ ref, dict, fai -> ref })  // [interval, map_file, ref_file]

    GATK_JOINT_GENOTYPE(ch_joint_genotype_input)

    // --- Step 10: PLINK2 Analysis ---
    // Collect all interval VCFs to analyze with PLINK2
    // For simplicity, use the first interval's VCF for demonstration
    // In a real scenario, you might want to merge all interval VCFs first
    PLINK2_ANALYSIS(GATK_JOINT_GENOTYPE.out.filtered_vcf_interval.map { interval, vcf, tbi -> vcf }.first())

    // --- Workflow Output ---
    ch_final_vcfs = GATK_JOINT_GENOTYPE.out.filtered_vcf_interval
    ch_final_vcfs.view { interval, vcf, tbi -> "Final VCF for interval ${interval}: ${vcf}" }
}

// Workflow Completion Message
workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed successfully at: ${workflow.complete}
        Duration               : ${workflow.duration}
        CPU hours              : ${workflow.cpuHours}
        Memory Mb-hours        : ${workflow.memoryMbHours}
        Work directory         : ${workflow.workDir}
        Output directory       : ${params.outdir}
        ================================================
        Pipeline execution finished!
        """ : """
        Pipeline execution failed!
        ---------------------------
        Error summary: ${workflow.errorReport}
        Work directory: ${workflow.workDir}
        Check the log file for details: .nextflow.log
        ================================================
        """
    )
}