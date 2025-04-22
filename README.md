# nf_BFFcent: Black-Footed Ferret Variant Calling Pipeline

A Nextflow pipeline for variant discovery and population genomics analysis of Black-Footed Ferret samples.

## Overview

This pipeline performs the following steps:
1. Reference Genome Indexing - creates indexes if not already present
2. FASTQ Metadata Extraction - extract read group information
3. Quality Control and Trimming - FASTP removes adapters and low-quality bases
4. Read Alignment - BWA-MEM2 aligns reads to the reference genome
5. SAM/BAM Processing - Samtools handles sorting, fixmate, and duplicate marking
6. Variant Calling - GATK HaplotypeCaller generates GVCFs per sample, per interval
7. GVCF Concatenation - bcftools combines per-interval GVCFs into per-sample GVCFs
8. Sample Map Creation - prepares sample maps required for joint genotyping
9. Joint Genotyping - GATK GenotypeGVCFs combines and genotypes per-interval samples
10. Population Analysis - PLINK2 performs principal component analysis

## Requirements

- Nextflow 21.10.0 or later
- Conda or Mamba (recommended for faster environment creation)
- SLURM workload manager (optional, for HPC execution)

## Usage

### 1. Installation

Clone this repository:
```bash
git clone https://github.com/yourusername/nf_BFFcent.git
cd nf_BFFcent
```

### 2. Conda Environment

The pipeline will automatically create and manage the required Conda environment using the provided `environment.yml` file. The environment will be created during the first run and reused in subsequent runs, so you don't need to manually create it.

If you want to make changes to the environment, simply edit the `environment.yml` file before running the pipeline.

### 3. Sample file preparation

Prepare your sample file (`samples.txt`) with tab-separated values:
```
sample_id    read1    read2
sample1    /path/to/sample1_R1.fastq.gz    /path/to/sample1_R2.fastq.gz
sample2    /path/to/sample2_R1.fastq.gz    /path/to/sample2_R2.fastq.gz
```

### 4. Reference genome

The pipeline is configured to use the Black-footed ferret reference genome. The path is set in the `nextflow.config` file. Before running, make sure to update this path to point to your copy of the reference genome:

```groovy
// In nextflow.config
params {
    reference_genome = "/path/to/blackfooted_ferret_genome.fasta"  // Update this path
    // ...
}
```

If needed, you can temporarily override this setting at runtime:
```bash
nextflow run main.nf --reference_genome /alternate/path/to/genome.fasta
```

The pipeline will automatically index the reference genome if the index files don't already exist. You can skip the index check with the `--skip_index_check` parameter if you know the indexes already exist and want a faster startup:

```bash
nextflow run main.nf --skip_index_check
```

### 5. Running the pipeline

Basic usage:
```bash
nextflow run main.nf \
  --input samples.txt \
  --outdir results
```

Using SLURM with resume capability:
```bash
nextflow run main.nf \
  -profile slurm \
  -resume \
  --outdir results
```

#### Additional Parameters

```bash
# With SLURM account
nextflow run main.nf \
  --input samples.txt \
  --outdir results \
  --slurm_account your_account \
  --slurm_partition partition_name \
  --slurm_qos qos_name

# With custom resources
nextflow run main.nf \
  --input samples.txt \
  --outdir results \
  --fastp_partition short-cpu \
  --alignment_partition long-highmem
```

## Pipeline Parameters

### Main Parameters
- `--input`: Path to samples.txt file (default: 'samples.txt')
- `--outdir`: Output directory (default: './results')
- `--temp_dir`: Directory for temporary files (default: './temp')
- `--skip_index_check`: Skip checking for reference genome indexes (default: false)

### SLURM Parameters
- `--slurm_account`: SLURM account (optional)
- `--slurm_partition`: Default SLURM partition (default: 'amilan')
- `--slurm_qos`: Default SLURM QOS (optional)

### Process-specific Parameters
Each process can have a custom partition and QOS defined:
- `--fastp_partition`, `--fastp_qos`
- `--alignment_partition`, `--alignment_qos`
- `--gatk_haplotype_partition`, `--gatk_haplotype_qos`
- `--concat_gvcfs_partition`, `--concat_gvcfs_qos`
- etc.

## Output Structure

```
results/
├── ref_index/              # Indexed reference genome files
├── trimming/               # Trimmed FASTQ files and QC reports
├── alignment/              # Aligned BAM files
├── readgroups/             # Read group information and summary
├── samtools/               # Processed BAM files
├── gatk_haplotype/         # Individual GVCF files per interval
├── concatenated_gvcfs/     # Concatenated GVCF files per sample
├── sample_maps/            # Interval-specific sample maps
├── gatk_pipe/              # Joint-called VCF files
├── plink2/                 # Population genetics outputs
└── pipeline_info/          # Execution reports and logs
```

## Read Group Extraction

The pipeline automatically extracts read group information from FASTQ headers to ensure proper sample tracking throughout the analysis. This information is used for the alignment step and is also saved in tabular format for reference.

The read group information includes:
- Sample ID
- Read Group ID (constructed from sample name, flowcell, and lane)
- Sample Name (SM tag)
- Library (LB tag)
- Platform Unit (PU tag, constructed from flowcell, lane, and barcode)
- Platform (PL tag, set to "Illumina")

## License

This pipeline is licensed under the MIT License.

## Contact

For questions or issues, please create a GitHub issue or contact the author.
