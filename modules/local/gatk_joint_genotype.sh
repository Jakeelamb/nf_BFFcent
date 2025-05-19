#!/bin/bash

#SBATCH --job-name=gatk_joint_genotype
#SBATCH --time=156:00:00
#SBATCH --partition=amem
#SBATCH --cpus-per-task=32
#SBATCH --qos=mem
#SBATCH --output=gatk_joint_genotype_%j.out
#SBATCH --error=gatk_joint_genotype_%j.err
#SBATCH --mem=120G

source ~/.bashrc
conda activate maize_vc

# Check if required arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <interval> <sample_map> <reference_genome>"
    echo "  interval: Genomic interval to process (e.g. chr1:1-1000000)"
    echo "  sample_map: Path to the sample map file"
    echo "  reference_genome: Path to the reference genome FASTA file"
    exit 1
fi

# Get command line arguments
interval=$1
sample_map=$2
reference_genome=$3

# Set up variables
mem=120
cpus=32
outdir="gatk_pipe"

# Create output directory
mkdir -p "$outdir"

# Define reference path variables
REF_PATH="$reference_genome"
REF_DIR=$(dirname "$REF_PATH")
REF_NAME=$(basename "$REF_PATH")
REF_BASE=$(echo "$REF_NAME" | sed 's/\.[^.]*$//')

echo "Processing interval: $interval"
echo "Sample map: $sample_map"
echo "Reference genome: $reference_genome"

# List current directory for debugging
echo "Contents of working directory before symlinks:"
ls -la > workdir_contents_before.txt

# Create symlink for reference
ln -sf "$reference_genome" "genome.fa"

# Create symlink for .fai index
ln -sf "${REF_PATH}.fai" "genome.fa.fai"
if [ ! -e "genome.fa.fai" ]; then
    echo "Error: genome.fa.fai not found. Check if ${REF_PATH}.fai exists."
    exit 1
fi

# Create symlink for .dict file
ln -sf "${REF_DIR}/${REF_BASE}.dict" "genome.dict"
if [ ! -e "genome.dict" ]; then
    echo "Error: genome.dict not found. Check if ${REF_DIR}/${REF_BASE}.dict exists."
    exit 1
fi

# List current directory after linking for debugging
echo "Contents of working directory after symlinks:"
ls -la > workdir_contents_after.txt

# Create a unique workspace for this interval
WORKSPACE="genomicsdb_${interval//[:\/-]/_}"

# Create tmp directory
mkdir -p ./tmp

# Import GVCFs to GenomicsDB
echo "Starting GenomicsDBImport for interval ${interval}..."
gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" GenomicsDBImport \
    --genomicsdb-workspace-path $WORKSPACE \
    --batch-size 50 \
    -L ${interval} \
    --sample-name-map ${sample_map} \
    --reader-threads $(( cpus / 4 )) \
    --max-num-intervals-to-import-in-parallel $(( cpus / 8 ))

# Run joint genotyping on the interval
echo "Starting GenotypeGVCFs for interval ${interval}..."
gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" GenotypeGVCFs \
    -R genome.fa \
    -V gendb://$WORKSPACE \
    -O ${interval//[:\/-]/_}.vcf.gz \
    -L ${interval} \
    --tmp-dir=./tmp

# Apply variant filtration
echo "Starting VariantFiltration for interval ${interval}..."
gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R genome.fa \
    -V ${interval//[:\/-]/_}.vcf.gz \
    -O ${outdir}/${interval//[:\/-]/_}.filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "basic_snp_filter"
    
echo "Completed processing for interval ${interval}"
echo "Output files are in ${outdir}/${interval//[:\/-]/_}.filtered.vcf.gz"


