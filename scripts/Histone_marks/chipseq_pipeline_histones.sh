#!/bin/bash -ex
# ─────────────────────────────────────────────────────────────────────────────
#  Script : chipseq_pipeline.sh
#  Author : Arlet Lara
#  Updated: 2025-11-03
#  Purpose: Automate the processing of Illumina H3K27ac ChIP-seq paired-end data.
#
#  Overview:
#    Each job runs the end-to-end ChIP-seq workflow on that sample.
#    This script submits the jobs so you just have to run with bash script.sh
#
#  Steps each per-sample job performs:
#    1) Adapter trimming with TrimGalore
#    2) Quality control with FastQC (on trimmed reads)
#    2.1) Plasmid contamination screen (external helper script)
#    3) Alignment to the mouse genome (mm10) with Bowtie2
#    4) Paired-end BAM preparation for MACS2 (mark duplicates; name-sort → fixmate → coord-sort)
#    5) Blacklist filtering with bedtools (remove mm10 ENCODE blacklist intervals) + Filtering of chrUN chrM regions
#
#  Inputs:
#    raw_data/${sample}_R{1,2}.fastq.gz
#    (Assumes `${sample}` matches the basename before `_R1.fastq.gz` / `_R2.fastq.gz`.)
#
#  Outputs:
#    - Per-sample FastQC reports (trimmed reads) in download/
#    - Coordinate-sorted + indexed BAM after blacklist filtering in filtered/
#    - BigWig coverage tracks (CPM-normalized) in bigwigs/
#
#  Requires:
#    - TrimGalore ≥ 0.6.10
#    - FastQC ≥ 0.11.9
#    - Bowtie2 ≥ 2.5
#    - samtools ≥ 1.10
#    - bedtools ≥ 2.29
#    - deepTools ≥ 3.5
#    - SLURM, mm10 Bowtie2 index, mm10 blacklist BED
#
#  Notes:
#    - Duplicate reads are not removed
#    - The heredoc uses `>` (overwrite) when writing job files; re-running this controller
#      will append to existing `.qsh` files—delete them first to avoid duplication.
#    - `--cpus-per-task=8` aligns with Bowtie2 `-p 8` and samtools `-@ 8`.
#    - The plasmid contamination helper is invoked on filenames ending `_R1.fastq.gz`
#      and `_R2.fastq.gz`. Adjust its `-1/-2` suffixes if your raw reads follow the
#      `_R1.fastq.gz` / `_R2.fastq.gz` convention or live outside the working dir.
#    - Designed for SLURM-managed HPC environments; paths assume shared filesystems.
# ─────────────────────────────────────────────────────────────────────────────


# ── User configuration section ───────────────────────────────────────────────
# Root path for this project 
path_to_anaconda=/path/to/anaconda3
path_to_my_profile=/path/to/my_profile
project=project_name

# Path to source directory containing input data and project folders
# Assumes raw fastqs exist under $base/raw_data/
base=$path_to_my_profile/data/$project

# List of samples (update with your experiment names)
# These basenames must match raw_data/${sample}_R{1,2}.fastq.gz
# List of samples (update with your experiment names)
Files=('sample1_chip' 'sample2_chip' 'sample1_IgG' 'samplsample2_IgGe2')


# Create a working directory for temporary outputs
# why: Centralizes SLURM stdout/err and large intermediates on scratch storage.
mkdir -p /path/to/scratch/$project
cd $base

# Subdirectories for pipeline stages (created once here; jobs write into them)
mkdir -p scripts trimmed aligned filtered bigwigs plasmid_contamination download/fastqc marked_dup 
# Note: Ensure the Bowtie2 index and blacklist paths referenced below are accessible.

# Step 1 ▸ Generate one SLURM job script per sample
# why: Produces reproducible, inspectable job files; submission is a separate step.
for sample in "${Files[@]}"; do
qsh=scripts/${sample}_ChIPseq.qsh

cat <<EOF> $qsh
#!/bin/bash

#SBATCH --job-name=ChIPseq_$sample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --time=20:00:00
#SBATCH --output=/path/to/scratch/$project/ChIPseq_$sample.out
#SBATCH --error=/path/to/scratch/$project/ChIPseq_$sample.err
#SBATCH --mail-type=END

# Initialize conda for this non-interactive shell and activate the environment
# why: Ensures required bioinformatics tools and versions are available in PATH.
env_path=$path_to_anaconda/envs/chipseq_env
export PATH="\$env_path/bin:\$PATH"

# Move to trimming directory
cd $base/trimmed

# Step 1 ▸ Adapter and quality trimming
# why: Removes adapters/low-quality bases; outputs *_val_{1,2}.fq.gz used downstream.
$path_to_my_profile/apps/TrimGalore-0.6.10/trim_galore --paired \
../raw_data/${sample}_R1_001.fastq.gz ../raw_data/${sample}_R2_001.fastq.gz

# Step 2 ▸ Quality control on trimmed reads
# why: Verify trimming outcomes and per-base quality before alignment.
fastqc ${sample}_R1_001_val_1.fq.gz -o ../download/fastqc/
fastqc ${sample}_R2_001_val_2.fq.gz -o ../download/fastqc/

# Step 2.1 ▸  Plasmid contamination
cd ..
# why: External helper runs from project root; expects raw FASTQ naming.
# Adjust -1/-2 suffixes in the call below if your raw files
$path_to_my_profile/data/Common_files/plasmid_contamination_files/plasmid_contamination.sh \
-s ${sample} -1 "_R1_001_val_1.fq.gz" -2 "_R2_001_val_2.fq.gz"

# Step 3 ▸ Alignment to mouse genome (mm10)
# why: Bowtie2 alignment using pre-built mm10 index; threads match SLURM cpus.
$path_to_my_profile/apps/bowtie2-2.5.3-linux-x86_64/bowtie2 -p 8 \
-x $path_to_my_profile/data/Common_files/mm10_index_bt2/mm10 \
-1 trimmed/${sample}_R1_001_val_1.fq.gz -2 trimmed/${sample}_R2_001_val_2.fq.gz \
-S aligned/${sample}.sam

# Convert SAM to BAM and cleanup
# why: BAM is compressed and indexable; SAM removed to save space.
samtools view -@ 8 -bS aligned/${sample}.sam > aligned/${sample}.bam
rm aligned/${sample}.sam


# Step 4 ▸ Prepare paired-end BAM for MACS2 (retain duplicates)
### fixmate requires name-sorted input; then coordinate-sort for indexing/IO efficiency.
samtools sort -@ 8 -n -o aligned/${sample}_namesorted.bam aligned/${sample}.bam
samtools fixmate -m aligned/${sample}_namesorted.bam aligned/${sample}_fixmate.bam
samtools sort -@ 8 -o aligned/${sample}_coordsorted.bam aligned/${sample}_fixmate.bam
samtools index aligned/${sample}_coordsorted.bam

### Mark duplicates
samtools markdup aligned/${sample}_coordsorted.bam marked_dup/${sample}.bam


# Step 5 ▸ Remove blacklisted regions
# why: Exclude systematic artifacts (mm10 blacklist) before track generation/peak calling.
bedtools intersect -v -abam marked_dup/${sample}.bam \
-b $path_to_my_profile/data/Common_files/blacklist_mm10.bed \
> marked_dup/${sample}_blacklisted_with_chrMUncl.bam

bedtools intersect -v -abam marked_dup/${sample}_blacklisted_with_chrMUncl.bam \
-b $path_to_my_profile/data/Common_files/exclude_chrM__UnclChr_mm10.bed \
> filtered/${sample}_blacklisted.bam

# Index 
samtools index filtered/${sample}_blacklisted.bam


EOF
done

for file in scripts/*_ChIPseq.qsh; do sbatch $file;done
