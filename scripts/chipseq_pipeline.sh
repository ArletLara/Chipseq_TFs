#!/bin/bash -ex
# ─────────────────────────────────────────────────────────────────────────────
#  Script : chipseq_pipeline.sh
#  Author : Arlet Lara-Custodio
#  Updated: 2025-11-13
#  Purpose: Automate the processing of Illumina TF ChIP-seq paired-end data
#           on a SLURM-managed HPC cluster.
#
#  Overview:
#    This "controller" script:
#      - Generates one SLURM job script per sample.
#      - Submits all jobs via `sbatch`.
#    Each generated job script runs the ChIP-seq preprocessing workflow
#    for that sample.
#
#  Per-sample workflow:
#    1) Adapter trimming with TrimGalore
#    2) Quality control with FastQC (on trimmed reads)
#    2.1) Plasmid contamination screen (external helper script)
#    3) Alignment to the mouse genome (mm10) with Bowtie2
#    4) Paired-end BAM preparation for MACS2 (remove duplicates; name-sort → fixmate → coord-sort)
#    5) Blacklist filtering with bedtools (remove mm10 ENCODE blacklist intervals) + Filtering of chrUN chrM regions
#
#  Inputs:
#    - raw_data/${sample}_R1.fastq.gz
#    - raw_data/${sample}_R2.fastq.gz
#    (Assumes `${sample}` matches the basename before `_R1.fastq.gz` / `_R2.fastq.gz`.)
#
#  Outputs:
#    - Per-sample FastQC reports (trimmed reads) in download/fastqc/
#    - Duplicate-filtered, blacklist-filtered BAM in filtered/${sample}_blacklisted.bam
#
#  Requires (software):
#    - TrimGalore ≥ 0.6.10
#    - FastQC ≥ 0.11.9
#    - Bowtie2 ≥ 2.5
#    - samtools ≥ 1.10
#    - bedtools ≥ 2.29
#    - deepTools ≥ 3.5 (for downstream coverage tracks; not used directly here)
#    - SLURM
#
#  Requires (reference files):
#    - mm10 Bowtie2 index
#    - mm10 blacklist BED
#    - BED file listing chrM and unplaced/chrUn regions to exclude
#
#
#  Notes:
#    - Duplicate reads (SAM flag 1024) are removed before blacklist filtering.
#    - The heredoc uses `>` (overwrite) when writing job files; re-running this
#      controller will REPLACE existing `.qsh` files.
#    - `--cpus-per-task=8` aligns with Bowtie2 `-p 8` and samtools `-@ 8`.
#    - The plasmid contamination helper is invoked on filenames ending `_R1.fastq.gz`
#      and `_R2.fastq.gz`. Adjust its `-1/-2` suffixes if your raw reads follow the
#      `_R1.fastq.gz` / `_R2.fastq.gz` convention or are placed in another working dir.
#    - Designed for SLURM-managed HPC environments; paths assume shared filesystems.
# ─────────────────────────────────────────────────────────────────────────────

# ── User configuration section ───────────────────────────────────────────────
# Name project as specified by initial folder name
project=project_name

# Path to source directory containing input data and project folders
# Assumes raw fastqs exist under $base/raw_data/
base=/path/to/my_profile/data/$project

# List of samples (update with your experiment names)
# These basenames must match raw_data/${sample}_R{1,2}.fastq.gz exactly
Files=('sample1' 'sample2' 'sample3' 'sample4')

# Directory for temporary outputs (scratch)
mkdir -p /path/to/scratch/$project
cd $base

# Subdirectories (created once here; jobs write into them).
mkdir -p scripts trimmed aligned filtered bigwigs_nodup plasmid_contamination download/fastqc marked_dup removed_dup


# Step 1 ▸ Generate one SLURM job script per sample
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
source /path/to/my_profile/anaconda3/etc/profile.d/conda.sh
conda activate chipseq_env

# Move to trimming directory
cd $base/trimmed

# Step 1 ▸ Adapter and quality trimming
# why: Removes adapters/low-quality bases; outputs *_val_{1,2}.fq.gz used downstream.
trim_galore --paired \
../raw_data/${sample}_R1.fastq.gz ../raw_data/${sample}_R2.fastq.gz


# Step 2 ▸ Quality control on trimmed reads
fastqc ${sample}_R1_val_1.fq.gz -o ../download/fastqc/
fastqc ${sample}_R2_val_2.fq.gz -o ../download/fastqc/

# Step 2.1 ▸ Plasmid contamination check on trimmed reads; external helper
cd ..
# Adjust -1/-2 suffixes in the call below if your trimmed file names differ.
/path/to/my_profile/data/Common_files/plasmid_contamination_files/plasmid_contamination.sh \
-s ${sample} -1 "_R1_val_1.fq.gz" -2 "_R2_val_2.fq.gz"

# Step 3 ▸ Alignment to mouse genome (mm10)
# Bowtie2 alignment using pre-built mm10 index.
bowtie2 -p 8 \
-x /path/to/my_profile/data/Common_files/mm10_index_bt2/mm10 \
-1 trimmed/${sample}_R1_val_1.fq.gz -2 trimmed/${sample}_R2_val_2.fq.gz \
-S aligned/${sample}.sam

# Convert SAM to BAM and cleanup
# BAM is compressed and indexable; SAM is removed to save space.
samtools view -@ 8 -bS aligned/${sample}.sam > aligned/${sample}.bam
rm aligned/${sample}.sam


# Step 4 ▸ Prepare paired-end BAM for MACS2 (remove duplicates)
# why: fixmate requires name-sorted input; then coordinate-sort for indexing/IO efficiency.
samtools sort -@ 8 -n -o aligned/${sample}_namesorted.bam aligned/${sample}.bam

### Fixmate
samtools fixmate -m aligned/${sample}_namesorted.bam aligned/${sample}_fixmate.bam

### Coordinate-sort for markdup
samtools sort -o aligned/${sample}_coordsorted.bam aligned/${sample}_fixmate.bam

### Mark (and remove) duplicates
samtools markdup aligned/${sample}_coordsorted.bam marked_dup/${sample}.bam
samtools view -b -F 1024 marked_dup/${sample}.bam > removed_dup/${sample}.bam

### Index 
samtools index removed_dup/${sample}.bam


# Step 5 ▸ Remove blacklisted regions
# why: Exclude systematic artifacts (ENCODE mm10 blacklist) and regions that are not of 
# interest before downstream visualization/analyses
bedtools intersect -v -abam removed_dup/${sample}.bam \
-b /mnt/BioAdHoc/Users/alara/data/Common_files/blacklist_mm10.bed \
> removed_dup/${sample}_blacklisted_with_chrMUncl.bam

bedtools intersect -v -abam removed_dup/${sample}_blacklisted_with_chrMUncl.bam \
-b /path/to/my_profile/data/Common_files/exclude_chrM__UnclChr_mm10.bed \
> filtered/${sample}_blacklisted.bam

# Index 
samtools index filtered/${sample}_blacklisted.bam

EOF
done

# Submit all generated jobs
for file in scripts/*_ChIPseq.qsh; do sbatch $file;done
