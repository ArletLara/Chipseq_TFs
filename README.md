# Chipseq_TFs
Modular ChIP-seq workflow: FASTQ preprocessing, peak calling &amp; motifs, peak selection &amp; counting, and R downstream analysis

# Description
This repository contains an end-to-end ChIP-seq analysis workflow for TF ChIP-Seq data. It covers preprocessing of sequencing files, peak calling and motif discovery, candidate peak selection and counts datasets. Will be updated with  downstream analyses in R, which can be useful in a case-by-case basis. The workflow is designed to run on SLURM-managed HPC systems and can be adapted to different projects.

## Module 1: Preprocessing of the sequencing files
- Per-sample job generation (`scripts/chipseq_pipeline.sh`)
- TrimGalore → FastQC → Bowtie2 → samtools → bedtools
- Duplicate removal and blacklist filtering for mm10

## Module 2: Peak calling and Motif discovery


## Module 3: Peak selection and Counts


## Module 4: Downstream analyses in R
