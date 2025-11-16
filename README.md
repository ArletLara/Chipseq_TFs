# Chipseq_TFs
Modular ChIP-seq workflow: FASTQ preprocessing, peak calling &amp; motifs, peak selection &amp; counting, and R downstream analysis

# Description
This repository contains an end-to-end ChIP-seq analysis workflow for TF ChIP-Seq data. It covers preprocessing of sequencing files, peak calling and motif discovery, candidate peak selection and counts datasets. Will be updated with  downstream analyses in R, which can be useful in a case-by-case basis. The workflow is designed to run on SLURM-managed HPC systems and can be adapted to different projects.

## Module 1: Preprocessing of the sequencing files
- Per-sample job generation (`scripts/chipseq_pipeline.sh`)
- TrimGalore → FastQC → Bowtie2 → samtools → bedtools
- Duplicate removal and blacklist filtering for mm10

## Module 2: Peak calling and Motif discovery
This module takes the filtered, duplicate-removed BAM files from Module 1 and performs TF ChIP-seq peak calling and motif analysis. For each treated/input pair, it runs MACS2 narrow peak calling, filters high-confidence peaks (score > 70), computes FRiP scores for both treated and input, and uses HOMER (mm10) to discover enriched motifs in the filtered peaks. It also generates browser-ready tracks (BigBed peaks, raw and CPM-normalized bigWigs, and fold-change bigWigs), and summarizes QC metrics such as blacklist removal, mapping flag stats, chrM/chrUn composition, library scaling factors, and a multi-sample fingerprint plot.

## Module 3: Peak selection and Counts
This module takes the MACS2 peak sets and summit files generated in Module 2 and constructs a reproducible, nonredundant peak universe for the selected ChIP mark. It merges peak intervals across condition replicates, filters for peaks supported by at least two contributing intervals, and then combines the different condition peak sets to produce a unified candidate peak list. From this candidate set, it identifies a single representative summit per peak (highest MACS2 score).

Using these finalized regions, the module quantifies read enrichment across all samples by running multiBamSummary in three modes:
Candidate peaks – coverage over merged reproducible peak intervals
Unique summits – 1-bp apex positions representing each candidate peak
Non-peak genomic bins – fixed-size background bins outside candidate peaks for normalization/QC

## Module 4: Downstream analyses in R
