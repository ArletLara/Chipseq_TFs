#!/bin/bash -ex
#SBATCH --job-name=Peaks_QC_Motifs_fingerprint
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --output=/path/to/scratch/name_of_the_project/Peaks_QC_Motifs_fingerprint.out
#SBATCH --error=/path/to/scratch/name_of_the_project/Peaks_QC_Motifs_fingerprint.err
#SBATCH --mail-type=END

################################################################################
# Script : Peaks_QC_Motifs_fingerprint.sh
# Author : Arlet Lara
# Updated: 2025-11-13
# Description: Driver script to run peak calling, QC summaries, motif discovery,
#          and fingerprinting for ChIP-seq treated vs input pairs on SLURM.
#
# Overview:
#     1) Activates a conda environment on the submission node.
#     2) Defines project/base paths and creates top-level output dirs.
#     3) Generates per-sample SLURM child jobs (.qsh) that run, per treated/input pair:
#          - MACS2 narrow peak calling (paired-end) vs matched input.
#          - Peak filtering by score (>70 on col 5 of narrowPeak).
#          - FRiP for treated and input against filtered peaks.
#          - Motif discovery with HOMER (mm10)from filtered (>70 score) MACS2 peaks.
#          - Summit BigBed and fold-change (treated/input) bigWig.
#     4) Submits all generated child jobs in batch.
#     5) Summarizes blacklist removal stats (before/after/percent removed).
#     6) Summarizes chrM and chrUn alignment composition per BAM.
#     7) Computes library scaling factors from raw BAMs.
#     8) Plots a multi-sample fingerprint using treated + input BAMs.
#
# Inputs (assumed to exist under $base):
#   - filtered/*_blacklisted.bam - Deduplicated, blacklist-filtered BAMs (treated & input).
#   - removed_dup/*.bam - Deduplicated BAMs before blacklist filtering.
#   - removed_dup/*_blacklisted_with_chrMUncl.bam - Intermediate BAMs after blacklist + chrM/chrUn filtering.
#   - External helper scripts (absolute paths used below):
#       frip_calculation.sh
#       get_scaling_factors.sh
#
# Outputs:
#   - $base/peak_calling_macs2/*                      : MACS2 peak sets + filtered peaks.
#   - $base/download/Allsamples_FRiP_raw.txt          : Appended FRiP table.
#   - $base/download/scaling_factors_fromfilteredbams.txt
#   - $base/download/blacklist_removal_stats.txt      : Before/after/percent removed.
#   - $base/download/alignment_stats.txt              : chrM / chrUn content.
#   - $base/download/mapped_stats.txt                 : Flag-based mapping stats.
#   - $base/download/fingerprints_all.png             : Multi-sample fingerprint plot.
#   - $base/download/motifs_in_macs2_peaks/<sample>/  : HOMER motif results.
#   - SLURM stdout/err logs for per-sample child jobs.
#   - bigwigs_nodup/<sample>_macs2peaks.bb            : Summit BigBed (see Notes).
#   - bigwigs_nodup/<sample>.bw / <input>.bw          : Raw coverage bigWigs.
#   - bigwigs_nodup/<sample>_normalized.bw            : CPM-normalized coverage.
#   - bigwigs_nodup/FC_<sample>_vs_IgG.bw             : FC bigWig (treated vs input).
#
# Requires:
#   - SLURM.
#   - Conda env "chipseq_env". Check env folder for details
#       MACS2, HOMER, deepTools (bamCoverage, bigwigCompare, plotFingerprint),
#       samtools, awk, perl, bedtools, UCSC bedToBigBed.
#   - HOMER mm10 genome.
#   - mm10 chrom.sizes and bigNarrowPeak.as for BigBed conversion.
#               
# Notes:
#   - treat_files[] and input_files[] must be index-aligned in pairs (same order & length).
#   - Score threshold (>70) is a heuristic for peak stringency; can be adjusted.
#   - Child jobs enforce `set -euo pipefail` and `module purge` for reproducibility.
#   - The directory and files in bigwigs_nodup/ must exist and be writable,
#     as BigBed and bigWig tracks are written there.
#   - The chrom sizes file must match the genome build used by MACS2 (mm10 here).
################################################################################


# Activate environment (per-node) and set project paths
source /path/to/myprofile/anaconda3/etc/profile.d/conda.sh
conda activate chipseq_env

# Project identifier
project=name_of_the_project

# Root path for this project 
base=/path/to/myprofile/data/$project

# Create directories to save outputs 
cd $base
mkdir -p peak_calling_macs2 download


###############################################################################
# Step 1 ▸ Generate per-sample .qsh jobs for MACS2 + FRiP + HOMER
###############################################################################

# List of treated and input samples (index-aligned)
treat_files=('sample1_HA' 'sample2_HA')
input_files=('sample1_IgG' 'sample2_IgG')


# Emit one child .qsh per treated/input pair; each job calls peaks, filters,
#      computes FRiP for treated & input vs treated peaks >70, finds motifs,
#      and derives browser tracks (summit BigBed, FC bigWig).
for i in "${!treat_files[@]}"; do
    sample=${treat_files[i]}
    input=${input_files[i]}
    qsh=$base/scripts/${sample}_Peak_Calling_QC_Motifs.qsh

cat <<EOF> $qsh
#!/bin/bash -ex
#SBATCH --job-name=Peaks_${sample}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80GB
#SBATCH --time=20:00:00
#SBATCH --output=/path/to/scratch/$project/Peaks_QC_Motifs_${sample}.out
#SBATCH --error=/path/to/scratch/$project/Peaks_QC_Motifs_${sample}.err
#SBATCH --mail-type=END

# Safety & reproducibility toggles: Fail fast on any error/undefined var; avoid hidden environment state.
set -euo pipefail
module purge

# Activate environment
source /path/to/myprofile/anaconda3/etc/profile.d/conda.sh
conda activate chipseq_env

cd $base

# Ensure HOMER mm10 support is present (no-op if already installed).
# For me I had to install every time otherwise would fail to run.
perl /path/to/myprofile/anaconda3/envs/chipseq_env/share/homer/.//configureHomer.pl -install mm10

# -----------------------------
# 1.1) Peak Calling with MACS2 (narrow peaks for TFs, paired-end)
# -----------------------------
# Call narrow peaks with matched control; q-value 0.01 is standard for TF ChIP.
macs2 callpeak \
    -t filtered/${sample}_blacklisted.bam \
    -c filtered/${input}_blacklisted.bam \
    -f BAMPE \
    -g mm \
    -n ${sample} \
    --outdir peak_calling_macs2 \
    --keep-dup auto \
    -q 0.01

# Retain high-confidence peaks: score (col 5) > 70. 
#Filters to a stringent subset for downstream FRiP/motifs
awk '\$5 > 70' peak_calling_macs2/${sample}_peaks.narrowPeak > peak_calling_macs2/${sample}_peaks_above70.narrowPeak

# -----------------------------
# 1.2) Peak Calling Browser tracks
# -----------------------------
### Cap scores
awk 'BEGIN{OFS="\t"}{ if(\$5>1000)\$5=1000; if(\$5<0)\$5=0; print }' "peak_calling_macs2/${sample}_peaks_above70.narrowPeak" > "peak_calling_macs2/${sample}_peaks_above70_fixed.narrowPeak"

### Sort 
cat peak_calling_macs2/${sample}_peaks_above70_fixed.narrowPeak | sort -k1,1 -k2,2n > peak_calling_macs2/${sample}_peaks_above70_fixed_sorted.narrowPeak

### bed to BigBed to visualize peak regions
bedToBigBed \
-as=/path/to/myprofile/data/Common_files/bigNarrowPeak.as  \
-type=bed6+4 \
peak_calling_macs2/${sample}_peaks_above70_fixed_sorted.narrowPeak \
/path/to/myprofile/data/Common_files/mm10.chrom.sizes \
bigwigs_nodup/${sample}_macs2peaks.bb


# -----------------------------
# 2) FRiP (treated and input). Evaluate signal enrichment over the filtered peak set (>70).
# -----------------------------
bash /path/to/myprofile/data/Common_files/frip_calculation.sh \
    filtered/${sample}_blacklisted.bam \
    peak_calling_macs2/${sample}_peaks_above70.narrowPeak \
    download/Allsamples_FRiP_raw.txt

bash /path/to/myprofile/data/Common_files/frip_calculation.sh \
    filtered/${input}_blacklisted.bam \
    peak_calling_macs2/${sample}_peaks_above70.narrowPeak \
    download/Allsamples_FRiP_raw.txt


# -----------------------------
# 3) Motif discovery with HOMER (mm10; default peak size 500 bp)
# -----------------------------
findMotifsGenome.pl peak_calling_macs2/${sample}_peaks_above70.narrowPeak \
    mm10 \
    download/motifs_in_macs2_peaks/${sample}

# -----------------------------
# 4) UCSC track lines for BigBed and BigWigs
# -----------------------------
# Collect browser-ready track definitions for this sample in a unified file

# Peak regions (MACS2 peaks as BigBed)
echo "track name=${sample}_macs2_peak_regions bigDataUrl=\${base#/mnt/}/bigwigs_nodup/${sample}_macs2peaks.bb type=bigBed color=0,0,0 visibility=pack" >> download/ucsc_tracks_bigwigs.txt

# CPM-normalized coverage
echo "track name=${sample}_CPM bigDataUrl=https://informaticsdata.liai.org/\${base#/mnt/}/bigwigs_nodup/${sample}_normalized.bw type=bigWig color=31,77,92 visibility=hide" >> download/ucsc_tracks_bigwigs.txt
echo "track name=${input}_CPM bigDataUrl=https://informaticsdata.liai.org/\${base#/mnt/}/bigwigs_nodup/${input}_normalized.bw type=bigWig color=59,125,118 visibility=hide" >> download/ucsc_tracks_bigwigs.txt

# Fold-change bigWig (treated vs input)
echo "track name=FC_${sample}_vs_IgG bigDataUrl=\${base#/mnt/}/bigwigs_nodup/FC_${sample}_vs_IgG.bw type=bigWig color=0,109,119 visibility=hide" >> download/ucsc_tracks_bigwigs.txt

# Raw coverage
echo "track name=${sample} bigDataUrl=https://informaticsdata.liai.org/\${base#/mnt/}/bigwigs_nodup/${sample}.bw type=bigWig color=113,60,124 visibility=hide" >> download/ucsc_tracks_bigwigs.txt
echo "track name=${input} bigDataUrl=https://informaticsdata.liai.org/\${base#/mnt/}/bigwigs_nodup/${input}.bw type=bigWig color=129,109,156 visibility=hide" >> download/ucsc_tracks_bigwigs.txt


# -----------------------------
# 5) Bigwig coverage and fold-change tracks
# -----------------------------
# BigWig coverage tracks unnormalized, no adjustments
bamCoverage -p 16 -b filtered/${sample}_blacklisted.bam -o bigwigs_nodup/${sample}.bw \
--binSize 1 \

# BigWig coverage tracks no normalized or smoothed
bamCoverage -p 16 -b filtered/${input}_blacklisted.bam -o bigwigs_nodup/${input}.bw \
--binSize 1 \

# BigWig coverage tracks normalized, centered and smoothed reads
bamCoverage -p 16 -b filtered/${sample}_blacklisted.bam -o bigwigs_nodup/${sample}_normalized.bw \
--normalizeUsing CPM \
--binSize 5 \
--smoothLength 60 \
--centerReads

bamCoverage -p 16 -b filtered/${input}_blacklisted.bam -o bigwigs_nodup/${input}_normalized.bw \
--normalizeUsing CPM \
--binSize 5 \
--smoothLength 60 \
--centerReads


# Fold-change (ChIP vs control) bigWig to directly visualize enrichment
bigwigCompare \
  -b1 bigwigs_nodup/${sample}.bw \
  -b2 bigwigs_nodup/${input}.bw \
  --operation ratio \
  -o bigwigs_nodup/FC_${sample}_vs_IgG.bw \
  --binSize 1

EOF

done

###############################################################################
# Step 2 ▸ Submit all .qsh jobs
###############################################################################
cd $base
for file in scripts/*_Peak_Calling_QC_Motifs.qsh; do
    sbatch $file;
done

###############################################################################
# Step 3 ▸ Blacklist removal stats. Quantify reads removed by blacklist filtering (QC check)
###############################################################################

echo -e "Sample\tBefore\tAfter\tRemoved\tPercent_removed" > download/blacklist_removal_stats.txt
for bam in removed_dup/*[0-9].bam; do
    sample=$(basename "$bam" .bam)
    before=$(samtools view -c "$bam")
    after=$(samtools view -c "removed_dup/${sample}_blacklisted_with_chrMUncl.bam")
    removed=$((before - after))
    perc=$(awk -v b=$before -v r=$removed 'BEGIN{printf "%.2f", (r/b)*100}')
    echo -e "${sample}\t${before}\t${after}\t${removed}\t${perc}" >> download/blacklist_removal_stats.txt
done


###############################################################################
# Step 4 ▸ Mapping flag stats. Summarize read categories (primary/secondary/supplementary/unmapped)
###############################################################################
echo -e "file\tall\tprimary_mapped\tprimary_unmapped\tsecondary\tsupplementary" > "download/mapped_stats.txt"

for b in removed_dup/*_blacklisted_with_chrMUncl.bam; do
  ALL=$(samtools view -c "$b")
  PRIMARY_MAPPED=$(samtools view -c -F 0x904 "$b")         # exclude 0x4,0x100,0x800
  PRIMARY_UNMAPPED=$(samtools view -c -f 0x4 -F 0x900 "$b") # unmapped, not sec/supp
  SECONDARY=$(samtools view -c -f 0x100 "$b")
  SUPPLEMENTARY=$(samtools view -c -f 0x800 "$b")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$b" "$ALL" "$PRIMARY_MAPPED" "$PRIMARY_UNMAPPED" "$SECONDARY" "$SUPPLEMENTARY" >> "download/mapped_stats.txt"
done


##############################################################################
# Step 5 ▸ chrM / chrUn composition. Mitochondrial and unplaced reads are common QC metrics.
###############################################################################

echo -e "Sample\tTotal\tchrM\tchrM_pct\tchrUn\tchrUn_pct" > download/alignment_stats.txt

for bam in removed_dup/*_blacklisted_with_chrMUncl.bam; do
  sample=$(basename "$bam" .bam)
  samtools index -@ 8 "$bam"

  total=$(samtools view -c -F 0x904 "$bam")
  chrM=$(samtools view -c -F 0x904 "$bam" chrM)
  chrUn=$(samtools idxstats "$bam" | awk '$1 ~ /^chrUn/ {s+=$3} END{print s+0}')

  chrM_pct=$(awk -v t=$total -v m=$chrM 'BEGIN{printf "%.3f", 100*(m/t)}')
  chrUn_pct=$(awk -v t=$total -v u=$chrUn 'BEGIN{printf "%.3f", 100*(u/t)}')

  echo -e "${sample}\t${total}\t${chrM}\t${chrM_pct}\t${chrUn}\t${chrUn_pct}" >> download/alignment_stats.txt
done

###############################################################################
# Step 6 ▸ Library Size Factors. Derive per-sample scaling factors for downstream normalization.
###############################################################################

# Path to filtered BAMs
bam_dir=$base/filtered

# Output text file for scaling factors
scaling_file=$base/download/scaling_factors_fromfilteredbams.txt

# BAM suffix
bam_suffix="_blacklisted.bam"

# Run the scaling factor script (external helper)
/path/to/myprofile/data/Common_files/get_scaling_factors.sh "$bam_dir" "$scaling_file" "$bam_suffix"

###############################################################################
# Step 7 ▸ Fingerprint plotting. Assess relative enrichment across all samples
###############################################################################

# Build the BAM paths from the basenames:
bams=(
  "${treat_files[@]/%/_blacklisted.bam}"
  "${input_files[@]/%/_blacklisted.bam}" )     # 1) append suffix to each element
bams=("${bams[@]/#/filtered/}")                # 2) prepend directory to each element

plotFingerprint \
-b "${bams[@]}" \
--smartLabels \
--plotTitle "Fingerprint: aHA + IgG" \
-plot download/fingerprints_all.png












