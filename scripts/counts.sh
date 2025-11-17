#!/bin/bash -ex
#SBATCH --job-name=counts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --time=24:00:00
#SBATCH --output=/path/to/scratch/name_of_the_project/counts.out
#SBATCH --error=/path/to/scratch/name_of_the_project/counts.err
#SBATCH --mail-type=END

###############################################################################
# counts.sh ▸ Summarize ChIP-seq peak sizes and extract counts per mark
#
# Overview:
#   1) Activates a conda environment on the compute node.
#   2) Defines project/base paths and the ChIP mark (${mark}), and creates
#      top-level output directories.
#   3) Computes per-sample peak-width statistics for MACS2 peaks (score ≥ 70)
#      for the specified mark.
#   4) Builds a nonredundant candidate peak set by merging overlapping peaks
#      across condition1 and condition2 replicates for the mark.
#   5) Runs multiBamSummary on:
#         - merged candidate peaks (region-level coverage).
#         - one representative 1-bp summit per candidate peak.
#         - fixed-size genomic bins outside the candidate
#          peak regions to obtain background coverage.
#    6) Calculate FRiP in candidate peaks
#
# Inputs (assumed to exist under $base):
#   - peak_calling_macs2/*${mark}*_peaks_above70.*Peak
#       MACS2 peak calls (score-filtered upstream) for the chosen mark.
#   - peak_calling_macs2/*_summits.bed
#       MACS2 summits (1-bp) corresponding to filtered peaks.
#   - filtered/*blacklisted.bam
#       Deduplicated, blacklist-filtered BAMs (treated & input).
#   - Tools that must be available in the active environment:
#       bedtools
#       multiBamSummary (deepTools)
#     (can just use chipseq_env)
#
#
# Outputs:
#   - download/peak_size_summary.txt
#       Per-peak-file width statistics (mean, min, max, SD, count).
#   - peak_calling_macs2/condition_${mark}_reps_peaks_merged_withNamesScore.bed
#       Merged peaks for each condition with collapsed names and max score.
#   - peak_calling_macs2/condition_${mark}_reps_peaks_overlapOnly_withNamesScore.bed
#       Subset of merged peaks whithin each condition supported by ≥2 peaks.
#   - peak_calling_macs2/candidate_${mark}_peaks_merged_withNamesScore.bed
#       Final candidate peak set used for counting and background definition.
#   - counts/multibamsummary_coverage_${mark}.npz
#   - download/counts/multibamsummary_coverage_${mark}.tab
#       multiBamSummary results on candidate peaks.
#   - peak_calling_macs2/all_summits.bed
#       Concatenated summits (with normalized name/score columns).
#   - peak_calling_macs2/summits_in_candidate_${mark}_peaks_merged.tsv
#       Summits intersecting candidate peaks, with candidate information.
#   - peak_calling_macs2/unique_summits_in_candidate_${mark}_peaks.bed
#       One highest-score summit per candidate peak.
#   - counts/multibamsummary_summits.npz
#   - download/counts/multibamsummary_summits.tab
#       multiBamSummary results at unique 1-bp summits.
#   - counts/multibamsummary_bins_noPeaks.npz
#   - download/counts/multibamsummary_bins_noPeaks.tab
#       multiBamSummary results on fixed-size bins outside candidate peaks.
#
#
# Notes:
#   - The `mark` variable controls which ChIP target (aHA/aFLAG/H3k27ac/etc.)
#     is used in file-matching patterns and output filenames.
#   - Project directory and SLURM output/error paths are specific to the
#     project and should be updated.
###############################################################################

# Activate environment
source /path/to/myprofile/anaconda3/etc/profile.d/conda.sh
conda activate chipseq_env

# Project identifier 
project=name_of_the_project

# Mark (epitope / histone modification) to process.
mark="aHA" # Replace with aFLAG/aHA/H3k27ac depending the dataset

# Assumes raw fastqs exist under $source/raw_data/
base=/path/to/myprofile/data/$project

# Create directories
cd $base
mkdir -p download/counts/ counts/


###############################################################################
# Step 1 ▸ Peak stats. Summarizing peak width distributions per sample.
###############################################################################

# Initialize summary table with header
echo -e "Sample\tMean_bp\tMin_bp\tMax_bp\tSD_bp\tCount" > "download/peak_size_summary.txt"

# Loop over all peak files for the selected mark with score ≥ 70
for np in peak_calling_macs2/*${mark}*_peaks_above70.*Peak; do
  sample=$(basename "$np")

  # Compute stats in awk and capture as tab-separated fields:
  #   - width = (end - start)
  #   - min/max, mean, and sample SD of widths
  #   - handle empty files by outputting NA and count=0.
  stats=$(awk '
    BEGIN{min=""; max=0}
    { w=$3-$2; sum+=w; sumsq+=w*w; n++; if(min==""||w<min) min=w; if(w>max) max=w }
    END{
      if(n==0){ print "NA\tNA\tNA\tNA\t0"; exit }
      mean = sum/n
      sd = (n>1)? sqrt((sumsq - (sum*sum)/n)/(n-1)) : 0
      printf "%.2f\t%d\t%d\t%.2f\t%d", mean, min, max, sd, n
    }' "$np")

  # Append per-sample stats row to the summary file
  echo -e "${sample}\t${stats}" >> "download/peak_size_summary.txt"
done


###############################################################################
# Step 2 ▸ Select peaks of interest (candidates)
###############################################################################

#### Find common (replicable) peaks within conditions to avoid false positives
# Step 2.1 ▸ Merge files whithin conditions 
#   1. cat: concatenate all MACS2 Peak files from a condition for this
#      mark with score ≥ 70 (score filter applied upstream).
#   2. sort: order by chromosome (field 1) then numeric start (field 2),
#      as required by bedtools merge.
#   3. bedtools merge:
#        - merges overlapping intervals across all conditionX files into
#          nonredundant intervals.
#        - -c 4,5 -o collapse,max:
#            • col4 (peak name): collapse contributing names into comma list.
#            • col5 (score): keep the maximum score among merged peaks.
#   4. awk -F"\t" '$4 ~ /,/':
#        - keep only rows where col4 contains a comma (≥2 peaks contributed).
#        - This enforces that at least two intervals overlapped, but it does
#          NOT ensure they are from distinct files/replicates.
#
# Important:
#   - If you have four replicates and at least two have overlapping peaks, that
#     region is kept.
#   - If a single file has two overlapping intervals, that region is also kept.
#   - The merge step does not track which input file each interval came from,  
#     it just requires two overlapping ranges from the merge of the peaks across the files.


# Considition 1
cat peak_calling_macs2/condition1*${mark}*above70.*Peak | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5 -o collapse,max > peak_calling_macs2/condition1_${mark}_reps_peaks_merged_withNamesScore.bed
awk -F"\t" '$4 ~ /,/' peak_calling_macs2/condition1_${mark}_reps_peaks_merged_withNamesScore.bed > peak_calling_macs2/condition1_${mark}_reps_peaks_overlapOnly_withNamesScore.bed

# Considition 2
cat peak_calling_macs2/condition2_${mark}*above70.*Peak | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5 -o collapse,max > peak_calling_macs2/condition2_${mark}_reps_peaks_merged_withNamesScore.bed
awk -F"\t" '$4 ~ /,/' peak_calling_macs2/condition2_${mark}_reps_peaks_merged_withNamesScore.bed > peak_calling_macs2/condition2_${mark}_reps_peaks_overlapOnly_withNamesScore.bed

# Step 2.3 ▸ Build final candidate peak set
#   - Combine overlap-only peaks across conditions.
#   - Sort and merge again to obtain nonredundant candidate intervals with
#     collapsed names and max scores.
cat \
peak_calling_macs2/condition1_${mark}_reps_peaks_overlapOnly_withNamesScore.bed \
peak_calling_macs2/condition2_${mark}_reps_peaks_overlapOnly_withNamesScore.bed \
| awk 'BEGIN{OFS="\t"} {sub(/\r$/,"")} NF>=5 {print $1,$2,$3,$4,$5}' \
| bedtools sort -i - \
| bedtools merge -i - -c 4,5 -o collapse,max \
> peak_calling_macs2/candidate_${mark}_peaks_merged_withNamesScore.bed

###############################################################################
# Step 3 ▸ Get counts; quantify coverage
###############################################################################
# Set sets of files
candidates="peak_calling_macs2/candidate_${mark}_peaks_merged_withNamesScore.bed"
summits=(peak_calling_macs2/*_summits.bed)
bams=(filtered/*blacklisted.bam)

# ---------- A. On candidate peaks. Count reads per BAM in each merged candidate peak interval
multiBamSummary BED-file -p 8 \
--BED peak_calling_macs2/candidate_${mark}_peaks_merged_withNamesScore.bed \
--bamfiles "${bams[@]}" \
--outRawCounts download/counts/multibamsummary_coverage_${mark}.tab \
-o counts/multibamsummary_coverage_${mark}.npz

# ---------- B. On summits. Count at 1-bp summits (peak apex) restricted to candidate peaks.
# 1) Combine all 1-bp summits using columns: chr start end name score
: > peak_calling_macs2/all_summits.bed
for f in "${summits[@]}"; do
  awk 'BEGIN{OFS="\t"}{
         n = ($4?$4:"summit_"NR);     # keep name if present
         s = ($5?$5:0);               # keep score if present; else 0
         print $1,$2,$3,n,s
       }' "$f" >> peak_calling_macs2/all_summits.bed
done

# 2) Keep only summits that fall inside candidate ranges
bedtools intersect -wa -wb \
  -a peak_calling_macs2/all_summits.bed \
  -b "$candidates" \
  > peak_calling_macs2/summits_in_candidate_${mark}_peaks_merged.tsv


# 3) For each candidate range (c_chr,c_start,c_end), keep ONE summit, choosing 
#    the highest 's_score' one (if scores absent, this is arbitrary)
awk 'BEGIN{OFS="\t"}
     {
       key = $6 ":" $7 "-" $8;   # candidate range label
       score = ($5+0);
       if(!(key in best) || score > best[key]) {
         best[key] = score;
         # output: summit chr/start/end; name = candidate range string
         line[key] = $1 "\t" $2 "\t" $3 "\t" key;
       }
     }
     END{for(k in line) print line[k]}' \
  peak_calling_macs2/summits_in_candidate_${mark}_peaks_merged.tsv \
  > peak_calling_macs2/unique_summits_in_candidate_${mark}_peaks.bed

# 4) multiBamSummary on the unique 1-bp summits (counts at that base)
multiBamSummary BED-file -p 8 \
  --BED peak_calling_macs2/unique_summits_in_candidate_${mark}_peaks.bed \
  --bamfiles "${bams[@]}" \
  --outRawCounts download/counts/multibamsummary_summits.tab \
  -o counts/multibamsummary_summits.npz \
  --centerReads 


# ---------- C. On non-peak regions. Segment the genome in 10kb bins and exclude  candidate peak regions
multiBamSummary bins -p 8 \
  --bamfiles "${bams[@]}" \
  --blackListFileName peak_calling_macs2/candidate_${mark}_peaks_merged_withNamesScore.bed  \
  --binSize 10000 --distanceBetweenBins 0 \
  --outRawCounts download/counts/multibamsummary_bins_noPeaks.tab \
  -o counts/multibamsummary_bins_noPeaks.npz


###############################################################################
# Step 4 ▸ Calculate FRiP over peaks of interest
###############################################################################
for bam in filtered/*_blacklisted.bam; do
  
bash path/to/myprofile/data/Common_files/frip_calculation.sh \
    ${bam} \
    peak_calling_macs2/candidate_aFLAG_peaks_merged_withNamesScore.bed \
    download/Allsamples_FRiP_afterpeakselection.txt
done



















