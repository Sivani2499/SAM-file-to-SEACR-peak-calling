# Input:  *_clean.sam files (human-classified reads from Kraken2)
# Output: SEACR peak BED files (.relaxed.bed) + FRiP QC scores
# =============================================================================

# =============================================================================
# CONFIG — UPDATE THESE PATHS BEFORE RUNNING
# =============================================================================

INPUT_DIR="/path/to/kraken/human_classified_sam"        # SAM files from Kraken2
BAM_DIR="/path/to/output/bam"                           # BAM output directory
BEDGRAPH_DIR="/path/to/output/bedgraph"                 # bedGraph output directory
SEACR_OUT="/path/to/output/seacr"                       # SEACR peak output
FRIP_OUT="/path/to/output/qc_frip"                      # FRiP QC output
SEACR_SCRIPT="/path/to/tools/SEACR/SEACR_1.3.sh"       # Path to SEACR_1.3.sh - the directory where you have cloned SEACR
IGG_BEDGRAPH="${BEDGRAPH_DIR}/IgG_sample.norm.bedgraph" # IgG control bedGraph

# =============================================================================
# STEP 1: SAM to BAM — Convert, Sort, and Index
# =============================================================================
# Why: SAM files are human-readable but large and slow. BAM is the compressed
#      binary equivalent. Coordinate sorting is required for downstream tools.
#
# Input:  *_clean.sam (human-classified reads from Kraken2)
# Output: *.sorted.bam + *.sorted.bam.bai (index)

echo "=== STEP 1: SAM to sorted BAM ==="

mkdir -p "${BAM_DIR}"

for sam in "${INPUT_DIR}"/*_clean.sam; do
  base=$(basename "$sam" _clean.sam)

  # Convert SAM to BAM
  samtools view -bS "$sam" -o "${BAM_DIR}/${base}.bam"

  # Sort by genomic coordinate
  samtools sort -o "${BAM_DIR}/${base}.sorted.bam" "${BAM_DIR}/${base}.bam"

  # Index for random access
  samtools index "${BAM_DIR}/${base}.sorted.bam"

  # Remove unsorted BAM to save space
  rm "${BAM_DIR}/${base}.bam"

  echo "Done: ${base}"
done

# QC check
echo "--- Flagstat check after Step 1 ---"
for bam in "${BAM_DIR}"/*.sorted.bam; do
  echo "=== $(basename $bam) ==="
  samtools flagstat "$bam" | grep -E "in total|mapped \(|properly paired"
  echo ""
done

# =============================================================================
# STEP 2: Duplicate Removal (Mate-Aware)
# =============================================================================
# Why: PCR duplicates artificially inflate signal. SAMtools markdup is used
#      here because it is MATE-AWARE — it uses MC (mate CIGAR) and MS (mate
#      score) tags added by fixmate to correctly identify duplicate read pairs.
#      This is more accurate than Picard MarkDuplicates for paired-end data.
#
# The process requires 4 sub-steps in this exact order:
#   sorted.bam → namesort.bam → fixmate.bam → fixmate.sorted.bam → rmdup.bam
#
# Input:  *.sorted.bam
# Output: dedup/*.rmdup.bam + dedup/*.markdup.stats.txt

echo "=== STEP 2: Duplicate Removal ==="

mkdir -p "${BAM_DIR}/dedup"

# --- Step 2a: Name Sort ---
# Why: fixmate requires reads sorted by NAME so read pairs are adjacent

echo "--- Step 2a: Name sort ---"
for bam in "${BAM_DIR}"/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  samtools sort -n \
    -o "${BAM_DIR}/${base}.namesort.bam" \
    "$bam"
  echo "Name sorted: ${base}"
done

# --- Step 2b: Fixmate ---
# Why: Adds MC and MS tags so markdup can identify duplicates at the
#      pair level, not just the individual read level

echo "--- Step 2b: Fixmate ---"
for bam in "${BAM_DIR}"/*.namesort.bam; do
  base=$(basename "$bam" .namesort.bam)
  samtools fixmate -m \
    "$bam" \
    "${BAM_DIR}/${base}.fixmate.bam"
  echo "Fixmate done: ${base}"
done

# --- Step 2c: Coordinate Sort (after fixmate) ---
# Why: markdup requires coordinate-sorted input

echo "--- Step 2c: Coordinate sort after fixmate ---"
for bam in "${BAM_DIR}"/*.fixmate.bam; do
  base=$(basename "$bam" .fixmate.bam)
  samtools sort \
    -o "${BAM_DIR}/${base}.fixmate.sorted.bam" \
    "$bam"
  echo "Coord sorted: ${base}"
done

# --- Step 2d: Mark and Remove Duplicates ---
# Why: -r physically removes duplicates (not just flags them)
#      -s prints statistics to a .txt file for QC

echo "--- Step 2d: Mark and remove duplicates ---"
for bam in "${BAM_DIR}"/*.fixmate.sorted.bam; do
  base=$(basename "$bam" .fixmate.sorted.bam)

  samtools markdup -r -s \
    "$bam" \
    "${BAM_DIR}/dedup/${base}.rmdup.bam" \
    2> "${BAM_DIR}/dedup/${base}.markdup.stats.txt"

  samtools index "${BAM_DIR}/dedup/${base}.rmdup.bam"

  echo "Dedup done: ${base}"
done

# QC: Check duplicate stats
echo "--- Duplicate removal stats ---"
for f in "${BAM_DIR}/dedup"/*.markdup.stats.txt; do
  echo "=== $(basename $f .markdup.stats.txt) ==="
  grep -E "READ|WRITTEN|DUPLICATE PRIMARY TOTAL" "$f"
  echo ""
done

# QC: Confirm 0 duplicates remain in output BAM
echo "--- Flagstat check after dedup (duplicates should be 0) ---"
for bam in "${BAM_DIR}/dedup"/*.rmdup.bam; do
  echo "=== $(basename $bam) ==="
  samtools flagstat "$bam" | grep -E "in total|duplicates|mapped \(|properly paired"
  echo ""
done

# =============================================================================
# STEP 3: BAM to Fragment BED
# =============================================================================
# Why: CUT&RUN signal should be counted per FRAGMENT (the actual DNA molecule
#      cut by MNase), not per READ. Each fragment = 2 reads in paired-end
#      sequencing. Using bedtools bamtobed -bedpe extracts the full fragment
#      coordinates so each fragment is counted exactly once. This avoids the
#      2x signal inflation that occurs when counting reads individually.
#
# Pipeline:
#   rmdup.bam → namesort.bam → BEDPE → filter artifacts
#   → extract chr/start/end → remove bad rows → fragments.clean.bed
#
# Filters applied:
#   1. Both reads on same chromosome ($1==$4) — removes chimeric pairs
#   2. Fragment size < 1000bp ($6-$2 < 1000) — removes incomplete digestion
#   3. No "." chromosome, no negative coordinates — removes unmapped artifacts
#
# Input:  dedup/*.rmdup.bam
# Output: bedgraph/*.fragments.clean.bed

echo "=== STEP 3: BAM to Fragment BED ==="

mkdir -p "${BEDGRAPH_DIR}"

for bam in "${BAM_DIR}/dedup"/*.rmdup.bam; do
  base=$(basename "$bam" .rmdup.bam)

  # Name sort required for BEDPE conversion
  samtools sort -n "$bam" \
    -o "${BEDGRAPH_DIR}/${base}.namesort.bam"

  # Convert paired BAM to BEDPE (one line per fragment pair)
  bedtools bamtobed -bedpe \
    -i "${BEDGRAPH_DIR}/${base}.namesort.bam" \
    > "${BEDGRAPH_DIR}/${base}.bed"

  # Filter artifacts
  awk '$1==$4 && $6-$2 < 1000 {print $0}' \
    "${BEDGRAPH_DIR}/${base}.bed" \
    > "${BEDGRAPH_DIR}/${base}.clean.bed"

  # Extract chr, fragment_start, fragment_end and sort
  cut -f 1,2,6 "${BEDGRAPH_DIR}/${base}.clean.bed" \
    | sort -k1,1 -k2,2n \
    > "${BEDGRAPH_DIR}/${base}.fragments.bed"

  # Remove rows where chromosome is "." or coordinates are negative
  awk '$1 != "." && $2 >= 0 && $3 >= 0' \
    "${BEDGRAPH_DIR}/${base}.fragments.bed" \
    > "${BEDGRAPH_DIR}/${base}.fragments.clean.bed"

  # Clean up intermediate files
  rm "${BEDGRAPH_DIR}/${base}.namesort.bam"
  rm "${BEDGRAPH_DIR}/${base}.bed"
  rm "${BEDGRAPH_DIR}/${base}.clean.bed"
  rm "${BEDGRAPH_DIR}/${base}.fragments.bed"

  echo "Done: ${base}"
done

# QC: Fragment counts per sample
echo "--- Fragment counts per sample ---"
for f in "${BEDGRAPH_DIR}"/*.fragments.clean.bed; do
  echo "$(basename $f .fragments.clean.bed): $(wc -l < $f) fragments"
done

# =============================================================================
# STEP 4: Generate Chromosome Sizes
# =============================================================================
# Why: bedtools genomecov requires a chromosome sizes file to know the
#      boundaries of each chromosome. We extract this directly from the BAM
#      header to ensure it matches the exact reference used in alignment.
#      Important: this pipeline uses RefSeq NC_ chromosome naming (not chr1),
#      so we cannot use a standard UCSC chrom.sizes file.
#
# Input:  Any one rmdup.bam (all were aligned to the same reference)
# Output: hg38.refseq.chrom.sizes

echo "=== STEP 4: Generate chromosome sizes ==="

REF_BAM=$(ls "${BAM_DIR}/dedup"/*.rmdup.bam | head -1)

samtools view -H "${REF_BAM}" \
  | grep '^@SQ' \
  | sed 's/@SQ\tSN://; s/\tLN:/\t/' \
  | awk '{print $1"\t"$2}' \
  > "${BEDGRAPH_DIR}/hg38.refseq.chrom.sizes"

echo "Chromosome sizes saved to: ${BEDGRAPH_DIR}/hg38.refseq.chrom.sizes"
echo "First 5 entries:"
head -5 "${BEDGRAPH_DIR}/hg38.refseq.chrom.sizes"

# =============================================================================
# STEP 5: bedGraph Generation + CPM Normalization
# =============================================================================
# Why: bedtools genomecov converts fragment positions into a genome-wide
#      coverage track (bedGraph format). CPM normalization (Counts Per Million
#      fragments) scales the signal so samples with different sequencing depths
#      can be directly compared.
#
# CPM formula:
#   scale_factor = 1,000,000 / total_fragments
#   normalized_signal = raw_signal x scale_factor
#
#
# Input:  *.fragments.clean.bed
# Output: *.raw.bedgraph + *.norm.bedgraph + fragment_counts.txt

echo "=== STEP 5: bedGraph generation and CPM normalization ==="

CHROMSIZES="${BEDGRAPH_DIR}/hg38.refseq.chrom.sizes"

# Count fragments per sample and save to file for use in normalization loop
echo "--- Counting fragments per sample ---"
for f in "${BEDGRAPH_DIR}"/*.fragments.clean.bed; do
  base=$(basename "$f" .fragments.clean.bed)
  count=$(wc -l < "$f")
  echo -e "${base}\t${count}"
done > "${BEDGRAPH_DIR}/fragment_counts.txt"

cat "${BEDGRAPH_DIR}/fragment_counts.txt"

# Generate raw and CPM-normalized bedGraphs
echo "--- Generating bedGraphs ---"
while read sample fragments; do

  # Raw bedGraph — one line per genomic interval with coverage > 0
  bedtools genomecov \
    -bg \
    -i "${BEDGRAPH_DIR}/${sample}.fragments.clean.bed" \
    -g "${CHROMSIZES}" \
    > "${BEDGRAPH_DIR}/${sample}.raw.bedgraph"

  # CPM scale factor
  scale=$(echo "1000000 / $fragments" | bc -l)

  # Apply CPM normalization to column 4 (the signal/count column)
  awk -v s="$scale" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*s}' \
    "${BEDGRAPH_DIR}/${sample}.raw.bedgraph" \
    > "${BEDGRAPH_DIR}/${sample}.norm.bedgraph"

  echo "Done: ${sample} | fragments: ${fragments} | CPM scale factor: ${scale}"

done < "${BEDGRAPH_DIR}/fragment_counts.txt"

# QC: Verify all bedGraphs were created
echo "--- bedGraph file counts (should match number of samples) ---"
echo "Raw bedGraphs:  $(ls ${BEDGRAPH_DIR}/*.raw.bedgraph | wc -l)"
echo "Norm bedGraphs: $(ls ${BEDGRAPH_DIR}/*.norm.bedgraph | wc -l)"

# QC: Max CPM signal per sample
# Note: IgG max signal should be much lower than STAT3 samples
echo "--- Max CPM signal per sample ---"
for f in "${BEDGRAPH_DIR}"/*.norm.bedgraph; do
  max=$(awk 'BEGIN{max=0} {if($4>max) max=$4} END{print max}' "$f")
  echo "$(basename $f .norm.bedgraph): max CPM = $max"
done

# =============================================================================
# STEP 6: SEACR Peak Calling
# =============================================================================
# Why: SEACR (Sparse Enrichment Analysis for CUT&RUN) is designed specifically
#      for CUT&RUN data.This pipeline uses the
#      IgG control bedGraph as background — only calling peaks where signal
#      genuinely exceeds non-specific antibody binding.
#
# Parameters:
#   norm     : normalize signal between treatment and IgG before comparison
#   relaxed  : more inclusive peak calling — recommended for low-depth TF
#              samples and for FRiP calculation
#
# Note on IgG:
#   IgG is used as the CONTROL for all samples — it is NOT run through SEACR
#   as a sample itself. If IgG produces no .relaxed.bed file, that is expected
#   and confirms background subtraction is working correctly.
#
# Input:  bedgraph/*.norm.bedgraph (treatment + IgG control)
# Output: seacr/*.relaxed.bed + seacr/*.stringent.bed

echo "=== STEP 6: SEACR Peak Calling ==="

mkdir -p "${SEACR_OUT}"

IGG_BEDGRAPH="${BEDGRAPH_DIR}/MCIA0354-Total-T-IgG.norm.bedgraph"

# --- Run SEACR for all STAT3 samples ---
# Uses glob *STAT3* to match all IL6 and Non-stimulated STAT3 samples
# across both patients (MCIA0354 and MCUI2217A)

echo "--- Running SEACR for all STAT3 samples ---"
for sig in "${BEDGRAPH_DIR}"/*STAT3*.norm.bedgraph; do
  base=$(basename "$sig" .norm.bedgraph)

  bash "${SEACR_SCRIPT}" \
    "$sig" \
    "${IGG_BEDGRAPH}" \
    norm \
    relaxed \
    "${SEACR_OUT}/${base}"

  echo "Done: ${base}"
done

# --- Run SEACR for H3 positive control ---
# H3 is a histone mark — expected to show broad, high-signal peaks
# High FRiP for H3 confirms the wet lab CUT&RUN protocol worked

echo "--- Running SEACR for H3 positive control ---"
bash "${SEACR_SCRIPT}" \
  "${BEDGRAPH_DIR}/MCIA0354-Total-T-H3.norm.bedgraph" \
  "${IGG_BEDGRAPH}" \
  norm \
  relaxed \
  "${SEACR_OUT}/MCIA0354-Total-T-H3"

echo "Done: MCIA0354-Total-T-H3"

# --- Run SEACR for Input control ---
# Input represents total open chromatin background
# Expected to have moderate peaks — not zero, not as specific as STAT3

echo "--- Running SEACR for Input control ---"
bash "${SEACR_SCRIPT}" \
  "${BEDGRAPH_DIR}/MCIA0354-Total-T-Input.norm.bedgraph" \
  "${IGG_BEDGRAPH}" \
  norm \
  relaxed \
  "${SEACR_OUT}/MCIA0354-Total-T-Input"

echo "Done: MCIA0354-Total-T-Input"

# QC: Count peaks per sample
echo "--- Peak counts per sample ---"
for bed in "${SEACR_OUT}"/*.relaxed.bed; do
  echo "$(basename $bed .relaxed.bed): $(wc -l < $bed) peaks"
done


# =============================================================================
# STEP 7: FRiP Score Calculation
# =============================================================================
# Why: FRiP (Fraction of Reads in Peaks) is a standard QC metric for
#      chromatin profiling. It measures what proportion of sequencing signal
#      falls inside called peaks (signal) vs. the rest of the genome (noise).
#
# Formula:
#   FRiP = reads overlapping peaks / total mapped reads
#
# Expected values for CUT&RUN:
#   > 0.30       : Excellent
#   0.10 - 0.30  : Acceptable
#   0.05 - 0.10  : Low but usable
#   < 0.05       : Poor — high background
#
# Key biological check:
#   IL6-stimulated FRiP > Non-stimulated FRiP
#   → confirms STAT3 is specifically activated by IL6 signaling
#   → if this direction is wrong, the experiment has failed
#
# Note on IgG:
#   IgG will have no .relaxed.bed file (0 peaks called).
#   The script handles this gracefully and reports FRiP = 0 for IgG.
#   
#
# Note 
#   Some samples may have low sequencing depth after Kraken2 filtering
#   (0.4M - 1.5M fragments). SEACR may call 0 peaks for the lowest
#   depth samples. FRiP = 0 for these does not mean the experiment
#   failed — it reflects insufficient depth for reliable peak calling.
#
# Input:  dedup/*.rmdup.bam + seacr/*.relaxed.bed
# Output: qc_frip/frip_scores.txt

echo "=== STEP 7: FRiP Score Calculation ==="

mkdir -p "${FRIP_OUT}"
FRIP_FILE="${FRIP_OUT}/frip_scores.txt"

echo -e "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP" > "${FRIP_FILE}"

for bam in "${BAM_DIR}/dedup"/*.rmdup.bam; do
  base=$(basename "$bam" .rmdup.bam)

  # Total mapped reads
  total=$(samtools flagstat "$bam" \
    | grep "mapped (" \
    | head -1 \
    | awk '{print $1}')

  peak="${SEACR_OUT}/${base}.relaxed.bed"

   # Handle missing peak file gracefully
  # This is expected for:
  #   - IgG (used as control, not as sample — 0 peaks is correct)
  #   - Very low depth samples where SEACR found no enriched regions
  if [ ! -f "$peak" ]; then
    echo "NOTE: No peak file for ${base} — 0 peaks called (expected for IgG or very low depth samples)"
    echo -e "${base}\t${total}\t0\t0.0000" | tee -a "${FRIP_FILE}"
    echo ""
    continue
  fi

  # Count reads overlapping any peak region
  # -u : count each read only ONCE even if it overlaps multiple peaks
  in_peaks=$(bedtools intersect \
    -a "$bam" \
    -b "$peak" \
    -u | samtools view -c)

  # FRiP to 4 decimal places
  frip=$(echo "scale=4; $in_peaks / $total" | bc)

  echo -e "${base}\t${total}\t${in_peaks}\t${frip}" | tee -a "${FRIP_FILE}"

done
echo ""
echo "FRiP scores saved to: ${FRIP_FILE}"

echo ""
echo "================================================================="
echo "  PIPELINE COMPLETE"
echo "================================================================="
echo ""
echo "Output summary:"
echo "  Deduplicated BAMs    : ${BAM_DIR}/dedup/*.rmdup.bam"
echo "  Dup stats            : ${BAM_DIR}/dedup/*.markdup.stats.txt"
echo "  Fragment BEDs        : ${BEDGRAPH_DIR}/*.fragments.clean.bed"
echo "  Fragment counts      : ${BEDGRAPH_DIR}/fragment_counts.txt"
echo "  Raw bedGraphs        : ${BEDGRAPH_DIR}/*.raw.bedgraph"
echo "  Norm bedGraphs       : ${BEDGRAPH_DIR}/*.norm.bedgraph"
echo "  Chrom sizes          : ${BEDGRAPH_DIR}/hg38.refseq.chrom.sizes"
echo "  SEACR peaks          : ${SEACR_OUT}/*.relaxed.bed"
echo "  FRiP scores          : ${FRIP_FILE}"
echo "================================================================="
