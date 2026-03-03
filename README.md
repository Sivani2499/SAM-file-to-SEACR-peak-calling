# SAM-file-to-SEACR-peak-calling
# 
# CUT&RUN Post-Alignment Pipeline: SAM to Peak Calling
# 
# Author: Sivani Ravindran
# Description:
   This script takes human-classified SAM files from Kraken2 output and
   processes them through BAM conversion, duplicate removal, fragment BED
   generation, CPM-normalized bedGraph creation, and SEACR peak calling.

# Usage:
   bash cutandrun_peak_calling_pipeline.sh
   (update the directory variables in the CONFIG section before running)
#
# Input:  *_clean.sam files (human-classified reads from Kraken2)
# Output: SEACR peak BED files (.relaxed.bed) + FRiP QC scores
