# CUT&RUN Post-Alignment Pipeline
### SAM (Kraken2 output) → BAM → Deduplication → Fragment BED → bedGraph → SEACR Peak Calling → FRiP QC

---
Author : Sivani Ravindran
## Overview

This pipeline processes CUT&RUN paired-end sequencing data starting from
human-classified SAM files (output of Kraken2 contamination filtering)
through to peak calling and quality control.

The pipeline is designed around fragment-based signal quantification and
IgG-normalized peak calling — ensuring that signal reflects true protein
binding rather than sequencing artefacts or non-specific background.

---

## When to Use This Pipeline

Use this pipeline if you have:

- **CUT&RUN** or **CUT&Tag** paired-end sequencing data
- Already run **Kraken2** to remove non-human contamination
  (e.g. E. coli carry-over from pA-MNase enzyme)
- Human-classified reads saved as `*_clean.sam` files
- An **IgG negative control** sample for peak calling background
- Samples aligned to **hg38** with **RefSeq chromosome naming**
  (e.g. `NC_000001.11`, not `chr1`)

---

## Pipeline Steps

```
*_clean.sam
(Kraken2 human-classified reads)
      │
      ▼ STEP 1 ── SAM → BAM
      │           samtools view + sort + index
      │
      ▼ STEP 2 ── Duplicate Removal (mate-aware)
      │           namesort → fixmate -m → coord sort → markdup -r -s
      │
      ▼ STEP 3 ── BAM → Fragment BED
      │           bamtobed -bedpe → filter chimeric/large fragments
      │           → extract chr/start/end → remove coordinate artifacts
      │
      ▼ STEP 4 ── Chromosome Sizes
      │           extracted directly from BAM header (RefSeq naming)
      │
      ▼ STEP 5 ── bedGraph + CPM Normalization
      │           bedtools genomecov → CPM scale factor applied
      │           (1 count per fragment)
      │
      ▼ STEP 6 ── SEACR Peak Calling
      │           IgG as background control | norm | relaxed
      │
      ▼ STEP 7 ── FRiP Score QC
                  bedtools intersect | reads in peaks / total reads
```

---

## Dependencies

| Tool | Minimum Version | Purpose |
|---|---|---|
| `samtools` | >= 1.13 | BAM processing, duplicate removal |
| `bedtools` | >= 2.30 | Fragment BED, bedGraph, FRiP |
| `SEACR` | 1.3 | Peak calling |
| `R` | >= 3.5 | Required internally by SEACR |
| `bc` | any | CPM and FRiP arithmetic |

### Installing Dependencies

```bash
# samtools
conda install -c bioconda samtools

# bedtools
conda install -c bioconda bedtools

# SEACR (download from GitHub)
git clone https://github.com/FredHutch/SEACR.git
chmod +x SEACR/SEACR_1.3.sh

# R (if not already installed)
conda install -c conda-forge r-base
```

---

## Input Files

| File | Description |
|---|---|
| `*_clean.sam` | Human-classified SAM files from Kraken2 |
| IgG sample | One IgG negative control sample processed alongside others |

### Expected SAM naming convention
```
SAMPLENAME_clean.sam
```

---

## How to Run

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/cutandrun-pipeline.git
cd cutandrun-pipeline
```

### 2. Update the CONFIG section
Open `cutandrun_peak_calling_pipeline.sh` and update the paths
at the top of the script:

```bash
# CONFIG — UPDATE THESE PATHS BEFORE RUNNING
INPUT_DIR="/path/to/kraken/human_classified_sam"
BAM_DIR="/path/to/output/bam"
BEDGRAPH_DIR="/path/to/output/bedgraph"
SEACR_OUT="/path/to/output/seacr"
FRIP_OUT="/path/to/output/qc_frip"
SEACR_SCRIPT="/path/to/SEACR/SEACR_1.3.sh"
IGG_BEDGRAPH="${BEDGRAPH_DIR}/YOUR-IgG-SAMPLE.norm.bedgraph"
```

### 3. Run
```bash
bash cutandrun_peak_calling_pipeline.sh
```

### 4. Monitor progress
Each step prints progress to the terminal:
```
=== STEP 1: SAM to sorted BAM ===
Done: MCIA0354-CD4-IL6-STAT3
...
=== STEP 2: Duplicate Removal ===
...
=== PIPELINE COMPLETE ===
```

---

## Output Files

```
bam/
├── *.sorted.bam                  # Coordinate-sorted BAMs (Step 1)
├── *.namesort.bam                # Intermediate — name sorted (Step 2)
├── *.fixmate.bam                 # Intermediate — fixmate output (Step 2)
├── *.fixmate.sorted.bam          # Intermediate — coord sorted (Step 2)
└── dedup/
    ├── *.rmdup.bam               # Final deduplicated BAMs (Step 2)
    ├── *.rmdup.bam.bai           # BAM index files
    └── *.markdup.stats.txt       # Duplicate removal statistics

bedgraph/
├── *.fragments.clean.bed         # Clean fragment BED files (Step 3)
├── fragment_counts.txt           # Fragment counts per sample (Step 5)
├── hg38.refseq.chrom.sizes       # Chromosome sizes from BAM header (Step 4)
├── *.raw.bedgraph                # Raw fragment coverage (Step 5)
└── *.norm.bedgraph               # CPM-normalized coverage (Step 5)

seacr/
├── *.relaxed.bed                 # SEACR peaks — relaxed mode (Step 6)
└── *.stringent.bed               # SEACR peaks — stringent mode (Step 6)

qc_frip/
└── frip_scores.txt               # FRiP scores for all samples (Step 7)
```

---

## Quality Control

### Duplicate Statistics (from `*.markdup.stats.txt`)

```
READ                    = total reads entering duplicate removal
WRITTEN                 = reads kept after duplicate removal
DUPLICATE PRIMARY TOTAL = reads removed as PCR duplicates

Dup Rate = DUPLICATE PRIMARY TOTAL / READ × 100
```

---

### FRiP Score (from `frip_scores.txt`)

```
FRiP = reads overlapping peaks / total mapped reads
```

---

## Important Notes

### IgG produces no peak file — this is correct
If your IgG sample has no `.relaxed.bed` file after SEACR runs, this is
**expected**. It means SEACR found no regions in the IgG sample that
exceeded background — confirming your negative control is clean.
The FRiP script handles this automatically and records FRiP = 0 for IgG.

### Low-depth samples may produce 0 peaks
Samples with fewer than ~1 million fragments after deduplication may
produce 0 peaks from SEACR. This reflects insufficient sequencing depth
rather than a pipeline failure. These samples should be flagged in your
QC report and interpreted with caution.

### RefSeq chromosome naming
This pipeline uses RefSeq chromosome names (e.g. `NC_000001.11`).
The chromosome sizes file is extracted directly from your BAM header.
Do **not** substitute a standard UCSC `hg38.chrom.sizes` file as it
uses different chromosome names and will cause `bedtools genomecov`
to fail silently.

### Fragment vs read counting
Each CUT&RUN fragment is sequenced as two reads (R1 + R2). This pipeline
counts each fragment **once** using `bedtools bamtobed -bedpe`, which
extracts the full insert coordinates. This ensures signal accurately
reflects the number of DNA molecules captured, not the number of reads.

---

## References

- **CUT&RUN protocol**: Skene PJ & Henikoff S. (2017).
  An efficient targeted nuclease strategy for high-resolution mapping
  of DNA binding sites. *eLife* 6:e21856.

- **SEACR**: Meers MP, Tenenbaum D, Henikoff S. (2019).
  Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin
  profiling. *Epigenetics and Chromatin* 12(1):42.

- **SAMtools**: Danecek P et al. (2021).
  Twelve years of SAMtools and BCFtools.
  *GigaScience* 10(2):giab008.

- **BEDTools**: Quinlan AR & Hall IM. (2010).
  BEDTools: a flexible suite of utilities for comparing genomic features.
  *Bioinformatics* 26(6):841–842.

- **E. coli contamination in CUT&RUN**: Kaya-Okur HS et al. (2019).
  CUT&Tag for efficient epigenomic profiling of small samples and
  single cells. *Nature Communications* 10:1930.

---


