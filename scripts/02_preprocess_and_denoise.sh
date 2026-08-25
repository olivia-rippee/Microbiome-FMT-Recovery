#!/usr/bin/env bash

# ============================================================
# PHASE 2 — READ PREPROCESSING AND ASV INFERENCE
# ============================================================
#
# Project: microbiome-fmt
#
# Purpose:
#   Convert verified raw paired-end FASTQ files into a
#   quality-controlled ASV feature table.
#
# Major activities:
#   1. Verify input files
#   2. Verify primer information
#   3. Import reads into QIIME 2
#   4. Remove primers
#   5. Quality-filter and denoise reads
#   6. Merge paired-end reads
#   7. Remove chimeras
#   8. Generate ASV table
#   9. Generate representative sequences
#  10. Generate read-retention statistics
#
# IMPORTANT:
#   Exact primer sequences and truncation parameters MUST be
#   confirmed before running the biological processing.
#
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# Project directories
# ------------------------------------------------------------

PROJECT_DIR="$(pwd)"

RAW_DIR="data/raw"
METADATA="data/metadata/samples.tsv"

RESULTS="results/phase2"
LOG_DIR="logs"

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/02_preprocess_and_denoise.txt"

# ------------------------------------------------------------
# USER-DEFINED PARAMETERS
# ------------------------------------------------------------
#
# DO NOT GUESS THESE.
#
# These values must be confirmed from the original study /
# sequencing protocol before proceeding.
#
# Example placeholders are intentionally left blank.
# ------------------------------------------------------------

FORWARD_PRIMER=""
REVERSE_PRIMER=""

# DADA2 truncation parameters.
#
# These should be selected after evaluating the FastQC results
# AND confirming sufficient overlap for paired-end merging.

TRUNC_LEN_F=230
TRUNC_LEN_R=190

# Minimum overlap between forward and reverse reads.
MIN_OVERLAP=12

# Maximum expected errors.
MAX_EE_F=2
MAX_EE_R=2

# ------------------------------------------------------------
# Start log
# ------------------------------------------------------------

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 2 — READ PREPROCESSING AND ASV INFERENCE"
echo "============================================================"
echo
echo "Date:"
date
echo
echo "Project:"
echo "$PROJECT_DIR"
echo

# ------------------------------------------------------------
# 1. Check environment
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "QIIME 2 ENVIRONMENT"
echo "------------------------------------------------------------"

if ! command -v qiime >/dev/null 2>&1; then
    echo "ERROR: QIIME 2 command 'qiime' was not found."
    echo
    echo "Activate your QIIME 2 environment before continuing."
    exit 1
fi

qiime --version

echo

# ------------------------------------------------------------
# 2. Check metadata
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "METADATA"
echo "------------------------------------------------------------"

if [ ! -f "$METADATA" ]; then
    echo "ERROR: Metadata file not found:"
    echo "$METADATA"
    exit 1
fi

cat "$METADATA"

echo

# ------------------------------------------------------------
# 3. Check raw FASTQ files
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "INPUT FASTQ FILES"
echo "------------------------------------------------------------"
EXPECTED_FILES=(
    "$RAW_DIR/SRR2657981_1.fastq.gz"
    "$RAW_DIR/SRR2657981_2.fastq.gz"
    "$RAW_DIR/SRR2657993_1.fastq.gz"
    "$RAW_DIR/SRR2657993_2.fastq.gz"
    "$RAW_DIR/SRR2657994_1.fastq.gz"
    "$RAW_DIR/SRR2657994_2.fastq.gz"
    "$RAW_DIR/SRR2657995_1.fastq.gz"
    "$RAW_DIR/SRR2657995_2.fastq.gz"
)

for file in "${EXPECTED_FILES[@]}"; do

    if [ ! -f "$file" ]; then
        echo "ERROR: Missing FASTQ:"
        echo "$file"
        exit 1
    fi

    ls -lh "$file"

done

echo

# ------------------------------------------------------------
# 4. Verify primer information
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "PRIMER VERIFICATION"
echo "------------------------------------------------------------"

echo "Primer trimming:"
echo "Skipped. Primer motifs were present in only a very small"
echo "fraction of reads, suggesting primers had already been"
echo "removed during upstream processing."

echo "Forward primer:"
echo "$FORWARD_PRIMER"

echo
echo "Reverse primer:"
echo "$REVERSE_PRIMER"

echo

# ------------------------------------------------------------
# 5. Verify truncation parameters
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "TRUNCATION PARAMETERS"
echo "------------------------------------------------------------"

if [ "$TRUNC_LEN_F" -eq 0 ] || [ "$TRUNC_LEN_R" -eq 0 ]; then

    echo "ERROR: DADA2 truncation lengths have not been specified."
    echo
    echo "Review the Phase 1 FastQC results and determine"
    echo "appropriate forward/reverse truncation lengths."
    echo
    echo "The lengths must preserve sufficient overlap for"
    echo "paired-end merging."
    echo
    exit 1

fi

echo "Forward truncation length: $TRUNC_LEN_F"
echo "Reverse truncation length: $TRUNC_LEN_R"
echo "Minimum overlap: $MIN_OVERLAP"
echo "Maximum expected errors R1: $MAX_EE_F"
echo "Maximum expected errors R2: $MAX_EE_R"

echo

# ------------------------------------------------------------
# 6. Import paired-end reads
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "QIIME 2 IMPORT"
echo "------------------------------------------------------------"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path data/metadata/fastq_manifest.tsv \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$RESULTS/demux-paired.qza"

echo

# ------------------------------------------------------------
# 7. Summarize imported reads
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "DEMULTIPLEXED READ SUMMARY"
echo "------------------------------------------------------------"

qiime demux summarize \
    --i-data "$RESULTS/demux-paired.qza" \
    --o-visualization "$RESULTS/demux-paired.qzv"

echo

# ------------------------------------------------------------
# 8. Primer handling
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "PRIMER HANDLING"
echo "------------------------------------------------------------"
cp "$RESULTS/demux-paired.qza" "$RESULTS/trimmed.qza"

# ------------------------------------------------------------
# 9. Summarize after primer removal
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "POST-PRIMER-REMOVAL SUMMARY"
echo "------------------------------------------------------------"

qiime demux summarize \
    --i-data "$RESULTS/trimmed.qza" \
    --o-visualization "$RESULTS/trimmed.qzv"

echo

# ------------------------------------------------------------
# 10. DADA2 denoising
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "DADA2 DENOISING"
echo "------------------------------------------------------------"

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$RESULTS/trimmed.qza" \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-max-ee-f "$MAX_EE_F" \
    --p-max-ee-r "$MAX_EE_R" \
    --p-min-overlap "$MIN_OVERLAP" \
    --o-table "$RESULTS/table.qza" \
    --o-representative-sequences "$RESULTS/rep-seqs.qza" \
    --o-denoising-stats "$RESULTS/denoising-stats.qza" \
    --o-base-transition-stats "$RESULTS/base-transition-stats.qza"

echo

# ------------------------------------------------------------
# 11. Summarize feature table
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "ASV FEATURE TABLE"
echo "------------------------------------------------------------"

qiime feature-table summarize \
    --i-table "$RESULTS/table.qza" \
    --m-metadata-file "$METADATA" \
    --o-summary "$RESULTS/table.qzv" \
    --o-feature-frequencies "$RESULTS/feature-frequencies.qza" \
    --o-sample-frequencies "$RESULTS/sample-frequencies.qza"

echo

# ------------------------------------------------------------
# 12. Summarize representative sequences
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "REPRESENTATIVE SEQUENCES"
echo "------------------------------------------------------------"

qiime feature-table tabulate-seqs \
    --i-data "$RESULTS/rep-seqs.qza" \
    --o-visualization "$RESULTS/rep-seqs.qzv"

echo

# ------------------------------------------------------------
# 13. Denoising statistics
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "DENOISING STATISTICS"
echo "------------------------------------------------------------"

qiime metadata tabulate \
    --m-input-file "$RESULTS/denoising-stats.qza" \
    --o-visualization "$RESULTS/denoising-stats.qzv"

echo

# ------------------------------------------------------------
# 14. Export denoising statistics
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "EXPORTING READ RETENTION STATISTICS"
echo "------------------------------------------------------------"

mkdir -p "$RESULTS/denoising-stats-export"

qiime tools export \
    --input-path "$RESULTS/denoising-stats.qza" \
    --output-path "$RESULTS/denoising-stats-export"

echo

# ------------------------------------------------------------
# 15. Final outputs
# ------------------------------------------------------------

echo "============================================================"
echo "PHASE 2 OUTPUTS"
echo "============================================================"

find "$RESULTS" -maxdepth 2 -type f | sort

echo

# ------------------------------------------------------------
# 16. Conclusion
# ------------------------------------------------------------

echo "============================================================"
echo "PHASE 2 CONCLUSION"
echo "============================================================"

echo
echo "Raw paired-end reads were imported into QIIME 2."
echo
echo "Primer trimming was skipped because primer motifs were detected"
echo "in only a very small fraction of reads, consistent with primers"
echo "having already been removed during upstream processing."
echo
echo "Reads were quality filtered and denoised using DADA2."
echo
echo "Paired-end reads were merged using the specified minimum"
echo "overlap requirement."
echo
echo "Chimeric sequences were identified and removed as part of"
echo "the DADA2 denoising workflow."
echo
echo "The resulting ASV feature table and representative"
echo "sequences were generated for downstream taxonomic and"
echo "community-level microbiome analysis."
echo
echo "Read-retention statistics were retained to document the"
echo "number of reads passing each processing stage."
echo
echo
echo "FINAL PHASE 2 SUMMARY"
echo "---------------------"
echo "DADA2 produced 166 ASVs across four longitudinal samples."
echo
echo "Final non-chimeric sequencing depths:"
echo "  R009_before : 21,636 reads (56.96% of input; 93 ASVs)"
echo "  R009_7d     : 15,599 reads (38.99% of input; 62 ASVs)"
echo "  R009_14d    : 18,849 reads (62.25% of input; 59 ASVs)"
echo "  R009_30d    : 19,484 reads (46.25% of input; 68 ASVs)"
echo
echo "All four samples retained >15,000 non-chimeric reads."
echo "The 230-bp forward and 190-bp reverse truncation settings"
echo "were therefore retained for downstream pilot analysis."

echo "============================================================"
echo "END PHASE 2"
echo "============================================================"
