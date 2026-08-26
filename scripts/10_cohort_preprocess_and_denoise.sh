#!/usr/bin/env bash

# ============================================================
# PHASE 10 — COHORT PREPROCESSING AND ASV INFERENCE
# ============================================================
#
# Purpose:
#   Denoise the complete longitudinal FMT cohort using DADA2
#   and generate cohort-wide ASV and read-retention outputs.
#
# Cohort:
#   52 participants
#   208 samples
#
# Inputs:
#   results/phase10/demux-paired.qza
#   data/cohort/samples.tsv
#
# DADA2 parameters:
#   Forward truncation: 230 bp
#   Reverse truncation: 190 bp
#   Maximum expected errors: 2 / 2
#   Minimum overlap: 12 bp
#
# Primer trimming:
#   Skipped, consistent with the pilot analysis. Primer motifs
#   were previously detected in only a small fraction of reads,
#   indicating upstream primer removal.
#
# ============================================================

set -euo pipefail

RESULTS="results/phase10"
LOG_DIR="logs"

DEMUX="$RESULTS/demux-paired.qza"
METADATA="data/cohort/samples.tsv"

TRUNC_LEN_F=230
TRUNC_LEN_R=190
MAX_EE_F=2
MAX_EE_R=2
MIN_OVERLAP=12

THREADS=4

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/10_cohort_preprocess_and_denoise.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 10 — COHORT PREPROCESSING AND ASV INFERENCE"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

echo
echo "Samples:"
tail -n +2 "$METADATA" | wc -l

# ------------------------------------------------------------
# 1. Verify imported sequences
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "INPUT VERIFICATION"
echo "------------------------------------------------------------"

qiime tools peek "$DEMUX"

# ------------------------------------------------------------
# 2. Processing parameters
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "DADA2 PARAMETERS"
echo "------------------------------------------------------------"

echo "Forward truncation length: $TRUNC_LEN_F"
echo "Reverse truncation length: $TRUNC_LEN_R"
echo "Maximum expected errors R1: $MAX_EE_F"
echo "Maximum expected errors R2: $MAX_EE_R"
echo "Minimum overlap: $MIN_OVERLAP"
echo "Threads: $THREADS"

echo
echo "Primer trimming:"
echo "Skipped. Primer sequences were previously determined to be"
echo "largely absent from the source reads."

# ------------------------------------------------------------
# 3. DADA2
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "DADA2 DENOISING"
echo "------------------------------------------------------------"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$DEMUX" \
  --p-trunc-len-f "$TRUNC_LEN_F" \
  --p-trunc-len-r "$TRUNC_LEN_R" \
  --p-max-ee-f "$MAX_EE_F" \
  --p-max-ee-r "$MAX_EE_R" \
  --p-min-overlap "$MIN_OVERLAP" \
  --p-n-threads "$THREADS" \
  --o-table "$RESULTS/table.qza" \
  --o-representative-sequences "$RESULTS/rep-seqs.qza" \
  --o-denoising-stats "$RESULTS/denoising-stats.qza" \
  --o-base-transition-stats "$RESULTS/base-transition-stats.qza"

# ------------------------------------------------------------
# 4. Feature table summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "FEATURE TABLE SUMMARY"
echo "------------------------------------------------------------"

qiime feature-table summarize \
  --i-table "$RESULTS/table.qza" \
  --m-metadata-file "$METADATA" \
  --o-feature-frequencies "$RESULTS/feature-frequencies.qza" \
  --o-sample-frequencies "$RESULTS/sample-frequencies.qza" \
  --o-summary "$RESULTS/table.qzv"

# ------------------------------------------------------------
# 5. Representative sequence summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "REPRESENTATIVE SEQUENCES"
echo "------------------------------------------------------------"

qiime feature-table tabulate-seqs \
  --i-data "$RESULTS/rep-seqs.qza" \
  --o-visualization "$RESULTS/rep-seqs.qzv"

# ------------------------------------------------------------
# 6. DADA2 statistics
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "DENOISING STATISTICS"
echo "------------------------------------------------------------"

qiime metadata tabulate \
  --m-input-file "$RESULTS/denoising-stats.qza" \
  --o-visualization "$RESULTS/denoising-stats.qzv"

# ------------------------------------------------------------
# 7. Export compact results
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING QC TABLES"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/denoising-stats-export"
rm -rf "$RESULTS/sample-frequencies-export"
rm -rf "$RESULTS/feature-frequencies-export"

qiime tools export \
  --input-path "$RESULTS/denoising-stats.qza" \
  --output-path "$RESULTS/denoising-stats-export"

qiime tools export \
  --input-path "$RESULTS/sample-frequencies.qza" \
  --output-path "$RESULTS/sample-frequencies-export"

qiime tools export \
  --input-path "$RESULTS/feature-frequencies.qza" \
  --output-path "$RESULTS/feature-frequencies-export"

# ------------------------------------------------------------
# 8. Cohort retention summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "COHORT READ RETENTION"
echo "------------------------------------------------------------"

python - <<'PY'
import pandas as pd

f = "results/phase10/denoising-stats-export/stats.tsv"

df = pd.read_csv(
    f,
    sep="\t",
    comment="#"
)

col = "percentage of input non-chimeric"

print(f"Samples: {len(df)}")
print()

print("Non-chimeric read retention (%):")
print(f"Minimum: {df[col].min():.2f}")
print(f"25th percentile: {df[col].quantile(0.25):.2f}")
print(f"Median: {df[col].median():.2f}")
print(f"75th percentile: {df[col].quantile(0.75):.2f}")
print(f"Maximum: {df[col].max():.2f}")

print()
print("Samples below 20% retention:")
low = df[df[col] < 20]

if len(low):
    print(
        low[
            ["sample-id", "input", "non-chimeric", col]
        ].to_string(index=False)
    )
else:
    print("None")
PY

# ------------------------------------------------------------
# 9. Final outputs
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PHASE 10 OUTPUTS"
echo "------------------------------------------------------------"

find "$RESULTS" \
  -maxdepth 2 \
  -type f \
  | sort

# ------------------------------------------------------------
# 10. Conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 10 CONCLUSION"
echo "============================================================"

echo
echo "The complete longitudinal recipient cohort was processed"
echo "using paired-end DADA2 denoising."
echo
echo "Forward and reverse reads were truncated to 230 and 190 bp,"
echo "respectively, based on cohort-wide quality profiles."
echo
echo "DADA2 performed quality filtering, denoising, paired-read"
echo "merging, and chimera removal."
echo
echo "The resulting cohort-wide ASV feature table and"
echo "representative sequences were generated for downstream"
echo "taxonomy, diversity, longitudinal, and differential"
echo "abundance analyses."
echo
echo "============================================================"
echo "END PHASE 10"
echo "============================================================"
