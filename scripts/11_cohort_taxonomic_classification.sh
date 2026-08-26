#!/usr/bin/env bash

# ============================================================
# PHASE 11 — COHORT TAXONOMIC CLASSIFICATION
# ============================================================
#
# Purpose:
#   Assign taxonomy to cohort-wide DADA2 ASVs using the
#   locally trained SILVA 138.2 V4 classifier.
#
# Inputs:
#   results/phase10/rep-seqs.qza
#   results/phase10/table.qza
#   data/cohort/samples.tsv
#
# Reference:
#   reference/silva138.2/
#   silva-138.2-v4-classifier-qiime2026.7.qza
#
# Outputs:
#   results/phase11/taxonomy.qza
#   results/phase11/taxonomy.qzv
#   results/phase11/taxa-bar-plots.qzv
#   results/phase11/taxonomy-export/taxonomy.tsv
#
# ============================================================

set -euo pipefail

RESULTS="results/phase11"
LOG_DIR="logs"

CLASSIFIER="reference/silva138.2/silva-138.2-v4-classifier-qiime2026.7.qza"
REP_SEQS="results/phase10/rep-seqs.qza"
TABLE="results/phase10/table.qza"
METADATA="data/cohort/samples.tsv"

THREADS=4

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/11_cohort_taxonomic_classification.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 11 — COHORT TAXONOMIC CLASSIFICATION"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

echo
echo "Classifier:"
qiime tools peek "$CLASSIFIER"

echo
echo "Representative sequences:"
qiime tools peek "$REP_SEQS"

# ------------------------------------------------------------
# 1. Classify ASVs
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CLASSIFYING ASVs"
echo "------------------------------------------------------------"

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$REP_SEQS" \
  --p-n-jobs "$THREADS" \
  --o-classification "$RESULTS/taxonomy.qza"

# ------------------------------------------------------------
# 2. Taxonomy visualization
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "GENERATING TAXONOMY TABLE"
echo "------------------------------------------------------------"

qiime metadata tabulate \
  --m-input-file "$RESULTS/taxonomy.qza" \
  --o-visualization "$RESULTS/taxonomy.qzv"

# ------------------------------------------------------------
# 3. Taxonomic bar plots
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "GENERATING TAXONOMIC BAR PLOTS"
echo "------------------------------------------------------------"

qiime taxa barplot \
  --i-table "$TABLE" \
  --i-taxonomy "$RESULTS/taxonomy.qza" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$RESULTS/taxa-bar-plots.qzv"

# ------------------------------------------------------------
# 4. Export taxonomy
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING TAXONOMY"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/taxonomy-export"

qiime tools export \
  --input-path "$RESULTS/taxonomy.qza" \
  --output-path "$RESULTS/taxonomy-export"

# ------------------------------------------------------------
# 5. Taxonomy summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "TAXONOMY SUMMARY"
echo "------------------------------------------------------------"

python - <<'PY'
import pandas as pd

f = "results/phase11/taxonomy-export/taxonomy.tsv"

df = pd.read_csv(
    f,
    sep="\t"
)

print(f"Classified ASVs: {len(df)}")

print()
print("Confidence summary:")

if "Confidence" in df.columns:
    print(df["Confidence"].describe().to_string())

print()
print("First classified ASVs:")

print(
    df.head(10).to_string(index=False)
)
PY

# ------------------------------------------------------------
# 6. Final outputs
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PHASE 11 OUTPUTS"
echo "------------------------------------------------------------"

find "$RESULTS" \
  -maxdepth 2 \
  -type f \
  | sort

# ------------------------------------------------------------
# 7. Conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 11 CONCLUSION"
echo "============================================================"

echo
echo "Cohort-wide DADA2 ASVs were classified using the locally"
echo "trained SILVA 138.2 V4 Naive Bayes classifier."
echo
echo "Taxonomic assignments were saved as a QIIME 2 artifact,"
echo "tabular visualization, exported TSV file, and taxonomic"
echo "bar-plot visualization."
echo
echo "These taxonomy results will support downstream cohort-level"
echo "composition analysis, diversity interpretation, and"
echo "differential-abundance testing."
echo
echo "============================================================"
echo "END PHASE 11"
echo "============================================================"
