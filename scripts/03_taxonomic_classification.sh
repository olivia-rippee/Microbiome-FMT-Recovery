#!/usr/bin/env bash

# ============================================================
# PHASE 3 — TAXONOMIC CLASSIFICATION
# ============================================================
#
# Purpose:
#   Assign taxonomy to the DADA2-derived ASVs using a
#   locally trained SILVA V4 classifier and generate compact,
#   reproducible taxonomy outputs for downstream analysis.
#
# Input:
#   results/phase2/rep-seqs.qza
#   results/phase2/table.qza
#
# Reference:
#   Custom SILVA 138.2 SSURef NR99 V4-region classifier
#   trained in QIIME 2 / Rachis 2026.7
#
# V4 reference reads were extracted in silico using
# 515F/806R-compatible primers, dereplicated, and used
# to train a Naive Bayes classifier in the current
# scikit-learn/QIIME 2 environment.
#
# Output:
#   results/phase3/taxonomy.qza
#   results/phase3/taxonomy.qzv
#   results/phase3/taxa-bar-plots.qzv
#   results/phase3/taxonomy-export/taxonomy.tsv
#
# ============================================================

set -euo pipefail

RESULTS="results/phase3"
LOG_DIR="logs"

CLASSIFIER="reference/silva138.2/silva-138.2-v4-classifier-qiime2026.7.qza"
REP_SEQS="results/phase2/rep-seqs.qza"
TABLE="results/phase2/table.qza"
METADATA="data/metadata/samples.tsv"

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/03_taxonomic_classification.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 3 — TAXONOMIC CLASSIFICATION"
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
# 1. Taxonomic classification
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CLASSIFYING ASVs"
echo "------------------------------------------------------------"

qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$REP_SEQS" \
    --p-n-jobs 1 \
    --o-classification "$RESULTS/taxonomy.qza"

# ------------------------------------------------------------
# 2. Taxonomy table
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

mkdir -p "$RESULTS/taxonomy-export"

qiime tools export \
    --input-path "$RESULTS/taxonomy.qza" \
    --output-path "$RESULTS/taxonomy-export"

echo
echo "Exported taxonomy:"
ls -lh "$RESULTS/taxonomy-export"

# ------------------------------------------------------------
# 5. Basic summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "TAXONOMY SUMMARY"
echo "------------------------------------------------------------"

if [ -f "$RESULTS/taxonomy-export/taxonomy.tsv" ]; then

    echo "Number of classified ASVs:"
    awk 'NR > 1 {count++} END {print count}' \
        "$RESULTS/taxonomy-export/taxonomy.tsv"

    echo
    echo "First classified features:"
    head -n 10 "$RESULTS/taxonomy-export/taxonomy.tsv"

fi

# ------------------------------------------------------------
# Conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 3 CONCLUSION"
echo "============================================================"

echo
echo "DADA2-derived ASVs were classified using a custom"
echo "SILVA 138.2 SSURef NR99 V4-region Naive Bayes classifier"
echo "trained locally in the QIIME 2 / Rachis 2026.7 environment."
echo
echo "Taxonomic assignments were saved as a QIIME 2 artifact,"
echo "tabular visualization, exported TSV file, and longitudinal"
echo "taxonomic bar-plot visualization."
echo
echo "Taxonomic results will be used in subsequent phases to"
echo "characterize microbial community composition and changes"
echo "across the FMT time course."
echo
echo "============================================================"
echo "END PHASE 3"
echo "============================================================"
