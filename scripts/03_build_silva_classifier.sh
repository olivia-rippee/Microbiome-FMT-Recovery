#!/usr/bin/env bash

# ============================================================
# PHASE 3A — BUILD SILVA V4 TAXONOMIC CLASSIFIER
# ============================================================
#
# Purpose:
#   Build a SILVA 138.2 SSURef NR99 V4-region classifier
#   compatible with the current QIIME 2 / Rachis 2026.7
#   environment.
#
# Rationale:
#   A downloaded pre-trained SILVA classifier was incompatible
#   with the installed scikit-learn version.
#
#   Full-length classifier training exceeded available memory,
#   so a V4-region-specific reference set was extracted,
#   dereplicated, and used to train a smaller Naive Bayes
#   classifier locally.
#
# ============================================================

set -euo pipefail

REFERENCE="reference/silva138.2"
LOG_DIR="logs"

mkdir -p "$REFERENCE"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/03_build_silva_classifier.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 3A — BUILD SILVA V4 CLASSIFIER"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

# ------------------------------------------------------------
# 1. Retrieve SILVA reference data
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "RETRIEVING SILVA 138.2 SSURef NR99"
echo "------------------------------------------------------------"

if [ ! -f "$REFERENCE/silva-138.2-ssu-nr99-rna-seqs.qza" ]; then

    qiime rescript get-silva-data \
      --p-version '138.2' \
      --p-target 'SSURef_NR99' \
      --o-silva-sequences "$REFERENCE/silva-138.2-ssu-nr99-rna-seqs.qza" \
      --o-silva-taxonomy "$REFERENCE/silva-138.2-ssu-nr99-taxonomy.qza"

else
    echo "SILVA RNA reference already exists; skipping download."
fi

# ------------------------------------------------------------
# 2. Convert RNA sequences to DNA
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CONVERTING RNA REFERENCE TO DNA"
echo "------------------------------------------------------------"

if [ ! -f "$REFERENCE/silva-138.2-ssu-nr99-dna-seqs.qza" ]; then

    qiime rescript reverse-transcribe \
      --i-rna-sequences "$REFERENCE/silva-138.2-ssu-nr99-rna-seqs.qza" \
      --o-dna-sequences "$REFERENCE/silva-138.2-ssu-nr99-dna-seqs.qza"

else
    echo "DNA reference already exists; skipping conversion."
fi

# ------------------------------------------------------------
# 3. Extract V4 reference region
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXTRACTING V4 REGION"
echo "------------------------------------------------------------"

if [ ! -f "$REFERENCE/silva-138.2-v4-reads.qza" ]; then

    qiime feature-classifier extract-reads \
      --i-sequences "$REFERENCE/silva-138.2-ssu-nr99-dna-seqs.qza" \
      --p-f-primer GTGYCAGCMGCCGCGGTAA \
      --p-r-primer GGACTACNVGGGTWTCTAAT \
      --p-min-length 200 \
      --p-max-length 400 \
      --o-reads "$REFERENCE/silva-138.2-v4-reads.qza" \
      --o-read-extraction-stats "$REFERENCE/silva-138.2-v4-extraction-stats.qza"

else
    echo "V4 reference reads already exist; skipping extraction."
fi

# ------------------------------------------------------------
# 4. Dereplicate V4 reference
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "DEREPLICATING V4 REFERENCE"
echo "------------------------------------------------------------"

if [ ! -f "$REFERENCE/silva-138.2-v4-uniq-seqs.qza" ]; then

    qiime rescript dereplicate \
      --i-sequences "$REFERENCE/silva-138.2-v4-reads.qza" \
      --i-taxa "$REFERENCE/silva-138.2-ssu-nr99-taxonomy.qza" \
      --p-mode uniq \
      --o-dereplicated-sequences "$REFERENCE/silva-138.2-v4-uniq-seqs.qza" \
      --o-dereplicated-taxa "$REFERENCE/silva-138.2-v4-uniq-taxonomy.qza"

else
    echo "Dereplicated V4 reference already exists; skipping."
fi

# ------------------------------------------------------------
# 5. Train classifier
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "TRAINING NAIVE BAYES CLASSIFIER"
echo "------------------------------------------------------------"

if [ ! -f "$REFERENCE/silva-138.2-v4-classifier-qiime2026.7.qza" ]; then

    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "$REFERENCE/silva-138.2-v4-uniq-seqs.qza" \
      --i-reference-taxonomy "$REFERENCE/silva-138.2-v4-uniq-taxonomy.qza" \
      --o-classifier "$REFERENCE/silva-138.2-v4-classifier-qiime2026.7.qza"

else
    echo "Classifier already exists; skipping training."
fi

# ------------------------------------------------------------
# 6. Verify classifier
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CLASSIFIER VERIFICATION"
echo "------------------------------------------------------------"

qiime tools peek \
  "$REFERENCE/silva-138.2-v4-classifier-qiime2026.7.qza"

echo
echo "============================================================"
echo "PHASE 3A CONCLUSION"
echo "============================================================"

echo
echo "A SILVA 138.2 SSURef NR99 V4-region reference set was"
echo "generated using in-silico 515F/806R-compatible extraction."
echo
echo "The extracted sequences were dereplicated and used to train"
echo "a Naive Bayes classifier in the local QIIME 2 / Rachis"
echo "2026.7 environment."
echo
echo "This avoided scikit-learn incompatibility with externally"
echo "trained classifiers and reduced memory requirements relative"
echo "to training against the full-length SILVA reference."
echo
echo "============================================================"
echo "END PHASE 3A"
echo "============================================================"
