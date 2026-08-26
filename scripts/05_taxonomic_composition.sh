#!/usr/bin/env bash

# ============================================================
# PHASE 5 — LONGITUDINAL TAXONOMIC COMPOSITION
# ============================================================
#
# Purpose:
#   Summarize taxonomic composition across the longitudinal
#   FMT time course using the unrarefied ASV feature table.
#
# Inputs:
#   results/phase2/table.qza
#   results/phase3/taxonomy.qza
#   data/metadata/samples.tsv
#
# Outputs:
#   Taxonomy-collapsed tables at:
#     Level 2 = phylum
#     Level 5 = family
#     Level 6 = genus
#
#   Relative-frequency tables
#   Exported TSV files
#   Taxonomic bar plot visualization
#
# Important:
#   Relative-abundance profiling uses the unrarefied table.
#   Rarefaction was used for diversity analysis in Phase 4,
#   but is not required for descriptive composition summaries.
#
# ============================================================

set -euo pipefail

RESULTS="results/phase5"
LOG_DIR="logs"

TABLE="results/phase2/table.qza"
TAXONOMY="results/phase3/taxonomy.qza"
METADATA="data/metadata/samples.tsv"

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/05_taxonomic_composition.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 5 — LONGITUDINAL TAXONOMIC COMPOSITION"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

# ------------------------------------------------------------
# 1. Taxonomic bar plot
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "TAXONOMIC BAR PLOT"
echo "------------------------------------------------------------"

qiime taxa barplot \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$RESULTS/taxa-bar-plots.qzv"

# ------------------------------------------------------------
# 2. Collapse taxonomy by rank
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "COLLAPSING TAXONOMY"
echo "------------------------------------------------------------"

# SILVA-style taxonomy levels:
# 1 domain
# 2 phylum
# 3 class
# 4 order
# 5 family
# 6 genus
# 7 species

qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 2 \
  --o-collapsed-table "$RESULTS/phylum-table.qza"

qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 5 \
  --o-collapsed-table "$RESULTS/family-table.qza"

qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 6 \
  --o-collapsed-table "$RESULTS/genus-table.qza"

# ------------------------------------------------------------
# 3. Convert to relative abundance
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "RELATIVE ABUNDANCE"
echo "------------------------------------------------------------"

qiime feature-table relative-frequency \
  --i-table "$RESULTS/phylum-table.qza" \
  --o-relative-frequency-table "$RESULTS/phylum-relative.qza"

qiime feature-table relative-frequency \
  --i-table "$RESULTS/family-table.qza" \
  --o-relative-frequency-table "$RESULTS/family-relative.qza"

qiime feature-table relative-frequency \
  --i-table "$RESULTS/genus-table.qza" \
  --o-relative-frequency-table "$RESULTS/genus-relative.qza"

# ------------------------------------------------------------
# 4. Export relative-abundance tables
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING RELATIVE-ABUNDANCE TABLES"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/export"
mkdir -p "$RESULTS/export"

for rank in phylum family genus
do

  qiime tools export \
    --input-path "$RESULTS/${rank}-relative.qza" \
    --output-path "$RESULTS/export/${rank}"

done

# ------------------------------------------------------------
# 5. Convert BIOM exports to TSV
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CONVERTING BIOM TO TSV"
echo "------------------------------------------------------------"

for rank in phylum family genus
do

  biom convert \
    -i "$RESULTS/export/${rank}/feature-table.biom" \
    -o "$RESULTS/export/${rank}/${rank}-relative.tsv" \
    --to-tsv

done

# ------------------------------------------------------------
# 6. Basic summaries
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PHYLUM RELATIVE ABUNDANCE"
echo "------------------------------------------------------------"

cat "$RESULTS/export/phylum/phylum-relative.tsv"

echo
echo "------------------------------------------------------------"
echo "TOP GENUS-LEVEL FEATURES"
echo "------------------------------------------------------------"

head -n 25 "$RESULTS/export/genus/genus-relative.tsv"

# ------------------------------------------------------------
# 7. Conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 5 CONCLUSION"
echo "============================================================"

echo
echo "ASV counts were collapsed to phylum, family, and genus"
echo "levels using SILVA 138.2 taxonomic assignments."
echo
echo "Collapsed count tables were converted to relative"
echo "frequencies using the unrarefied ASV feature table."
echo
echo "Taxonomic bar plots and exported relative-abundance TSV"
echo "files were generated for longitudinal interpretation."
echo
echo "These results will be used to identify bacterial groups"
echo "that increase, decrease, appear, or disappear across"
echo "Before FMT -> 7d -> 14d -> 30d."
echo
echo "============================================================"
echo "END PHASE 5"
echo "============================================================"
