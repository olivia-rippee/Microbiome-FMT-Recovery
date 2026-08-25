#!/usr/bin/env bash

# ============================================================
# PHASE 4 — ALPHA AND BETA DIVERSITY ANALYSIS
# ============================================================
#
# Purpose:
#   Quantify within-sample and between-sample microbiome
#   diversity across the longitudinal FMT time course.
#
# Inputs:
#   results/phase2/table.qza
#   results/phase2/rep-seqs.qza
#   data/metadata/samples.tsv
#
# Outputs:
#   Alpha diversity metrics
#   Beta diversity distance matrices
#   PCoA results
#   Emperor visualizations
#   Phylogenetic tree
#
# Important:
#   This pilot includes one participant with one sample per
#   timepoint. Results are descriptive and should not be
#   interpreted as population-level statistical inference.
#
#   All downstream diversity metrics are calculated from the
#   same table rarefied to 15,000 reads/sample for consistency.
#
# ============================================================

set -euo pipefail

RESULTS="results/phase4"
LOG_DIR="logs"

TABLE="results/phase2/table.qza"
REP_SEQS="results/phase2/rep-seqs.qza"
METADATA="data/metadata/samples.tsv"

SAMPLING_DEPTH=15000

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/04_diversity_analysis.txt"

exec > >(tee "$LOG_FILE") 2>&1

echo "============================================================"
echo "PHASE 4 — ALPHA AND BETA DIVERSITY ANALYSIS"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

echo
echo "Rarefaction depth:"
echo "$SAMPLING_DEPTH"

# ------------------------------------------------------------
# 1. Build phylogenetic tree
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PHYLOGENETIC TREE"
echo "------------------------------------------------------------"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$REP_SEQS" \
  --o-alignment "$RESULTS/aligned-rep-seqs.qza" \
  --o-masked-alignment "$RESULTS/masked-aligned-rep-seqs.qza" \
  --o-tree "$RESULTS/unrooted-tree.qza" \
  --o-rooted-tree "$RESULTS/rooted-tree.qza"

# ------------------------------------------------------------
# 2. Core phylogenetic diversity metrics
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CORE DIVERSITY METRICS"
echo "------------------------------------------------------------"

# Remove previous core-metrics directory so reruns work cleanly
rm -rf "$RESULTS/core-metrics"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --i-table "$TABLE" \
  --p-sampling-depth "$SAMPLING_DEPTH" \
  --m-metadata-file "$METADATA" \
  --output-dir "$RESULTS/core-metrics"

RAREFIED_TABLE="$RESULTS/core-metrics/rarefied_table.qza"

# ------------------------------------------------------------
# 3. Alpha diversity on rarefied table
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "ALPHA DIVERSITY"
echo "------------------------------------------------------------"

qiime diversity alpha \
  --i-table "$RAREFIED_TABLE" \
  --p-metric observed_features \
  --o-alpha-diversity "$RESULTS/observed-features.qza"

qiime diversity alpha \
  --i-table "$RAREFIED_TABLE" \
  --p-metric shannon \
  --o-alpha-diversity "$RESULTS/shannon.qza"

qiime diversity alpha-phylogenetic \
  --i-table "$RAREFIED_TABLE" \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --p-metric faith_pd \
  --o-alpha-diversity "$RESULTS/faith-pd.qza"

# ------------------------------------------------------------
# 4. Export alpha diversity
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING ALPHA DIVERSITY"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/alpha-export"

mkdir -p "$RESULTS/alpha-export"

qiime tools export \
  --input-path "$RESULTS/observed-features.qza" \
  --output-path "$RESULTS/alpha-export/observed-features"

qiime tools export \
  --input-path "$RESULTS/shannon.qza" \
  --output-path "$RESULTS/alpha-export/shannon"

qiime tools export \
  --input-path "$RESULTS/faith-pd.qza" \
  --output-path "$RESULTS/alpha-export/faith-pd"

# ------------------------------------------------------------
# 5. Alpha rarefaction
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "ALPHA RAREFACTION"
echo "------------------------------------------------------------"

qiime diversity alpha-rarefaction \
  --i-table "$TABLE" \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --p-max-depth "$SAMPLING_DEPTH" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$RESULTS/alpha-rarefaction.qzv"

# ------------------------------------------------------------
# 6. Beta diversity on rarefied table
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "BETA DIVERSITY"
echo "------------------------------------------------------------"

qiime diversity beta \
  --i-table "$RAREFIED_TABLE" \
  --p-metric braycurtis \
  --o-distance-matrix "$RESULTS/bray-curtis.qza"

qiime diversity beta \
  --i-table "$RAREFIED_TABLE" \
  --p-metric jaccard \
  --o-distance-matrix "$RESULTS/jaccard.qza"

qiime diversity beta-phylogenetic \
  --i-table "$RAREFIED_TABLE" \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --p-metric weighted_unifrac \
  --o-distance-matrix "$RESULTS/weighted-unifrac.qza"

qiime diversity beta-phylogenetic \
  --i-table "$RAREFIED_TABLE" \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --p-metric unweighted_unifrac \
  --o-distance-matrix "$RESULTS/unweighted-unifrac.qza"

# ------------------------------------------------------------
# 7. PCoA
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PCOA"
echo "------------------------------------------------------------"

for metric in bray-curtis jaccard weighted-unifrac unweighted-unifrac
do

  qiime diversity pcoa \
    --i-distance-matrix "$RESULTS/${metric}.qza" \
    --o-pcoa "$RESULTS/${metric}-pcoa.qza"

  qiime emperor plot \
    --i-pcoa "$RESULTS/${metric}-pcoa.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$RESULTS/${metric}-emperor.qzv"

done

# ------------------------------------------------------------
# 8. Export beta diversity matrices
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING BETA DIVERSITY MATRICES"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/beta-export"

mkdir -p "$RESULTS/beta-export"

for metric in bray-curtis jaccard weighted-unifrac unweighted-unifrac
do

  qiime tools export \
    --input-path "$RESULTS/${metric}.qza" \
    --output-path "$RESULTS/beta-export/${metric}"

done

# ------------------------------------------------------------
# 9. Print alpha diversity results
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "ALPHA DIVERSITY RESULTS"
echo "------------------------------------------------------------"

echo
echo "Observed features:"
cat "$RESULTS/alpha-export/observed-features/alpha-diversity.tsv"

echo
echo "Shannon diversity:"
cat "$RESULTS/alpha-export/shannon/alpha-diversity.tsv"

echo
echo "Faith's PD:"
cat "$RESULTS/alpha-export/faith-pd/alpha-diversity.tsv"

# ------------------------------------------------------------
# 10. Print selected beta diversity results
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "SELECTED BETA DIVERSITY RESULTS"
echo "------------------------------------------------------------"

echo
echo "Bray-Curtis:"
cat "$RESULTS/beta-export/bray-curtis/distance-matrix.tsv"

echo
echo "Weighted UniFrac:"
cat "$RESULTS/beta-export/weighted-unifrac/distance-matrix.tsv"

# ------------------------------------------------------------
# 11. Final conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 4 CONCLUSION"
echo "============================================================"

echo
echo "Alpha and beta diversity metrics were calculated across"
echo "the four longitudinal FMT samples."
echo
echo "All downstream diversity metrics were calculated using the"
echo "same feature table rarefied to 15,000 reads per sample."
echo
echo "Alpha diversity included observed features, Shannon"
echo "diversity, and Faith's phylogenetic diversity."
echo
echo "Beta diversity included Bray-Curtis, Jaccard, weighted"
echo "UniFrac, and unweighted UniFrac distances."
echo
echo "PCoA and Emperor visualizations were generated for each"
echo "beta-diversity metric."
echo
echo "A rarefaction depth of 15,000 reads was selected so that"
echo "all four samples were retained in the pilot analysis."
echo
echo "Because only one participant is included, diversity"
echo "patterns are interpreted descriptively rather than as"
echo "population-level statistical inference."
echo
echo "============================================================"
echo "END PHASE 4"
echo "============================================================"
