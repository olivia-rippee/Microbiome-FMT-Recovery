#!/usr/bin/env bash

# ============================================================
# PHASE 13 — LONGITUDINAL COHORT STATISTICS
# ============================================================
#
# Purpose:
#   Perform repeated-measures statistical analysis of
#   cohort-wide alpha and beta diversity following FMT.
#
# Cohort:
#   52 participants
#   208 samples
#   before, 7d, 14d, 30d
#
# Primary analyses:
#   1. Paired alpha-diversity comparisons
#   2. Linear mixed-effects alpha-diversity trajectories
#   3. Beta-diversity distance from participant baseline
#   4. Longitudinal modeling of baseline displacement
#
# ============================================================

set -euo pipefail

RESULTS="results/phase13"
LOG_DIR="logs"

METADATA="data/cohort/samples.tsv"
ALPHA="results/phase12/alpha-diversity-summary.tsv"

mkdir -p "$RESULTS"
mkdir -p "$RESULTS/alpha-pairwise"
mkdir -p "$RESULTS/alpha-lme"
mkdir -p "$RESULTS/beta-first-distances"
mkdir -p "$RESULTS/beta-lme"
mkdir -p "$RESULTS/tables"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/13_longitudinal_statistics.txt"

exec > >(tee "$LOG_FILE") 2>&1


echo "============================================================"
echo "PHASE 13 — LONGITUDINAL COHORT STATISTICS"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version


# ============================================================
# 1. CREATE QIIME METADATA CONTAINING ALPHA DIVERSITY
# ============================================================

echo
echo "------------------------------------------------------------"
echo "PREPARING ALPHA-DIVERSITY METADATA"
echo "------------------------------------------------------------"

python - <<'PY'
import pandas as pd

f = "results/phase12/alpha-diversity-summary.tsv"

df = pd.read_csv(
    f,
    sep="\t"
)

required = [
    "sample-id",
    "participant",
    "timepoint",
    "days_after_fmt",
    "observed_features",
    "shannon",
    "faith_pd",
    "evenness",
]

missing = [
    x for x in required
    if x not in df.columns
]

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )

out = df[required].copy()

out.to_csv(
    "results/phase13/"
    "alpha-metadata.tsv",
    sep="\t",
    index=False
)

print(f"Samples: {len(out)}")
print(
    f"Participants: "
    f"{out['participant'].nunique()}"
)

print()
print("Samples per timepoint:")
print(
    out["timepoint"]
    .value_counts()
    .reindex(
        ["before", "7d", "14d", "30d"]
    )
    .to_string()
)
PY


ALPHA_META="$RESULTS/alpha-metadata.tsv"


# ============================================================
# 2. QIIME PAIRED ALPHA-DIVERSITY TESTS
# ============================================================

echo
echo "------------------------------------------------------------"
echo "PAIRED ALPHA-DIVERSITY TESTS"
echo "------------------------------------------------------------"

METRICS=(
    observed_features
    shannon
    faith_pd
    evenness
)

STATE1=(
    before
    before
    before
    7d
    7d
    14d
)

STATE2=(
    7d
    14d
    30d
    14d
    30d
    30d
)

for METRIC in "${METRICS[@]}"
do

    echo
    echo "Metric: $METRIC"

    for i in "${!STATE1[@]}"
    do

        A="${STATE1[$i]}"
        B="${STATE2[$i]}"

        NAME="${METRIC}_${A}_vs_${B}"

        echo "  $A vs $B"

        qiime longitudinal pairwise-differences \
          --m-metadata-file "$ALPHA_META" \
          --p-metric "$METRIC" \
          --p-state-column timepoint \
          --p-state-1 "$A" \
          --p-state-2 "$B" \
          --p-individual-id-column participant \
          --p-no-parametric \
          --p-replicate-handling error \
          --o-visualization \
          "$RESULTS/alpha-pairwise/${NAME}.qzv"

    done

done


# ============================================================
# 3. LINEAR MIXED-EFFECTS ALPHA MODELS
# ============================================================

echo
echo "------------------------------------------------------------"
echo "ALPHA-DIVERSITY LINEAR MIXED-EFFECTS MODELS"
echo "------------------------------------------------------------"

for METRIC in "${METRICS[@]}"
do

    echo
    echo "Modeling: $METRIC"

    qiime longitudinal linear-mixed-effects \
  	--m-metadata-file "$ALPHA_META" \
  	--p-state-column days_after_fmt \
  	--p-individual-id-column participant \
  	--p-metric "$METRIC" \
  	--o-visualization \
  	"$RESULTS/alpha-lme/${METRIC}_lme.qzv"

done


# ============================================================
# 4. CALCULATE BETA DISTANCE FROM BASELINE WITH QIIME
# ============================================================

echo
echo "------------------------------------------------------------"
echo "BETA-DIVERSITY DISTANCE FROM BASELINE"
echo "------------------------------------------------------------"

declare -A BETA

BETA[bray_curtis]="results/phase12/core-metrics/bray_curtis_distance_matrix.qza"
BETA[jaccard]="results/phase12/core-metrics/jaccard_distance_matrix.qza"
BETA[weighted_unifrac]="results/phase12/core-metrics/weighted_unifrac_distance_matrix.qza"
BETA[unweighted_unifrac]="results/phase12/core-metrics/unweighted_unifrac_distance_matrix.qza"

for METRIC in \
    bray_curtis \
    jaccard \
    weighted_unifrac \
    unweighted_unifrac

do

    echo
    echo "Metric: $METRIC"

    qiime longitudinal first-distances \
      --i-distance-matrix "${BETA[$METRIC]}" \
      --m-metadata-file "$METADATA" \
      --p-state-column days_after_fmt \
      --p-individual-id-column participant \
      --p-baseline 0 \
      --p-replicate-handling error \
      --o-first-distances \
      "$RESULTS/beta-first-distances/${METRIC}.qza"

    rm -rf \
      "$RESULTS/beta-first-distances/${METRIC}-export"

    qiime tools export \
      --input-path \
      "$RESULTS/beta-first-distances/${METRIC}.qza" \
      --output-path \
      "$RESULTS/beta-first-distances/${METRIC}-export"

done


# ============================================================
# 5. BETA-DIVERSITY LINEAR MIXED-EFFECTS MODELS
# ============================================================

echo
echo "------------------------------------------------------------"
echo "BETA-DIVERSITY LONGITUDINAL MODELS"
echo "------------------------------------------------------------"

for METRIC in \
    bray_curtis \
    jaccard \
    weighted_unifrac \
    unweighted_unifrac

do

    echo
    echo "Modeling baseline displacement: $METRIC"

    qiime longitudinal linear-mixed-effects \
  	--m-metadata-file "$METADATA" \
  	--m-metadata-file \
  	"$RESULTS/beta-first-distances/${METRIC}.qza" \
  	--p-state-column days_after_fmt \
  	--p-individual-id-column participant \
  	--p-metric Distance \
  	--o-visualization \
  	"$RESULTS/beta-lme/${METRIC}_lme.qzv"

done


# ============================================================
# 6. STATISTICAL SUMMARY TABLES
# ============================================================

echo
echo "------------------------------------------------------------"
echo "CALCULATING PAIRED STATISTICAL SUMMARY"
echo "------------------------------------------------------------"

python - <<'PY'

import pandas as pd
import numpy as np

from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


df = pd.read_csv(
    "results/phase13/alpha-metadata.tsv",
    sep="\t"
)

metrics = [
    "observed_features",
    "shannon",
    "faith_pd",
    "evenness",
]

comparisons = [
    ("before", "7d"),
    ("before", "14d"),
    ("before", "30d"),
    ("7d", "14d"),
    ("7d", "30d"),
    ("14d", "30d"),
]

records = []

for metric in metrics:

    wide = df.pivot(
        index="participant",
        columns="timepoint",
        values=metric
    )

    metric_records = []

    for a, b in comparisons:

        pair = wide[[a, b]].dropna()

        x = pair[a]
        y = pair[b]

        delta = y - x

        try:

            stat, p = wilcoxon(
                y,
                x,
                alternative="two-sided"
            )

        except ValueError:

            stat = np.nan
            p = 1.0

        rec = {
            "metric": metric,
            "state_1": a,
            "state_2": b,
            "n_pairs": len(pair),
            "median_state_1": x.median(),
            "median_state_2": y.median(),
            "median_change": delta.median(),
            "wilcoxon_statistic": stat,
            "p_value": p,
        }

        metric_records.append(rec)

    pvals = [
        r["p_value"]
        for r in metric_records
    ]

    reject, qvals, _, _ = multipletests(
        pvals,
        alpha=0.05,
        method="fdr_bh"
    )

    for rec, q, sig in zip(
        metric_records,
        qvals,
        reject
    ):

        rec["q_value_bh"] = q
        rec["significant_fdr_0.05"] = sig

        records.append(rec)


result = pd.DataFrame(records)

result.to_csv(
    "results/phase13/tables/"
    "alpha_pairwise_statistics.tsv",
    sep="\t",
    index=False
)


print()
print("PAIRED ALPHA-DIVERSITY RESULTS")
print()

display = result.copy()

for col in [
    "median_state_1",
    "median_state_2",
    "median_change",
]:
    display[col] = display[col].round(3)

display["p_value"] = (
    display["p_value"]
    .map(lambda x: f"{x:.3e}")
)

display["q_value_bh"] = (
    display["q_value_bh"]
    .map(lambda x: f"{x:.3e}")
)

print(
    display[
        [
            "metric",
            "state_1",
            "state_2",
            "n_pairs",
            "median_state_1",
            "median_state_2",
            "median_change",
            "p_value",
            "q_value_bh",
            "significant_fdr_0.05",
        ]
    ].to_string(index=False)
)

PY


# ============================================================
# 7. BETA BASELINE-DISTANCE STATISTICS
# ============================================================

echo
echo "------------------------------------------------------------"
echo "BETA BASELINE-DISTANCE STATISTICS"
echo "------------------------------------------------------------"

python - <<'PY'

import pandas as pd

from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


df = pd.read_csv(
    "results/phase12/distance-from-baseline.tsv",
    sep="\t"
)

metrics = [
    "bray_curtis",
    "jaccard",
    "weighted_unifrac",
    "unweighted_unifrac",
]

comparisons = [
    ("7d", "14d"),
    ("7d", "30d"),
    ("14d", "30d"),
]

records = []

for metric in metrics:

    wide = df.pivot(
        index="participant",
        columns="timepoint",
        values=metric
    )

    metric_records = []

    for a, b in comparisons:

        pair = wide[[a, b]].dropna()

        x = pair[a]
        y = pair[b]

        delta = y - x

        try:

            stat, p = wilcoxon(
                y,
                x,
                alternative="two-sided"
            )

        except ValueError:

            stat = float("nan")
            p = 1.0

        rec = {
            "metric": metric,
            "state_1": a,
            "state_2": b,
            "n_pairs": len(pair),
            "median_state_1": x.median(),
            "median_state_2": y.median(),
            "median_change": delta.median(),
            "wilcoxon_statistic": stat,
            "p_value": p,
        }

        metric_records.append(rec)

    pvals = [
        r["p_value"]
        for r in metric_records
    ]

    reject, qvals, _, _ = multipletests(
        pvals,
        alpha=0.05,
        method="fdr_bh"
    )

    for rec, q, sig in zip(
        metric_records,
        qvals,
        reject
    ):

        rec["q_value_bh"] = q
        rec["significant_fdr_0.05"] = sig

        records.append(rec)


result = pd.DataFrame(records)

result.to_csv(
    "results/phase13/tables/"
    "beta_baseline_distance_statistics.tsv",
    sep="\t",
    index=False
)


print()
print("PAIRED BETA-DISTANCE RESULTS")
print()

display = result.copy()

for col in [
    "median_state_1",
    "median_state_2",
    "median_change",
]:
    display[col] = display[col].round(3)

display["p_value"] = (
    display["p_value"]
    .map(lambda x: f"{x:.3e}")
)

display["q_value_bh"] = (
    display["q_value_bh"]
    .map(lambda x: f"{x:.3e}")
)

print(
    display[
        [
            "metric",
            "state_1",
            "state_2",
            "n_pairs",
            "median_state_1",
            "median_state_2",
            "median_change",
            "p_value",
            "q_value_bh",
            "significant_fdr_0.05",
        ]
    ].to_string(index=False)
)

PY


# ============================================================
# 8. OUTPUT INVENTORY
# ============================================================

echo
echo "------------------------------------------------------------"
echo "PHASE 13 OUTPUTS"
echo "------------------------------------------------------------"

find "$RESULTS" \
    -maxdepth 3 \
    -type f \
    | sort


# ============================================================
# 9. CONCLUSION
# ============================================================

echo
echo "============================================================"
echo "PHASE 13 CONCLUSION"
echo "============================================================"

echo
echo "Repeated-measures longitudinal statistical analyses were"
echo "performed across the complete 52-participant FMT cohort."
echo
echo "Alpha-diversity changes were evaluated using paired"
echo "non-parametric tests and participant-level mixed-effects"
echo "models."
echo
echo "Benjamini-Hochberg correction was applied across the six"
echo "pairwise timepoint comparisons within each alpha metric."
echo
echo "Beta-diversity displacement was calculated relative to"
echo "each participant's own pre-FMT baseline."
echo
echo "Post-FMT baseline distances were compared across 7-, 14-,"
echo "and 30-day timepoints using paired longitudinal tests."
echo
echo "These analyses account for the repeated-measures structure"
echo "of the longitudinal cohort."
echo
echo "============================================================"
echo "END PHASE 13"
echo "============================================================"
