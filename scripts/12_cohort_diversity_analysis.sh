#!/usr/bin/env bash

# ============================================================
# PHASE 12 — COHORT ALPHA AND BETA DIVERSITY ANALYSIS
# ============================================================
#
# Purpose:
#   Calculate cohort-wide alpha and beta diversity across
#   the longitudinal FMT recipient cohort.
#
# Cohort:
#   52 participants
#   208 samples
#   4 timepoints:
#     before
#     7d
#     14d
#     30d
#
# Inputs:
#   results/phase10/table.qza
#   results/phase10/rep-seqs.qza
#   data/cohort/samples.tsv
#
# Rarefaction depth:
#   7,500 reads/sample
#
# Rationale:
#   This depth retains all 208 samples and all 52 complete
#   longitudinal participants.
#
# Outputs:
#   results/phase12/
#
# ============================================================

set -euo pipefail

RESULTS="results/phase12"
LOG_DIR="logs"

TABLE="results/phase10/table.qza"
REP_SEQS="results/phase10/rep-seqs.qza"
METADATA="data/cohort/samples.tsv"

SAMPLING_DEPTH=7500
THREADS=4

mkdir -p "$RESULTS"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/12_cohort_diversity_analysis.txt"

exec > >(tee "$LOG_FILE") 2>&1


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

echo "============================================================"
echo "PHASE 12 — COHORT ALPHA AND BETA DIVERSITY ANALYSIS"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

echo
echo "Configured sampling depth:"
echo "$SAMPLING_DEPTH"

echo
echo "Samples:"
tail -n +2 "$METADATA" | wc -l


# ------------------------------------------------------------
# 1. Rarefaction-depth audit
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "RAREFACTION DEPTH AUDIT"
echo "------------------------------------------------------------"

export SAMPLING_DEPTH

python - <<'PY'
import os
from pathlib import Path

import pandas as pd


FREQUENCY_FILE = Path(
    "results/phase10/"
    "sample-frequencies-export/"
    "metadata.tsv"
)

METADATA_FILE = Path(
    "data/cohort/samples.tsv"
)

OUTPUT_FILE = Path(
    "results/phase12/"
    "rarefaction_depth_audit.tsv"
)

SELECTED_DEPTH = int(
    os.environ["SAMPLING_DEPTH"]
)


# ------------------------------------------------------------
# Read post-DADA2 sequencing depths
# ------------------------------------------------------------

freq = pd.read_csv(
    FREQUENCY_FILE,
    sep="\t",
    comment="#"
)

meta = pd.read_csv(
    METADATA_FILE,
    sep="\t"
)

# QIIME exports values such as "22,050.0" as text.
freq = freq.iloc[:, :2].copy()

freq.columns = [
    "sample-id",
    "frequency",
]

freq["frequency"] = pd.to_numeric(
    freq["frequency"]
    .astype(str)
    .str.replace(
        ",",
        "",
        regex=False
    )
)


# ------------------------------------------------------------
# Merge sequencing depth with cohort metadata
# ------------------------------------------------------------

df = meta.merge(
    freq,
    on="sample-id",
    how="inner",
    validate="one_to_one"
)

if len(df) != len(meta):

    raise RuntimeError(
        "Post-DADA2 frequency table does not contain "
        "all cohort metadata samples."
    )


# ------------------------------------------------------------
# Evaluate candidate rarefaction depths
# ------------------------------------------------------------

depths = [
    5000,
    7500,
    10000,
    12500,
    15000,
    17500,
    20000,
]

# Ensure configured depth is always audited.
if SELECTED_DEPTH not in depths:
    depths.append(SELECTED_DEPTH)
    depths.sort()

records = []

total_samples = len(df)

total_participants = (
    df["participant"]
    .nunique()
)

for depth in depths:

    kept = df[
        df["frequency"] >= depth
    ]

    participant_counts = (
        kept
        .groupby("participant")["timepoint"]
        .nunique()
    )

    complete_participants = int(
        (
            participant_counts == 4
        ).sum()
    )

    records.append(
        {
            "sampling_depth": depth,
            "samples_retained": len(kept),
            "samples_total": total_samples,
            "complete_participants":
                complete_participants,
            "participants_total":
                total_participants,
        }
    )


audit = pd.DataFrame(
    records
)

audit.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Report results
# ------------------------------------------------------------

print("Post-DADA2 sequencing depth:")

print(
    df["frequency"]
    .describe()
    .to_string()
)

print()
print("Candidate rarefaction depths:")
print()

for _, row in audit.iterrows():

    print(
        f"{int(row['sampling_depth']):>6,} reads: "
        f"{int(row['samples_retained']):>3}/"
        f"{int(row['samples_total'])} samples retained | "
        f"{int(row['complete_participants']):>2}/"
        f"{int(row['participants_total'])} "
        "complete participants"
    )


# ------------------------------------------------------------
# Validate configured sampling depth
# ------------------------------------------------------------

selected = audit[
    audit["sampling_depth"] == SELECTED_DEPTH
].iloc[0]

if (
    selected["samples_retained"]
    != selected["samples_total"]
):

    raise RuntimeError(
        f"Configured rarefaction depth "
        f"{SELECTED_DEPTH:,} does not retain "
        "all cohort samples."
    )

if (
    selected["complete_participants"]
    != selected["participants_total"]
):

    raise RuntimeError(
        f"Configured rarefaction depth "
        f"{SELECTED_DEPTH:,} does not retain "
        "all complete longitudinal participants."
    )

print()
print(
    f"Selected rarefaction depth: "
    f"{SELECTED_DEPTH:,} reads/sample"
)

print(
    "Validation: all samples and all complete "
    "four-timepoint participants are retained."
)

print()
print(
    f"Audit table written to: {OUTPUT_FILE}"
)
PY


# ------------------------------------------------------------
# 2. Build phylogenetic tree
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "PHYLOGENETIC TREE"
echo "------------------------------------------------------------"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$REP_SEQS" \
  --p-n-threads "$THREADS" \
  --o-alignment "$RESULTS/aligned-rep-seqs.qza" \
  --o-masked-alignment "$RESULTS/masked-aligned-rep-seqs.qza" \
  --o-tree "$RESULTS/unrooted-tree.qza" \
  --o-rooted-tree "$RESULTS/rooted-tree.qza"


# ------------------------------------------------------------
# 3. Core diversity metrics
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "CORE DIVERSITY METRICS"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/core-metrics"

qiime diversity core-metrics-phylogenetic \
  --i-table "$TABLE" \
  --i-phylogeny "$RESULTS/rooted-tree.qza" \
  --p-sampling-depth "$SAMPLING_DEPTH" \
  --m-metadata-file "$METADATA" \
  --p-n-jobs-or-threads "$THREADS" \
  --output-dir "$RESULTS/core-metrics"


# ------------------------------------------------------------
# 4. Export alpha diversity vectors
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING ALPHA DIVERSITY"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/alpha-export"

mkdir -p "$RESULTS/alpha-export"

qiime tools export \
  --input-path "$RESULTS/core-metrics/observed_features_vector.qza" \
  --output-path "$RESULTS/alpha-export/observed-features"

qiime tools export \
  --input-path "$RESULTS/core-metrics/shannon_vector.qza" \
  --output-path "$RESULTS/alpha-export/shannon"

qiime tools export \
  --input-path "$RESULTS/core-metrics/faith_pd_vector.qza" \
  --output-path "$RESULTS/alpha-export/faith-pd"

qiime tools export \
  --input-path "$RESULTS/core-metrics/evenness_vector.qza" \
  --output-path "$RESULTS/alpha-export/evenness"


# ------------------------------------------------------------
# 5. Export beta diversity matrices
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "EXPORTING BETA DIVERSITY"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/beta-export"

mkdir -p "$RESULTS/beta-export"

qiime tools export \
  --input-path "$RESULTS/core-metrics/bray_curtis_distance_matrix.qza" \
  --output-path "$RESULTS/beta-export/bray-curtis"

qiime tools export \
  --input-path "$RESULTS/core-metrics/jaccard_distance_matrix.qza" \
  --output-path "$RESULTS/beta-export/jaccard"

qiime tools export \
  --input-path "$RESULTS/core-metrics/weighted_unifrac_distance_matrix.qza" \
  --output-path "$RESULTS/beta-export/weighted-unifrac"

qiime tools export \
  --input-path "$RESULTS/core-metrics/unweighted_unifrac_distance_matrix.qza" \
  --output-path "$RESULTS/beta-export/unweighted-unifrac"


# ------------------------------------------------------------
# 6. Alpha rarefaction
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
# 7. Cohort alpha-diversity summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "COHORT ALPHA DIVERSITY SUMMARY"
echo "------------------------------------------------------------"

python - <<'PY'
from pathlib import Path
import pandas as pd

meta = pd.read_csv(
    "data/cohort/samples.tsv",
    sep="\t"
)

files = {
    "observed_features":
        "results/phase12/alpha-export/"
        "observed-features/alpha-diversity.tsv",

    "shannon":
        "results/phase12/alpha-export/"
        "shannon/alpha-diversity.tsv",

    "faith_pd":
        "results/phase12/alpha-export/"
        "faith-pd/alpha-diversity.tsv",

    "evenness":
        "results/phase12/alpha-export/"
        "evenness/alpha-diversity.tsv",
}

combined = None

for metric, path in files.items():

    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0
    )

    df.columns = [metric]

    if combined is None:
        combined = df
    else:
        combined = combined.join(
            df,
            how="outer"
        )

combined.index.name = "sample-id"

combined = (
    combined
    .reset_index()
    .merge(
        meta,
        on="sample-id",
        how="left"
    )
)

order = {
    "before": 0,
    "7d": 1,
    "14d": 2,
    "30d": 3,
}

combined["timepoint_order"] = (
    combined["timepoint"]
    .map(order)
)

combined = combined.sort_values(
    [
        "participant",
        "timepoint_order",
    ]
)

out = (
    "results/phase12/"
    "alpha-diversity-summary.tsv"
)

combined.to_csv(
    out,
    sep="\t",
    index=False
)

print(
    f"Alpha-diversity samples: {len(combined)}"
)

print()
print("Samples per timepoint:")

print(
    combined["timepoint"]
    .value_counts()
    .reindex(
        [
            "before",
            "7d",
            "14d",
            "30d",
        ]
    )
    .to_string()
)

print()
print("Median alpha diversity by timepoint:")

summary = (
    combined
    .groupby("timepoint")[
        [
            "observed_features",
            "shannon",
            "faith_pd",
            "evenness",
        ]
    ]
    .median()
    .reindex(
        [
            "before",
            "7d",
            "14d",
            "30d",
        ]
    )
)

print(
    summary.to_string(
        float_format=lambda x: f"{x:.3f}"
    )
)
PY


# ------------------------------------------------------------
# 8. Distance from each participant's own baseline
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "WITHIN-PARTICIPANT DISTANCE FROM BASELINE"
echo "------------------------------------------------------------"

python - <<'PY'
from pathlib import Path
import pandas as pd

meta = pd.read_csv(
    "data/cohort/samples.tsv",
    sep="\t"
)

metric_files = {
    "bray_curtis":
        "results/phase12/beta-export/"
        "bray-curtis/distance-matrix.tsv",

    "jaccard":
        "results/phase12/beta-export/"
        "jaccard/distance-matrix.tsv",

    "weighted_unifrac":
        "results/phase12/beta-export/"
        "weighted-unifrac/distance-matrix.tsv",

    "unweighted_unifrac":
        "results/phase12/beta-export/"
        "unweighted-unifrac/distance-matrix.tsv",
}

meta_lookup = (
    meta
    .set_index(
        [
            "participant",
            "timepoint"
        ]
    )["sample-id"]
)

records = []

for participant in sorted(
    meta["participant"].unique()
):

    baseline = meta_lookup.loc[
        (
            participant,
            "before"
        )
    ]

    for timepoint in [
        "7d",
        "14d",
        "30d",
    ]:

        sample = meta_lookup.loc[
            (
                participant,
                timepoint
            )
        ]

        record = {
            "participant": participant,
            "timepoint": timepoint,
            "baseline_sample": baseline,
            "sample-id": sample,
        }

        for metric, path in metric_files.items():

            matrix = pd.read_csv(
                path,
                sep="\t",
                index_col=0
            )

            record[metric] = matrix.loc[
                baseline,
                sample
            ]

        records.append(
            record
        )

result = pd.DataFrame(
    records
)

result.to_csv(
    "results/phase12/"
    "distance-from-baseline.tsv",
    sep="\t",
    index=False
)

print(
    f"Baseline-to-post-FMT comparisons: {len(result)}"
)

print()
print("Median distance from baseline:")

summary = (
    result
    .groupby("timepoint")[
        [
            "bray_curtis",
            "jaccard",
            "weighted_unifrac",
            "unweighted_unifrac",
        ]
    ]
    .median()
    .reindex(
        [
            "7d",
            "14d",
            "30d",
        ]
    )
)

print(
    summary.to_string(
        float_format=lambda x: f"{x:.3f}"
    )
)
PY


# ------------------------------------------------------------
# 9. Final outputs
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "PHASE 12 OUTPUTS"
echo "------------------------------------------------------------"

find "$RESULTS" \
  -maxdepth 3 \
  -type f \
  | sort


# ------------------------------------------------------------
# 10. Conclusion
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 12 CONCLUSION"
echo "============================================================"

echo
echo "Cohort-wide alpha and beta diversity metrics were calculated"
echo "after rarefaction to 7,500 reads per sample."
echo
echo "This sampling depth retained all 208 samples and all 52"
echo "participants with complete four-timepoint longitudinal data."
echo
echo "Alpha diversity included observed features, Shannon"
echo "diversity, Faith's phylogenetic diversity, and evenness."
echo
echo "Beta diversity included Bray-Curtis, Jaccard, weighted"
echo "UniFrac, and unweighted UniFrac distances."
echo
echo "Within-participant distances from each recipient's own"
echo "pre-FMT baseline were calculated for 7d, 14d, and 30d."
echo
echo "These outputs are ready for repeated-measures cohort-level"
echo "statistical analysis."
echo
echo "============================================================"
echo "END PHASE 12"
echo "============================================================"
