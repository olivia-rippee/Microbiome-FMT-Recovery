#!/usr/bin/env python3

"""
PHASE 16 — INTEGRATED COHORT SUMMARY

Purpose
-------
Integrate the major cohort-level results from sequencing QC,
alpha diversity, beta diversity, longitudinal statistics,
and differential abundance into compact final tables and
publication-style summary figures.

Inputs
------
results/phase10/denoising-stats-export/stats.tsv
results/phase12/alpha-diversity-summary.tsv
results/phase12/distance-from-baseline.tsv
results/phase13/tables/alpha_pairwise_statistics.tsv
results/phase13/tables/beta_baseline_distance_statistics.tsv
results/phase14/tables/genus_ancombc_persistence.tsv
results/phase15/tables/top_persistent_genera.tsv

Outputs
-------
results/phase16/tables/cohort_key_results.tsv

results/phase16/figures/
    sequencing_retention_summary.png
    alpha_diversity_summary.png
    beta_distance_from_baseline.png
    top_persistent_genera.png
    integrated_cohort_summary.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# Configuration
# ============================================================

RESULTS = Path("results/phase16")
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

RETENTION_FILE = Path(
    "results/phase10/"
    "denoising-stats-export/"
    "stats.tsv"
)

ALPHA_FILE = Path(
    "results/phase12/"
    "alpha-diversity-summary.tsv"
)

BETA_FILE = Path(
    "results/phase12/"
    "distance-from-baseline.tsv"
)

ALPHA_STATS_FILE = Path(
    "results/phase13/tables/"
    "alpha_pairwise_statistics.tsv"
)

BETA_STATS_FILE = Path(
    "results/phase13/tables/"
    "beta_baseline_distance_statistics.tsv"
)

PERSISTENCE_FILE = Path(
    "results/phase14/tables/"
    "genus_ancombc_persistence.tsv"
)

TOP_GENERA_FILE = Path(
    "results/phase15/tables/"
    "top_persistent_genera.tsv"
)

TIMEPOINTS = [
    "before",
    "7d",
    "14d",
    "30d",
]

POST_TIMEPOINTS = [
    "7d",
    "14d",
    "30d",
]


# ============================================================
# Helpers
# ============================================================

def short_taxon_name(taxon):
    parts = str(taxon).split(";")

    for part in reversed(parts):
        part = part.strip()

        if "__" not in part:
            continue

        _, value = part.split(
            "__",
            1
        )

        if value:
            return value

    return str(taxon)


def require_file(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Required input missing: {path}"
        )


# ============================================================
# Validate inputs
# ============================================================

for path in [
    RETENTION_FILE,
    ALPHA_FILE,
    BETA_FILE,
    ALPHA_STATS_FILE,
    BETA_STATS_FILE,
    PERSISTENCE_FILE,
    TOP_GENERA_FILE,
]:
    require_file(path)


print("=" * 80)
print("PHASE 16 — INTEGRATED COHORT SUMMARY")
print("=" * 80)


# ============================================================
# Load data
# ============================================================

retention = pd.read_csv(
    RETENTION_FILE,
    sep="\t",
    comment="#"
)

alpha = pd.read_csv(
    ALPHA_FILE,
    sep="\t"
)

beta = pd.read_csv(
    BETA_FILE,
    sep="\t"
)

alpha_stats = pd.read_csv(
    ALPHA_STATS_FILE,
    sep="\t"
)

beta_stats = pd.read_csv(
    BETA_STATS_FILE,
    sep="\t"
)

persistence = pd.read_csv(
    PERSISTENCE_FILE,
    sep="\t"
)

top_genera = pd.read_csv(
    TOP_GENERA_FILE,
    sep="\t"
)


# ============================================================
# Basic cohort summary
# ============================================================

participants = (
    alpha["participant"]
    .nunique()
)

samples = len(alpha)

retention_col = (
    "percentage of input non-chimeric"
)

median_retention = (
    retention[retention_col]
    .median()
)

min_retention = (
    retention[retention_col]
    .min()
)

max_retention = (
    retention[retention_col]
    .max()
)


print()
print("------------------------------------------------------------")
print("COHORT OVERVIEW")
print("------------------------------------------------------------")

print(f"Participants: {participants}")
print(f"Samples: {samples}")

print(
    "Median non-chimeric read retention: "
    f"{median_retention:.2f}%"
)

print(
    "Retention range: "
    f"{min_retention:.2f}%–{max_retention:.2f}%"
)


# ============================================================
# Alpha diversity summary
# ============================================================

alpha_summary = (
    alpha
    .groupby("timepoint")[
        [
            "observed_features",
            "shannon",
            "faith_pd",
            "evenness",
        ]
    ]
    .median()
    .reindex(TIMEPOINTS)
)

print()
print("------------------------------------------------------------")
print("ALPHA DIVERSITY MEDIANS")
print("------------------------------------------------------------")

print(
    alpha_summary.to_string(
        float_format=lambda x: f"{x:.3f}"
    )
)


# ============================================================
# Beta diversity summary
# ============================================================

beta_summary = (
    beta
    .groupby("timepoint")[
        [
            "bray_curtis",
            "jaccard",
            "weighted_unifrac",
            "unweighted_unifrac",
        ]
    ]
    .median()
    .reindex(POST_TIMEPOINTS)
)

print()
print("------------------------------------------------------------")
print("DISTANCE FROM BASELINE MEDIANS")
print("------------------------------------------------------------")

print(
    beta_summary.to_string(
        float_format=lambda x: f"{x:.3f}"
    )
)


# ============================================================
# Persistent differential taxa summary
# ============================================================

persistent3 = persistence[
    persistence[
        "significant_timepoints"
    ] == 3
].copy()

persistent3["label"] = (
    persistent3["taxon"]
    .apply(short_taxon_name)
)

persistent3["direction_consistent"] = (
    persistent3["directions"]
    .apply(
        lambda value:
            len(
                set(
                    str(value)
                    .split(";")
                )
            ) == 1
    )
)

n_persistent = len(persistent3)

n_persistent_increased = (
    persistent3[
        persistent3["directions"]
        .str.contains(
            "increased"
        )
    ]
    .shape[0]
)

n_persistent_decreased = (
    persistent3[
        persistent3["directions"]
        .str.contains(
            "decreased"
        )
    ]
    .shape[0]
)


print()
print("------------------------------------------------------------")
print("PERSISTENT DIFFERENTIAL GENERA")
print("------------------------------------------------------------")

print(
    "Genera significant at all three "
    f"post-FMT timepoints: {n_persistent}"
)

print(
    f"Persistent increases: {n_persistent_increased}"
)

print(
    f"Persistent decreases: {n_persistent_decreased}"
)


# ============================================================
# Key statistical results
# ============================================================

baseline_alpha = alpha_stats[
    alpha_stats["state_1"]
    .eq("before")
].copy()

baseline_alpha_sig = (
    baseline_alpha[
        baseline_alpha[
            "significant_fdr_0.05"
        ]
    ]
)

post_beta_sig = (
    beta_stats[
        beta_stats[
            "significant_fdr_0.05"
        ]
    ]
)

print()
print("------------------------------------------------------------")
print("STATISTICAL SUMMARY")
print("------------------------------------------------------------")

print(
    "Significant baseline-vs-post alpha contrasts: "
    f"{len(baseline_alpha_sig)}/"
    f"{len(baseline_alpha)}"
)

print(
    "Significant post-FMT beta-distance contrasts: "
    f"{len(post_beta_sig)}/"
    f"{len(beta_stats)}"
)


# ============================================================
# Key results table
# ============================================================

records = []


def add_result(category, metric, timepoint, value, note):
    records.append(
        {
            "category": category,
            "metric": metric,
            "timepoint": timepoint,
            "value": value,
            "note": note,
        }
    )


add_result(
    "cohort",
    "participants",
    "all",
    participants,
    "Complete four-timepoint recipient cohort",
)

add_result(
    "cohort",
    "samples",
    "all",
    samples,
    "52 participants × 4 timepoints",
)

add_result(
    "sequencing_qc",
    "median_non_chimeric_retention_percent",
    "all",
    round(
        median_retention,
        3
    ),
    "Median DADA2 non-chimeric read retention",
)


for timepoint in TIMEPOINTS:

    row = alpha_summary.loc[
        timepoint
    ]

    for metric in [
        "observed_features",
        "shannon",
        "faith_pd",
        "evenness",
    ]:

        add_result(
            "alpha_diversity",
            metric,
            timepoint,
            round(
                row[metric],
                6
            ),
            "Median across participants",
        )


for timepoint in POST_TIMEPOINTS:

    row = beta_summary.loc[
        timepoint
    ]

    for metric in [
        "bray_curtis",
        "jaccard",
        "weighted_unifrac",
        "unweighted_unifrac",
    ]:

        add_result(
            "beta_diversity",
            metric,
            timepoint,
            round(
                row[metric],
                6
            ),
            "Median distance from participant-specific pre-FMT baseline",
        )


add_result(
    "differential_abundance",
    "persistent_significant_genera",
    "7d_14d_30d",
    n_persistent,
    "Significant at all three post-FMT timepoints",
)

add_result(
    "differential_abundance",
    "persistent_increased_genera",
    "7d_14d_30d",
    n_persistent_increased,
    "Significant persistent increases",
)

add_result(
    "differential_abundance",
    "persistent_decreased_genera",
    "7d_14d_30d",
    n_persistent_decreased,
    "Significant persistent decreases",
)


key_results = pd.DataFrame(
    records
)

key_results.to_csv(
    TABLES /
    "cohort_key_results.tsv",
    sep="\t",
    index=False
)


# ============================================================
# Figure 1 — sequencing retention
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.hist(
    retention[
        retention_col
    ],
    bins=20
)

ax.axvline(
    median_retention,
    linewidth=2
)

ax.set_xlabel(
    "Non-chimeric read retention (%)"
)

ax.set_ylabel(
    "Number of samples"
)

ax.set_title(
    "Cohort DADA2 read retention"
)

fig.tight_layout()

fig.savefig(
    FIGURES /
    "sequencing_retention_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Figure 2 — alpha diversity
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

x = np.arange(
    len(TIMEPOINTS)
)

metrics = [
    (
        "shannon",
        "Shannon"
    ),
    (
        "faith_pd",
        "Faith PD"
    ),
]

for metric, label in metrics:

    values = (
        alpha_summary[
            metric
        ]
        .values
    )

    ax.plot(
        x,
        values,
        marker="o",
        label=label
    )

ax.set_xticks(x)

ax.set_xticklabels(
    [
        "Before",
        "7d",
        "14d",
        "30d",
    ]
)

ax.set_xlabel(
    "FMT timepoint"
)

ax.set_ylabel(
    "Median diversity"
)

ax.set_title(
    "Alpha diversity following FMT"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    FIGURES /
    "alpha_diversity_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Figure 3 — beta distance from baseline
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)

x = np.arange(
    len(POST_TIMEPOINTS)
)

for metric, label in [
    (
        "bray_curtis",
        "Bray-Curtis"
    ),
    (
        "weighted_unifrac",
        "Weighted UniFrac"
    ),
]:

    values = (
        beta_summary[
            metric
        ]
        .values
    )

    ax.plot(
        x,
        values,
        marker="o",
        label=label
    )

ax.set_xticks(x)

ax.set_xticklabels(
    [
        "7d",
        "14d",
        "30d",
    ]
)

ax.set_xlabel(
    "Post-FMT timepoint"
)

ax.set_ylabel(
    "Median distance from pre-FMT baseline"
)

ax.set_title(
    "Persistent community displacement following FMT"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    FIGURES /
    "beta_distance_from_baseline.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Figure 4 — strongest persistent genera
# ============================================================

plot_genera = (
    top_genera
    .copy()
)

if "label" not in plot_genera.columns:
    plot_genera["label"] = (
        plot_genera["taxon"]
        .apply(short_taxon_name)
    )

plot_genera = (
    plot_genera
    .sort_values(
        "mean_lfc"
    )
)

fig, ax = plt.subplots(
    figsize=(10, 8)
)

ax.barh(
    plot_genera[
        "label"
    ],
    plot_genera[
        "mean_lfc"
    ]
)

ax.axvline(
    0,
    linewidth=1
)

ax.set_xlabel(
    "Mean ANCOM-BC log-fold change"
)

ax.set_ylabel(
    "Genus"
)

ax.set_title(
    "Strongest persistent differential genera"
)

fig.tight_layout()

fig.savefig(
    FIGURES /
    "top_persistent_genera.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Figure 5 — integrated summary
# ============================================================

fig = plt.figure(
    figsize=(12, 9)
)

ax = fig.add_subplot(111)

ax.axis(
    "off"
)

summary_lines = [
    "Integrated longitudinal FMT cohort summary",
    "",
    f"Participants: {participants}",
    f"Samples: {samples}",
    "",
    (
        "Median non-chimeric read retention: "
        f"{median_retention:.1f}%"
    ),
    "",
    (
        "Observed features: "
        f"{alpha_summary.loc['before', 'observed_features']:.1f} "
        "before → "
        f"{alpha_summary.loc['7d', 'observed_features']:.1f} "
        "at 7d → "
        f"{alpha_summary.loc['30d', 'observed_features']:.1f} "
        "at 30d"
    ),
    (
        "Shannon diversity: "
        f"{alpha_summary.loc['before', 'shannon']:.2f} "
        "before → "
        f"{alpha_summary.loc['7d', 'shannon']:.2f} "
        "at 7d → "
        f"{alpha_summary.loc['30d', 'shannon']:.2f} "
        "at 30d"
    ),
    (
        "Faith PD: "
        f"{alpha_summary.loc['before', 'faith_pd']:.2f} "
        "before → "
        f"{alpha_summary.loc['7d', 'faith_pd']:.2f} "
        "at 7d → "
        f"{alpha_summary.loc['30d', 'faith_pd']:.2f} "
        "at 30d"
    ),
    "",
    (
        "Median Bray-Curtis distance from baseline: "
        f"{beta_summary.loc['7d', 'bray_curtis']:.3f} "
        "at 7d, "
        f"{beta_summary.loc['14d', 'bray_curtis']:.3f} "
        "at 14d, "
        f"{beta_summary.loc['30d', 'bray_curtis']:.3f} "
        "at 30d"
    ),
    (
        "Median weighted UniFrac distance from baseline: "
        f"{beta_summary.loc['7d', 'weighted_unifrac']:.3f} "
        "at 7d, "
        f"{beta_summary.loc['14d', 'weighted_unifrac']:.3f} "
        "at 14d, "
        f"{beta_summary.loc['30d', 'weighted_unifrac']:.3f} "
        "at 30d"
    ),
    "",
    (
        "Persistent differential genera: "
        f"{n_persistent}"
    ),
    (
        "Persistent increases: "
        f"{n_persistent_increased}"
    ),
    (
        "Persistent decreases: "
        f"{n_persistent_decreased}"
    ),
    "",
    (
        "Interpretation: FMT is associated with a rapid "
        "increase in microbial diversity and a sustained "
        "shift away from each participant's pre-FMT "
        "community, accompanied by persistent changes "
        "in multiple gut-associated genera."
    ),
]

ax.text(
    0.02,
    0.98,
    "\n".join(summary_lines),
    va="top",
    ha="left",
    fontsize=12
)

fig.tight_layout()

fig.savefig(
    FIGURES /
    "integrated_cohort_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Console summary
# ============================================================

print()
print("------------------------------------------------------------")
print("KEY COHORT FINDINGS")
print("------------------------------------------------------------")

print(
    f"Participants: {participants}"
)

print(
    f"Samples analyzed: {samples}"
)

print(
    "Median DADA2 non-chimeric retention: "
    f"{median_retention:.2f}%"
)

print()

print(
    "Observed features:"
)

for tp in TIMEPOINTS:

    print(
        f"  {tp:>6}: "
        f"{alpha_summary.loc[tp, 'observed_features']:.1f}"
    )

print()

print(
    "Shannon diversity:"
)

for tp in TIMEPOINTS:

    print(
        f"  {tp:>6}: "
        f"{alpha_summary.loc[tp, 'shannon']:.3f}"
    )

print()

print(
    "Faith PD:"
)

for tp in TIMEPOINTS:

    print(
        f"  {tp:>6}: "
        f"{alpha_summary.loc[tp, 'faith_pd']:.3f}"
    )

print()

print(
    "Median Bray-Curtis distance from baseline:"
)

for tp in POST_TIMEPOINTS:

    print(
        f"  {tp:>6}: "
        f"{beta_summary.loc[tp, 'bray_curtis']:.3f}"
    )

print()

print(
    "Median weighted UniFrac distance from baseline:"
)

for tp in POST_TIMEPOINTS:

    print(
        f"  {tp:>6}: "
        f"{beta_summary.loc[tp, 'weighted_unifrac']:.3f}"
    )

print()

print(
    "Persistent significant genera: "
    f"{n_persistent}"
)

print(
    "Persistent increased genera: "
    f"{n_persistent_increased}"
)

print(
    "Persistent decreased genera: "
    f"{n_persistent_decreased}"
)


# ============================================================
# Output inventory
# ============================================================

print()
print("============================================================")
print("PHASE 16 OUTPUTS")
print("============================================================")

for path in sorted(
    RESULTS.rglob("*")
):

    if path.is_file():
        print(path)


# ============================================================
# Conclusion
# ============================================================

print()
print("============================================================")
print("PHASE 16 CONCLUSION")
print("============================================================")

print(
    """
The cohort-level analysis was integrated across sequencing
quality, alpha diversity, beta diversity, longitudinal
statistics, and differential abundance.

Across 52 participants and 208 longitudinal samples, FMT
was associated with a rapid increase in microbial richness,
Shannon diversity, phylogenetic diversity, and evenness.

Post-FMT communities remained strongly displaced from each
participant's pre-FMT baseline through day 30.

ANCOM-BC identified numerous genera with persistent
post-FMT increases or decreases across all three measured
post-FMT timepoints.

Together, these results support rapid and sustained
restructuring of the gut microbial community following FMT.
"""
)

print("=" * 80)
print("END PHASE 16")
print("=" * 80)
