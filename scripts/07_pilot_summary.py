#!/usr/bin/env python3

"""
PHASE 7 — INTEGRATED PILOT SUMMARY

Purpose
-------
Integrate the major results from the R009 longitudinal FMT
pilot into compact summary tables and repository-ready figures.

Inputs
------
results/phase2/denoising-stats-export/stats.tsv
results/phase4/alpha-export/observed-features/alpha-diversity.tsv
results/phase4/alpha-export/shannon/alpha-diversity.tsv
results/phase4/alpha-export/faith-pd/alpha-diversity.tsv
results/phase4/beta-export/bray-curtis/distance-matrix.tsv
results/phase4/beta-export/weighted-unifrac/distance-matrix.tsv
results/phase6/tables/top_genera_relative_abundance.tsv

Outputs
-------
results/phase7/tables/
results/phase7/figures/

Important
---------
This pilot contains one participant sampled longitudinally.
Results are descriptive and are not population-level
statistical inference.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuration
# ============================================================

RESULTS = Path("results/phase7")
TABLE_DIR = RESULTS / "tables"
FIGURE_DIR = RESULTS / "figures"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

DENOISING = Path("results/phase2/denoising-stats-export/stats.tsv")

OBSERVED = Path("results/phase4/alpha-export/"
        "observed-features/alpha-diversity.tsv")

SHANNON = Path("results/phase4/alpha-export/"
               "shannon/alpha-diversity.tsv")

FAITH = Path("results/phase4/alpha-export/"
             "faith-pd/alpha-diversity.tsv")

BRAY = Path("results/phase4/beta-export/"
            "bray-curtis/distance-matrix.tsv")

WEIGHTED_UNIFRAC = Path("results/phase4/beta-export/"
                "weighted-unifrac/distance-matrix.tsv")

TOP_GENERA = Path("results/phase6/tables/"
              "top_genera_relative_abundance.tsv")

TIMEPOINT_ORDER = ["R009_before", "R009_7d", "R009_14d", "R009_30d"]

TIMEPOINT_LABELS = {
    "R009_before": "Before",
    "R009_7d": "7 days",
    "R009_14d": "14 days",
    "R009_30d": "30 days"}


# ============================================================
# Helpers
# ============================================================

def read_alpha(path, value_name):
    """Read an exported QIIME alpha-diversity vector."""

    df = pd.read_csv(path, sep="\t", index_col=0)

    df.index.name = "sample-id"
    df.columns = [value_name]

    return df


def check_samples(df, name):
    """Confirm all expected longitudinal samples are present."""

    missing = [sample
        for sample in TIMEPOINT_ORDER
        if sample not in df.index]

    if missing:
        raise ValueError(f"{name} is missing expected samples: {missing}")


# ============================================================
# Start
# ============================================================

print("=" * 60)
print("PHASE 7 — INTEGRATED PILOT SUMMARY")
print("=" * 60)


# ============================================================
# 1. Sequencing retention
# ============================================================

print("\nLoading sequencing-retention statistics...")

denoise = pd.read_csv(DENOISING, sep="\t", comment="#")

denoise = denoise.set_index("sample-id")

check_samples(denoise, "DADA2 statistics")

retention = denoise.loc[TIMEPOINT_ORDER,
    ["input", "filtered", "denoised", "merged",
     "non-chimeric", "percentage of input non-chimeric"]].copy()

retention.index = [TIMEPOINT_LABELS[x]
    for x in retention.index]

retention.index.name = "timepoint"

retention.to_csv(TABLE_DIR / "sequencing_retention.tsv", sep="\t")


# ============================================================
# 2. Alpha diversity
# ============================================================

print("Loading alpha-diversity results...")

observed = read_alpha(OBSERVED, "observed_features")

shannon = read_alpha(SHANNON, "shannon")

faith = read_alpha(FAITH, "faith_pd")

alpha = observed.join(shannon).join(faith)

check_samples(alpha, "Alpha diversity")

alpha = alpha.loc[TIMEPOINT_ORDER].copy()

alpha.index = [TIMEPOINT_LABELS[x]
    for x in alpha.index]

alpha.index.name = "timepoint"

alpha.to_csv(TABLE_DIR / "alpha_diversity_summary.tsv", sep="\t")


# ============================================================
# 3. Baseline beta diversity
# ============================================================

print("Loading beta-diversity results...")

bray = pd.read_csv(BRAY, sep="\t", index_col=0)

weighted = pd.read_csv(WEIGHTED_UNIFRAC, sep="\t", index_col=0)

check_samples(bray, "Bray-Curtis matrix")
check_samples(weighted, "Weighted UniFrac matrix")

post_samples = ["R009_7d", "R009_14d", "R009_30d"]

beta_summary = pd.DataFrame(
    {"timepoint": [
        TIMEPOINT_LABELS[x]
          for x in post_samples ],
     "bray_curtis_vs_before": [
        bray.loc["R009_before", x]
        for x in post_samples],
      "weighted_unifrac_vs_before": [
         weighted.loc["R009_before", x]
         for x in post_samples]})

beta_summary.to_csv(TABLE_DIR / "distance_from_baseline.tsv",
    sep="\t", index=False)


# ============================================================
# 4. Dominant genera
# ============================================================

print("Loading dominant-genus results...")

genera = pd.read_csv(TOP_GENERA, sep="\t")

required = ["Genus"] + TIMEPOINT_ORDER

missing_columns = [column
    for column in required
    if column not in genera.columns]

if missing_columns:
    raise ValueError("Top-genera table is missing columns: "
        f"{missing_columns}" )

genera_summary = genera[required].copy()

for sample in TIMEPOINT_ORDER:
    genera_summary[sample] *= 100

genera_summary = genera_summary.rename(
    columns={sample: TIMEPOINT_LABELS[sample]
        for sample in TIMEPOINT_ORDER})

genera_summary.to_csv( TABLE_DIR / "dominant_genera_percent.tsv",
    sep="\t", index=False)


# ============================================================
# 5. Figure — Sequencing retention
# ============================================================

fig, ax = plt.subplots(figsize=(8, 5))

ax.bar(retention.index, retention["percentage of input non-chimeric"])

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Input reads retained (%)")
ax.set_title("DADA2 Read Retention Across the FMT Time Course")
ax.set_ylim(0, 100)

plt.tight_layout()

plt.savefig(FIGURE_DIR / "sequencing_retention.png", dpi=300)

plt.close()


# ============================================================
# 6. Figure — Alpha diversity trajectories
# ============================================================

x = np.arange(len(alpha.index))

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x, alpha["observed_features"], marker="o")

ax.set_xticks(x)
ax.set_xticklabels(alpha.index)

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Observed features")
ax.set_title("Observed ASV Richness Across the FMT Time Course")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "observed_features_trajectory.png", dpi=300)

plt.close()


fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x, alpha["shannon"], marker="o")

ax.set_xticks(x)
ax.set_xticklabels(alpha.index)

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Shannon diversity")
ax.set_title("Shannon Diversity Across the FMT Time Course")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "shannon_trajectory.png", dpi=300)

plt.close()


fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x, alpha["faith_pd"], marker="o")

ax.set_xticks(x)
ax.set_xticklabels(alpha.index)

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Faith's phylogenetic diversity")
ax.set_title("Faith's PD Across the FMT Time Course")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "faith_pd_trajectory.png", dpi=300)

plt.close()


# ============================================================
# 7. Figure — Distance from baseline
# ============================================================

x = np.arange(len(beta_summary))

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x, beta_summary["bray_curtis_vs_before"], marker="o")

ax.set_xticks(x)
ax.set_xticklabels(beta_summary["timepoint"])

ax.set_xlabel("Post-FMT timepoint")
ax.set_ylabel("Bray-Curtis distance from baseline")
ax.set_title("Compositional Distance from Pre-FMT Baseline")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "bray_curtis_from_baseline.png", dpi=300)

plt.close()


fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x, beta_summary["weighted_unifrac_vs_before"], marker="o")

ax.set_xticks(x)
ax.set_xticklabels(beta_summary["timepoint"])

ax.set_xlabel("Post-FMT timepoint")
ax.set_ylabel("Weighted UniFrac distance from baseline")
ax.set_title("Phylogenetic Distance from Pre-FMT Baseline")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "weighted_unifrac_from_baseline.png", dpi=300)

plt.close()


# ============================================================
# 8. Figure — Dominant genera
# ============================================================

plot_genera = (genera_summary.set_index("Genus")[
        ["Before", "7 days", "14 days", "30 days"]].T)

fig, ax = plt.subplots(figsize=(10, 7))

plot_genera.plot(kind="bar", stacked=True, ax=ax)

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Relative abundance (%)")
ax.set_title("Dominant Genera Across the FMT Time Course")

ax.legend(title="Genus", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.xticks(rotation=0)
plt.tight_layout()

plt.savefig(FIGURE_DIR / "dominant_genera_summary.png", dpi=300)

plt.close()


# ============================================================
# 9. Print integrated results
# ============================================================

print("\n" + "-" * 60)
print("SEQUENCING RETENTION")
print("-" * 60)

print(retention[["percentage of input non-chimeric"]].to_string(
        float_format=lambda x: f"{x:.2f}"))

print("\n" + "-" * 60)
print("ALPHA DIVERSITY")
print("-" * 60)

print(alpha.to_string(float_format=lambda x: f"{x:.3f}"))

print("\n" + "-" * 60)
print("DISTANCE FROM PRE-FMT BASELINE")
print("-" * 60)

print(beta_summary.to_string(index=False,
        float_format=lambda x: f"{x:.3f}"))

print("\n" + "-" * 60)
print("DOMINANT GENERA (%)")
print("-" * 60)

print(genera_summary.to_string(index=False,
        float_format=lambda x: f"{x:.2f}"))


# ============================================================
# 10. Outputs
# ============================================================

print("\n" + "=" * 60)
print("PHASE 7 OUTPUTS")
print("=" * 60)

for path in sorted(RESULTS.rglob("*")):
    if path.is_file():
        print(path)


# ============================================================
# 11. Conclusion
# ============================================================

print("\n" + "=" * 60)
print("PHASE 7 CONCLUSION")
print("=" * 60)

print("""The R009 pilot workflow was integrated across sequencing
quality, alpha diversity, beta diversity, and taxonomic
composition.

Post-FMT samples showed reduced richness and phylogenetic
diversity relative to the pre-FMT sample.

Post-FMT communities remained compositionally distinct
from the pre-FMT baseline across the observed time course.

Taxonomic analysis identified substantial longitudinal
changes among dominant bacterial genera.

Because the pilot consists of a single participant, these
patterns are interpreted descriptively and are intended to
validate the computational workflow before expansion to a
larger longitudinal cohort.""")

print("=" * 60)
print("END PHASE 7")
print("=" * 60)
