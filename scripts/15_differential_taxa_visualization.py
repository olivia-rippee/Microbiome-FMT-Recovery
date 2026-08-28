#!/usr/bin/env python3

"""
PHASE 15 — DIFFERENTIAL TAXA VISUALIZATION

Purpose
-------
Visualize genus-level ANCOM-BC differential-abundance
results across the longitudinal FMT cohort.

Inputs
------
results/phase14/tables/genus_ancombc_all_results.tsv
results/phase14/tables/genus_ancombc_persistence.tsv

Outputs
-------
results/phase15/figures/
results/phase15/tables/

Figures
-------
1. top_persistent_differential_genera.png
2. genus_lfc_by_timepoint.png
3. persistent_genera_heatmap.png

Tables
------
1. top_persistent_genera.tsv
2. selected_genus_lfc_matrix.tsv
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# Configuration
# ============================================================

INPUT_ALL = Path("results/phase14/tables/"
	"genus_ancombc_all_results.tsv")

INPUT_PERSISTENCE = Path("results/phase14/tables/"
	"genus_ancombc_persistence.tsv")

RESULTS = Path("results/phase15")
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

TOP_N = 20

TIMEPOINT_ORDER = ["7d_vs_before", "14d_vs_before", "30d_vs_before"]


# ============================================================
# Helpers
# ============================================================

def short_taxon_name(taxon):
    """
    Extract the most specific available taxonomic label.
    """

    parts = str(taxon).split(";")

    for part in reversed(parts):

        part = part.strip()

        if "__" not in part:
            continue

        prefix, value = part.split("__", 1)

        if value:

            return value

    return str(taxon)


# ============================================================
# Load data
# ============================================================

print("=" * 80)
print("PHASE 15 — DIFFERENTIAL TAXA VISUALIZATION")
print("=" * 80)

all_results = pd.read_csv(INPUT_ALL, sep="\t")

persistence = pd.read_csv(INPUT_PERSISTENCE, sep="\t")

print()
print(f"ANCOM-BC result rows: {len(all_results)}")
print(f"Persistence rows: {len(persistence)}")


# ============================================================
# Keep taxa significant at all three timepoints
# ============================================================

persistent3 = persistence[
	persistence["significant_timepoints"] == 3].copy()

print()
print("Genera significant at all three post-FMT "
    f"timepoints: {len(persistent3)}")


# ============================================================
# Build per-taxon LFC summary
# ============================================================

lfc = all_results[
	all_results["comparison"].isin(TIMEPOINT_ORDER)].copy()

matrix = (lfc.pivot(index="taxon", columns="comparison", values="lfc")
    .reindex(columns=TIMEPOINT_ORDER))

matrix["mean_lfc"] = (matrix[TIMEPOINT_ORDER].mean(axis=1))

matrix = matrix.reset_index()

matrix = matrix.merge(
    persistent3[["taxon", "significant_timepoints", "min_q"]],
    on="taxon", how="inner")

matrix["label"] = (matrix["taxon"].apply(short_taxon_name))

matrix = matrix.sort_values("mean_lfc")


# ============================================================
# Select strongest increases and decreases
# ============================================================

decreased = (matrix[matrix["mean_lfc"] < 0].nsmallest(TOP_N // 2, "mean_lfc"))

increased = (matrix[matrix["mean_lfc"] > 0]
    .nlargest(TOP_N // 2, "mean_lfc"))

selected = pd.concat([decreased, increased], ignore_index=True)

selected = selected.sort_values("mean_lfc")

selected.to_csv( TABLES / "top_persistent_genera.tsv", sep="\t", index=False)

selected[["taxon", "label", *TIMEPOINT_ORDER, "mean_lfc", "min_q"]].to_csv(
	TABLES / "selected_genus_lfc_matrix.tsv", sep="\t", index=False)


# ============================================================
# Figure 1 — persistent genera mean LFC
# ============================================================

fig, ax = plt.subplots(figsize=(10, 8))

ax.barh(selected["label"], selected["mean_lfc"])

ax.axvline(0, linewidth=1)

ax.set_xlabel("Mean ANCOM-BC log-fold change\n"
    "(post-FMT vs before)")

ax.set_ylabel("Genus")

ax.set_title("Persistent differential genera after FMT")

fig.tight_layout()

fig.savefig(FIGURES / "top_persistent_differential_genera.png",
    dpi=300, bbox_inches="tight")

plt.close(fig)


# ============================================================
# Figure 2 — LFC by timepoint
# ============================================================

plot_data = selected.copy()

x = np.arange(len(plot_data))

width = 0.25

fig, ax = plt.subplots(figsize=(14, 7))

for i, comparison in enumerate(TIMEPOINT_ORDER):

    ax.bar(
        x + (i - 1) * width,
        plot_data[comparison],
        width,
        label=comparison.replace("_vs_before", ""))

ax.axhline(0, linewidth=1)

ax.set_xticks(x)

ax.set_xticklabels(plot_data["label"], rotation=60, ha="right")

ax.set_ylabel("ANCOM-BC log-fold change")

ax.set_xlabel("Genus")

ax.set_title("Differential abundance across post-FMT timepoints")

ax.legend(title="Timepoint")

fig.tight_layout()

fig.savefig(FIGURES / "genus_lfc_by_timepoint.png",
    dpi=300, bbox_inches="tight")

plt.close(fig)


# ============================================================
# Figure 3 — heatmap
# ============================================================

heat = selected.set_index("label")[TIMEPOINT_ORDER].copy()

fig, ax = plt.subplots(figsize=(8, 9))

image = ax.imshow(heat.values, aspect="auto")

ax.set_xticks(range(len(TIMEPOINT_ORDER)))

ax.set_xticklabels(["7d", "14d", "30d"])

ax.set_yticks(range(len(heat.index)))

ax.set_yticklabels(heat.index)

ax.set_xlabel("Post-FMT timepoint")

ax.set_ylabel("Genus")

ax.set_title("Persistent genus-level ANCOM-BC effects")

cbar = fig.colorbar(image, ax=ax)

cbar.set_label("Log-fold change")

fig.tight_layout()

fig.savefig(FIGURES / "persistent_genera_heatmap.png",
    dpi=300, bbox_inches="tight")

plt.close(fig)


# ============================================================
# Console summary
# ============================================================

print()
print("------------------------------------------------------------")
print("TOP PERSISTENT INCREASES")
print("------------------------------------------------------------")

print(selected[selected["mean_lfc"] > 0][["label", "mean_lfc", "min_q"]]
    .sort_values("mean_lfc", ascending=False).to_string(index=False))

print()
print("------------------------------------------------------------")
print("TOP PERSISTENT DECREASES")
print("------------------------------------------------------------")

print(selected[selected["mean_lfc"] < 0][["label", "mean_lfc", "min_q"]]
    .sort_values("mean_lfc").to_string( index=False))


# ============================================================
# Output inventory
# ============================================================

print()
print("============================================================")
print("PHASE 15 OUTPUTS")
print("============================================================")

for path in sorted(RESULTS.rglob("*")):

    if path.is_file():

        print(path)


# ============================================================
# Conclusion
# ============================================================

print()
print("============================================================")
print("PHASE 15 CONCLUSION")
print("============================================================")

print( """Genus-level ANCOM-BC results were summarized across the
three post-FMT timepoints.

Taxa that were significantly different from pre-FMT
baseline at all three timepoints were identified as
persistent differential genera.

The strongest persistent increases and decreases were
visualized using mean log-fold change, timepoint-specific
log-fold changes, and a heatmap.

These figures summarize sustained taxonomic restructuring
following FMT across the 52-participant longitudinal cohort.""")

print("=" * 80)
print("END PHASE 15")
print("=" * 80)
