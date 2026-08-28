#!/usr/bin/env python3

"""
PHASE 6 — LONGITUDINAL TAXONOMIC CHANGE AND VISUALIZATION

Purpose
-------
Analyze longitudinal taxonomic relative abundance for the
R009 FMT pilot and generate repository-ready tables and figures.

Inputs
------
results/phase5/export/genus/genus-relative.tsv
results/phase5/export/phylum/phylum-relative.tsv

Outputs
-------
results/phase6/tables/
results/phase6/figures/

Important
---------
This pilot contains one participant sampled at four timepoints.
Results are descriptive and are not population-level statistical
or differential-abundance inference.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuration
# ============================================================

GENUS_FILE = Path("results/phase5/export/genus/genus-relative.tsv")

PHYLUM_FILE = Path("results/phase5/export/phylum/phylum-relative.tsv")

RESULTS = Path("results/phase6")
TABLE_DIR = RESULTS / "tables"
FIGURE_DIR = RESULTS / "figures"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

TIMEPOINT_ORDER = ["R009_before", "R009_7d", "R009_14d", "R009_30d"]

TIMEPOINT_LABELS = {
    "R009_before": "Before",
    "R009_7d": "7 days",
    "R009_14d": "14 days",
    "R009_30d": "30 days"}

TOP_N_GENERA = 12

# ============================================================
# Helper functions
# ============================================================

def read_biom_tsv(path):
    """
    Read a BIOM-generated TSV while preserving the actual
    sample names and column order written to the file.
    """

    with open(path, "r") as handle:
        lines = handle.readlines()

    header_line = None

    for line in lines:
        if line.startswith("#OTU ID"):
            header_line = line.strip()
            break

    if header_line is None:
        raise ValueError(
            f"Could not find '#OTU ID' header in {path}" )

    columns = header_line.split("\t")
    columns[0] = "Taxon"

    df = pd.read_csv(path, sep="\t", comment="#",
        header=None, names=columns)

    missing = [sample
        for sample in TIMEPOINT_ORDER
        if sample not in df.columns ]

    if missing:
        raise ValueError(
            f"Missing expected sample columns: {missing}")

    for sample in TIMEPOINT_ORDER:
        df[sample] = pd.to_numeric(
            df[sample],
            errors="raise")

    return df


def clean_taxon_name(taxon):
    """Return the most specific available taxonomic label."""

    parts = str(taxon).split(";")

    for part in reversed(parts):
        if "__" in part:
            name = part.split("__", 1)[1].strip()

            if name:
                return name

    return "Unclassified"


# ============================================================
# Load data
# ============================================================

print("=" * 60)
print("PHASE 6 — LONGITUDINAL TAXONOMIC CHANGE")
print("=" * 60)

print("\nLoading Phase 5 relative-abundance tables...")

genus = read_biom_tsv(GENUS_FILE)
phylum = read_biom_tsv(PHYLUM_FILE)

genus["Genus"] = genus["Taxon"].apply(clean_taxon_name)
phylum["Phylum"] = phylum["Taxon"].apply(clean_taxon_name)

print(f"Genus-level rows: {len(genus)}")
print(f"Phylum-level rows: {len(phylum)}")

print("\nDetected sample columns:")
for sample in TIMEPOINT_ORDER:
    print(f"  {sample}")


# ============================================================
# Rank genera by mean relative abundance
# ============================================================

genus["mean_relative_abundance"] = (genus[TIMEPOINT_ORDER].mean(axis=1))

genus_ranked = genus.sort_values(
    "mean_relative_abundance", ascending=False).copy()

genus_ranked.to_csv( TABLE_DIR / "genus_relative_abundance_ranked.tsv",
    sep="\t", index=False)


# ============================================================
# Top genera
# ============================================================

top = genus_ranked.head(TOP_N_GENERA).copy()

top_output = top[["Genus"] + TIMEPOINT_ORDER
    + ["mean_relative_abundance"]].copy()

top_output.to_csv(TABLE_DIR / "top_genera_relative_abundance.tsv",
    sep="\t", index=False)

print("\nTop genera by mean relative abundance (%):")

display_table = top_output.copy()

for column in TIMEPOINT_ORDER + ["mean_relative_abundance"]:
    display_table[column] *= 100

print( display_table.to_string(index=False,
        float_format=lambda x: f"{x:.2f}"))


# ============================================================
# Baseline-to-post-FMT changes
# ============================================================

changes = genus.copy()

changes["change_7d_vs_before"] = (
	changes["R009_7d"] - changes["R009_before"])

changes["change_14d_vs_before"] = (
	changes["R009_14d"] - changes["R009_before"])

changes["change_30d_vs_before"] = (
	changes["R009_30d"] - changes["R009_before"])

changes["max_absolute_change"] = (
    changes[
        ["change_7d_vs_before",
         "change_14d_vs_before",
         "change_30d_vs_before",
        ]].abs().max(axis=1))

changes = changes.sort_values("max_absolute_change", ascending=False)

changes[["Genus", "R009_before", "R009_7d", "R009_14d", "R009_30d",
        "change_7d_vs_before", "change_14d_vs_before",
        "change_30d_vs_before", "max_absolute_change",  ]
].to_csv(TABLE_DIR / "genus_baseline_changes.tsv", sep="\t", index=False)


# ============================================================
# Phylum table
# ============================================================

phylum_output = phylum[["Phylum"] + TIMEPOINT_ORDER].copy()

phylum_output.to_csv( TABLE_DIR / "phylum_relative_abundance.tsv",
    sep="\t", index=False)


# ============================================================
# Figure 1 — Phylum composition
# ============================================================

plot_df = (phylum.set_index("Phylum")[TIMEPOINT_ORDER].T)

plot_df.index = [ TIMEPOINT_LABELS[x]
    for x in plot_df.index]

ax = plot_df.plot(kind="bar", stacked=True, figsize=(9, 6))

ax.set_ylabel("Relative abundance")
ax.set_xlabel("FMT timepoint")
ax.set_title("Phylum-Level Microbiome Composition")

ax.legend(title="Phylum", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.xticks(rotation=0)
plt.tight_layout()

plt.savefig(FIGURE_DIR / "phylum_relative_abundance.png", dpi=300)

plt.close()


# ============================================================
# Figure 2 — Top genera stacked composition
# ============================================================

top_names = set(top["Genus"])

genus_plot = genus.copy()

genus_plot["Plot_group"] = np.where(
    genus_plot["Genus"].isin(top_names),
    genus_plot["Genus"], "Other")

genus_plot = ( genus_plot .groupby("Plot_group")[TIMEPOINT_ORDER].sum().T)

genus_plot.index = [ TIMEPOINT_LABELS[x]
    for x in genus_plot.index]

ax = genus_plot.plot(kind="bar", stacked=True, figsize=(10, 7))

ax.set_ylabel("Relative abundance")
ax.set_xlabel("FMT timepoint")
ax.set_title(f"Top {TOP_N_GENERA} Genera Across the FMT Time Course")

ax.legend(title="Genus", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.xticks(rotation=0)
plt.tight_layout()

plt.savefig(FIGURE_DIR / "top_genera_relative_abundance.png", dpi=300)

plt.close()


# ============================================================
# Figure 3 — Genus trajectories
# ============================================================

trajectory = (top.set_index("Genus")[TIMEPOINT_ORDER])

x = np.arange(len(TIMEPOINT_ORDER))

plt.figure(figsize=(10, 7))

for genus_name, row in trajectory.iterrows():

    plt.plot(x, row.values * 100, marker="o",
        label=genus_name)

plt.xticks(x, [TIMEPOINT_LABELS[t]
        for t in TIMEPOINT_ORDER])

plt.xlabel("FMT timepoint")
plt.ylabel("Relative abundance (%)")
plt.title("Longitudinal Trajectories of Dominant Genera")

plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "top_genera_trajectories.png", dpi=300)

plt.close()


# ============================================================
# Figure 4 — Heatmap
# ============================================================

heatmap = (top.set_index("Genus")[TIMEPOINT_ORDER]* 100)

heatmap.columns = [TIMEPOINT_LABELS[x]
    for x in heatmap.columns]

fig, ax = plt.subplots(figsize=(7, max(5, TOP_N_GENERA * 0.45)))

image = ax.imshow(heatmap.values, aspect="auto")

ax.set_xticks(np.arange(len(heatmap.columns)))

ax.set_xticklabels(heatmap.columns)

ax.set_yticks(np.arange(len(heatmap.index)))

ax.set_yticklabels(heatmap.index)

ax.set_xlabel("FMT timepoint")
ax.set_ylabel("Genus")
ax.set_title("Relative Abundance of Dominant Genera")

cbar = fig.colorbar(image, ax=ax)

cbar.set_label("Relative abundance (%)")

plt.tight_layout()

plt.savefig(FIGURE_DIR / "top_genera_heatmap.png", dpi=300)

plt.close()


# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 60)
print("PHASE 6 OUTPUTS")
print("=" * 60)

for path in sorted(RESULTS.rglob("*")):
    if path.is_file():
        print(path)

print("\n" + "=" * 60)
print("PHASE 6 CONCLUSION")
print("=" * 60)

print("""Longitudinal taxonomic composition was summarized across
the four R009 FMT timepoints.

The original BIOM-exported sample labels were preserved and
then explicitly reordered as Before, 7 days, 14 days, and
30 days for longitudinal analysis.

Genera were ranked by mean relative abundance across the
entire time course.

Baseline-to-post-FMT abundance changes were calculated for
7-, 14-, and 30-day samples.

Phylum composition, dominant-genus composition, genus
trajectories, and a taxonomic heatmap were generated.

Because this pilot contains only one participant, these
results represent descriptive longitudinal changes rather
than population-level differential-abundance inference.""")

print("=" * 60)
print("END PHASE 6")
print("=" * 60)
