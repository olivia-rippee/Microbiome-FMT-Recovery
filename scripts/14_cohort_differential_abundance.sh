#!/usr/bin/env bash

# ============================================================
# PHASE 14 — COHORT DIFFERENTIAL ABUNDANCE
# ============================================================
#
# Purpose:
#   Identify taxa that change systematically after FMT across
#   the complete 52-participant longitudinal cohort.
#
# Method:
#   ANCOM-BC
#
# Model:
#   abundance ~ timepoint + participant
#
# Reference:
#   timepoint = before
#
# This design adjusts for participant-specific baseline
# differences while estimating changes at:
#   7d vs before
#   14d vs before
#   30d vs before
#
# Primary taxonomic level:
#   Genus (SILVA level 6)
#
# Secondary level:
#   Family (SILVA level 5)
#
# Inputs:
#   results/phase10/table.qza
#   results/phase11/taxonomy.qza
#   data/cohort/samples.tsv
#
# ============================================================

set -euo pipefail


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

RESULTS="results/phase14"
LOG_DIR="logs"

TABLE="results/phase10/table.qza"
TAXONOMY="results/phase11/taxonomy.qza"
METADATA="data/cohort/samples.tsv"

PRV_CUT=0.10
ALPHA=0.05

mkdir -p "$RESULTS"
mkdir -p "$RESULTS/genus"
mkdir -p "$RESULTS/family"
mkdir -p "$RESULTS/tables"
mkdir -p "$LOG_DIR"

LOG_FILE="$LOG_DIR/14_cohort_differential_abundance.txt"

exec > >(tee "$LOG_FILE") 2>&1


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "PHASE 14 — COHORT DIFFERENTIAL ABUNDANCE"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "QIIME 2 environment:"
qiime --version

echo
echo "Model:"
echo "timepoint + participant"

echo
echo "Reference timepoint:"
echo "before"

echo
echo "Prevalence cutoff:"
echo "$PRV_CUT"

echo
echo "FDR method:"
echo "Benjamini-Hochberg"

echo
echo "Significance threshold:"
echo "$ALPHA"


# ============================================================
# 1. Verify inputs
# ============================================================

echo
echo "------------------------------------------------------------"
echo "INPUT VERIFICATION"
echo "------------------------------------------------------------"

qiime tools peek "$TABLE"
qiime tools peek "$TAXONOMY"

echo
echo "Metadata samples:"
tail -n +2 "$METADATA" | wc -l


# ============================================================
# 2. Collapse ASVs to genus level
# ============================================================

echo
echo "------------------------------------------------------------"
echo "COLLAPSING TO GENUS LEVEL"
echo "------------------------------------------------------------"

qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 6 \
  --o-collapsed-table "$RESULTS/genus/genus-table.qza"


# ============================================================
# 3. Collapse ASVs to family level
# ============================================================

echo
echo "------------------------------------------------------------"
echo "COLLAPSING TO FAMILY LEVEL"
echo "------------------------------------------------------------"

qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 5 \
  --o-collapsed-table "$RESULTS/family/family-table.qza"


# ============================================================
# 4. Genus-level ANCOM-BC
# ============================================================

echo
echo "------------------------------------------------------------"
echo "GENUS-LEVEL ANCOM-BC"
echo "------------------------------------------------------------"

qiime composition ancombc \
  --i-table "$RESULTS/genus/genus-table.qza" \
  --m-metadata-file "$METADATA" \
  --p-formula 'timepoint + participant' \
  --p-reference-levels timepoint::before \
  --p-p-adj-method BH \
  --p-prv-cut "$PRV_CUT" \
  --p-lib-cut 0 \
  --p-alpha "$ALPHA" \
  --o-differentials "$RESULTS/genus/ancombc.qza"


# ============================================================
# 5. Family-level ANCOM-BC
# ============================================================

echo
echo "------------------------------------------------------------"
echo "FAMILY-LEVEL ANCOM-BC"
echo "------------------------------------------------------------"

qiime composition ancombc \
  --i-table "$RESULTS/family/family-table.qza" \
  --m-metadata-file "$METADATA" \
  --p-formula 'timepoint + participant' \
  --p-reference-levels timepoint::before \
  --p-p-adj-method BH \
  --p-prv-cut "$PRV_CUT" \
  --p-lib-cut 0 \
  --p-alpha "$ALPHA" \
  --o-differentials "$RESULTS/family/ancombc.qza"


# ============================================================
# 6. Export collapsed count tables
# ============================================================

echo
echo "------------------------------------------------------------"
echo "EXPORTING COLLAPSED TABLES"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/genus/table-export"
rm -rf "$RESULTS/family/table-export"

qiime tools export \
  --input-path "$RESULTS/genus/genus-table.qza" \
  --output-path "$RESULTS/genus/table-export"

qiime tools export \
  --input-path "$RESULTS/family/family-table.qza" \
  --output-path "$RESULTS/family/table-export"


# ============================================================
# 7. Export ANCOM-BC results
# ============================================================

echo
echo "------------------------------------------------------------"
echo "EXPORTING ANCOM-BC RESULTS"
echo "------------------------------------------------------------"

rm -rf "$RESULTS/genus/ancombc-export"
rm -rf "$RESULTS/family/ancombc-export"

qiime tools export \
  --input-path "$RESULTS/genus/ancombc.qza" \
  --output-path "$RESULTS/genus/ancombc-export"

qiime tools export \
  --input-path "$RESULTS/family/ancombc.qza" \
  --output-path "$RESULTS/family/ancombc-export"


# ============================================================
# 8. Inspect exported result structure
# ============================================================

echo
echo "------------------------------------------------------------"
echo "ANCOM-BC EXPORT STRUCTURE"
echo "------------------------------------------------------------"

echo
echo "Genus:"
find "$RESULTS/genus/ancombc-export" \
  -maxdepth 2 \
  -type f \
  | sort

echo
echo "Family:"
find "$RESULTS/family/ancombc-export" \
  -maxdepth 2 \
  -type f \
  | sort


# ============================================================
# 9. Build taxon prevalence summaries
# ============================================================

echo
echo "------------------------------------------------------------"
echo "TAXONOMIC FEATURE SUMMARY"
echo "------------------------------------------------------------"

python - <<'PY'
from pathlib import Path
import pandas as pd


def biom_to_tsv(biom_file, output_file):

    import subprocess

    subprocess.run(
        [
            "biom",
            "convert",
            "-i",
            str(biom_file),
            "-o",
            str(output_file),
            "--to-tsv",
        ],
        check=True,
    )


for rank in ["genus", "family"]:

    biom_file = Path(
        f"results/phase14/{rank}/"
        "table-export/feature-table.biom"
    )

    tsv_file = Path(
        f"results/phase14/tables/"
        f"{rank}_counts.tsv"
    )

    biom_to_tsv(
        biom_file,
        tsv_file
    )

    df = pd.read_csv(
        tsv_file,
        sep="\t",
        comment="#",
        header=None,
    )

    print()
    print(
        f"{rank.upper()} TABLE"
    )

    print(
        f"Taxonomic features before "
        f"ANCOM-BC prevalence filtering: "
        f"{len(df)}"
    )
PY


# ============================================================
# 10. Find and preview ANCOM-BC tables
# ============================================================

echo
echo "------------------------------------------------------------"
echo "ANCOM-BC RESULT PREVIEW"
echo "------------------------------------------------------------"

python - <<'PY'
from pathlib import Path
import pandas as pd


for rank in [
    "genus",
    "family",
]:

    directory = Path(
        f"results/phase14/"
        f"{rank}/ancombc-export"
    )

    files = list(
        directory.rglob("*.tsv")
    )

    if not files:

        files = list(
            directory.rglob("*.csv")
        )

    print()
    print("=" * 70)
    print(rank.upper())
    print("=" * 70)

    if not files:

        print(
            "No TSV/CSV file found in export."
        )

        continue

    for file in files:

        print()
        print(f"File: {file}")

        try:

            if file.suffix == ".csv":

                df = pd.read_csv(
                    file
                )

            else:

                df = pd.read_csv(
                    file,
                    sep="\t"
                )

            print(
                f"Rows: {len(df)}"
            )

            print(
                "Columns:"
            )

            print(
                list(df.columns)
            )

            print()

            print(
                df.head(10)
                .to_string(
                    index=False
                )
            )

        except Exception as error:

            print(
                f"Could not parse: {error}"
            )

PY


# ============================================================
# 11. Output inventory
# ============================================================

echo
echo "------------------------------------------------------------"
echo "PHASE 14 OUTPUTS"
echo "------------------------------------------------------------"

find "$RESULTS" \
  -maxdepth 4 \
  -type f \
  | sort


# ============================================================
# 12. Conclusion
# ============================================================

echo
echo "============================================================"
echo "PHASE 14 CONCLUSION"
echo "============================================================"

echo
echo "Cohort-wide differential abundance was evaluated using"
echo "ANCOM-BC at genus and family taxonomic levels."
echo
echo "The model included both timepoint and participant:"
echo
echo "    abundance ~ timepoint + participant"
echo
echo "Pre-FMT samples were used as the reference condition."
echo
echo "Participant was included as a blocking covariate to account"
echo "for persistent inter-individual differences in microbial"
echo "composition."
echo
echo "Taxa present in fewer than 10% of samples were excluded"
echo "using the ANCOM-BC prevalence filter."
echo
echo "Benjamini-Hochberg adjustment was used to control the"
echo "false-discovery rate."
echo
echo "Genus-level results represent the primary differential"
echo "abundance analysis; family-level results provide a"
echo "secondary taxonomic summary."
echo
echo "============================================================"
echo "END PHASE 14"
echo "============================================================"
