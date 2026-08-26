#!/usr/bin/env python3

"""
PHASE 14 — SUMMARIZE ANCOM-BC RESULTS

Purpose
-------
Convert exported QIIME 2 ANCOM-BC coefficient, p-value,
and q-value slices into compact interpretable tables.

For genus and family levels, this script generates:

  *_ancombc_all_results.tsv
  *_ancombc_significant.tsv
  *_ancombc_persistence.tsv

Comparisons
-----------
7d vs before
14d vs before
30d vs before

Significance
------------
Benjamini-Hochberg adjusted q < 0.05

Interpretation
--------------
Positive LFC:
    higher abundance post-FMT relative to before

Negative LFC:
    lower abundance post-FMT relative to before
"""

from pathlib import Path

import pandas as pd


# ============================================================
# Configuration
# ============================================================

BASE = Path("results/phase14")
OUT = BASE / "tables"

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

TIMEPOINTS = {
    "timepoint7d": "7d_vs_before",
    "timepoint14d": "14d_vs_before",
    "timepoint30d": "30d_vs_before",
}

ALPHA = 0.05


# ============================================================
# Process genus and family levels
# ============================================================

print("=" * 80)
print("PHASE 14 — ANCOM-BC RESULT SUMMARY")
print("=" * 80)


for rank in [
    "genus",
    "family",
]:

    directory = (
        BASE
        / rank
        / "ancombc-export"
    )

    lfc_file = (
        directory
        / "lfc_slice.csv"
    )

    p_file = (
        directory
        / "p_val_slice.csv"
    )

    q_file = (
        directory
        / "q_val_slice.csv"
    )

    for path in [
        lfc_file,
        p_file,
        q_file,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                f"Required ANCOM-BC export missing: {path}"
            )

    lfc = pd.read_csv(
        lfc_file
    )

    pval = pd.read_csv(
        p_file
    )

    qval = pd.read_csv(
        q_file
    )

    print()
    print("=" * 80)
    print(rank.upper())
    print("=" * 80)

    all_results = []

    # --------------------------------------------------------
    # Extract each timepoint coefficient
    # --------------------------------------------------------

    for coefficient, comparison in TIMEPOINTS.items():

        for table_name, table in [
            ("LFC", lfc),
            ("p-value", pval),
            ("q-value", qval),
        ]:

            if coefficient not in table.columns:

                raise RuntimeError(
                    f"{coefficient} not found in "
                    f"{rank} {table_name} output.\n"
                    f"Available columns: "
                    f"{table.columns.tolist()}"
                )

        result = pd.DataFrame(
            {
                "taxon": lfc["id"],
                "comparison": comparison,
                "lfc": pd.to_numeric(
                    lfc[coefficient],
                    errors="coerce",
                ),
                "p_value": pd.to_numeric(
                    pval[coefficient],
                    errors="coerce",
                ),
                "q_value": pd.to_numeric(
                    qval[coefficient],
                    errors="coerce",
                ),
            }
        )

        result["significant"] = (
            result["q_value"]
            < ALPHA
        )

        result["direction"] = (
            result["lfc"]
            .apply(
                lambda value:
                    "increased"
                    if value > 0
                    else "decreased"
                    if value < 0
                    else "no_change"
            )
        )

        result = result.sort_values(
            [
                "q_value",
                "p_value",
            ],
            na_position="last",
        )

        all_results.append(
            result
        )

        significant = result[
            result["significant"]
        ].copy()

        print()
        print("-" * 80)
        print(comparison)
        print("-" * 80)

        print(
            f"Taxa tested:      {len(result)}"
        )

        print(
            f"Significant taxa: {len(significant)}"
        )

        if len(significant):

            print()

            print(
                significant[
                    [
                        "taxon",
                        "lfc",
                        "p_value",
                        "q_value",
                        "direction",
                    ]
                ]
                .to_string(
                    index=False
                )
            )

        else:

            print(
                "No taxa significant at BH q < 0.05."
            )

    # --------------------------------------------------------
    # Combine comparisons
    # --------------------------------------------------------

    combined = pd.concat(
        all_results,
        ignore_index=True,
    )

    all_output = (
        OUT
        / f"{rank}_ancombc_all_results.tsv"
    )

    combined.to_csv(
        all_output,
        sep="\t",
        index=False,
    )

    significant = combined[
        combined["significant"]
    ].copy()

    significant_output = (
        OUT
        / f"{rank}_ancombc_significant.tsv"
    )

    significant.to_csv(
        significant_output,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Persistence summary
    # --------------------------------------------------------

    if len(significant):

        persistence = (
            significant
            .groupby("taxon")
            .agg(
                significant_timepoints=(
                    "comparison",
                    "nunique",
                ),
                comparisons=(
                    "comparison",
                    lambda values:
                        ";".join(
                            sorted(values)
                        ),
                ),
                directions=(
                    "direction",
                    lambda values:
                        ";".join(values),
                ),
                min_q=(
                    "q_value",
                    "min",
                ),
            )
            .reset_index()
            .sort_values(
                [
                    "significant_timepoints",
                    "min_q",
                ],
                ascending=[
                    False,
                    True,
                ],
            )
        )

    else:

        persistence = pd.DataFrame(
            columns=[
                "taxon",
                "significant_timepoints",
                "comparisons",
                "directions",
                "min_q",
            ]
        )

    persistence_output = (
        OUT
        / f"{rank}_ancombc_persistence.tsv"
    )

    persistence.to_csv(
        persistence_output,
        sep="\t",
        index=False,
    )

    print()
    print("-" * 80)
    print("PERSISTENT DIFFERENTIAL TAXA")
    print("-" * 80)

    persistent = persistence[
        persistence[
            "significant_timepoints"
        ] >= 2
    ]

    if len(persistent):

        print(
            persistent.to_string(
                index=False
            )
        )

    else:

        print(
            "No taxa significant at two "
            "or more timepoints."
        )


# ============================================================
# Outputs
# ============================================================

print()
print("=" * 80)
print("OUTPUTS")
print("=" * 80)

for path in sorted(
    OUT.glob("*ancombc*.tsv")
):

    print(path)


print()
print("=" * 80)
print("END PHASE 14 ANCOM-BC SUMMARY")
print("=" * 80)
