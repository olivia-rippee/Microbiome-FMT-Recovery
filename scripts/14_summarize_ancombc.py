from pathlib import Path
import pandas as pd

BASE = Path("results/phase14")
OUT = BASE / "tables"
OUT.mkdir(parents=True, exist_ok=True)

TIMEPOINTS = {
    "timepoint7d":  "7d_vs_before",
    "timepoint14d": "14d_vs_before",
    "timepoint30d": "30d_vs_before",
}

for rank in ["genus", "family"]:

    d = BASE / rank / "ancombc-export"

    lfc = pd.read_csv(d / "lfc_slice.csv")
    pval = pd.read_csv(d / "p_val_slice.csv")
    qval = pd.read_csv(d / "q_val_slice.csv")

    print()
    print("=" * 80)
    print(rank.upper())
    print("=" * 80)

    all_results = []

    for coef, comparison in TIMEPOINTS.items():

        if coef not in lfc.columns:
            raise RuntimeError(
                f"{coef} not found in {rank} ANCOM-BC output.\n"
                f"Available columns: {lfc.columns.tolist()}"
            )

        x = pd.DataFrame({
            "taxon": lfc["id"],
            "comparison": comparison,
            "lfc": pd.to_numeric(lfc[coef], errors="coerce"),
            "p_value": pd.to_numeric(pval[coef], errors="coerce"),
            "q_value": pd.to_numeric(qval[coef], errors="coerce"),
        })

        x["significant"] = x["q_value"] < 0.05

        x["direction"] = x["lfc"].apply(
            lambda z:
                "increased" if z > 0
                else "decreased" if z < 0
                else "no_change"
        )

        x = x.sort_values(
            ["q_value", "p_value"],
            na_position="last"
        )

        all_results.append(x)

        sig = x[x["significant"]].copy()

        print()
        print("-" * 80)
        print(comparison)
        print("-" * 80)

        print(f"Taxa tested:       {len(x)}")
        print(f"Significant taxa:  {len(sig)}")

        if len(sig):

            print()
            print(
                sig[
                    [
                        "taxon",
                        "lfc",
                        "p_value",
                        "q_value",
                        "direction",
                    ]
                ].to_string(index=False)
            )

        else:
            print("No taxa significant at BH q < 0.05.")

    combined = pd.concat(
        all_results,
        ignore_index=True
    )

    combined.to_csv(
        OUT / f"{rank}_ancombc_all_results.tsv",
        sep="\t",
        index=False
    )

    significant = combined[
        combined["significant"]
    ].copy()

    significant.to_csv(
        OUT / f"{rank}_ancombc_significant.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Taxa significant at multiple post-FMT timepoints
    # --------------------------------------------------------

    if len(significant):

        persistence = (
            significant
            .groupby("taxon")
            .agg(
                significant_timepoints=("comparison", "nunique"),
                comparisons=(
                    "comparison",
                    lambda x: ";".join(sorted(x))
                ),
                directions=(
                    "direction",
                    lambda x: ";".join(x)
                ),
                min_q=("q_value", "min"),
            )
            .reset_index()
            .sort_values(
                ["significant_timepoints", "min_q"],
                ascending=[False, True]
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

    persistence.to_csv(
        OUT / f"{rank}_ancombc_persistence.tsv",
        sep="\t",
        index=False
    )

    print()
    print("-" * 80)
    print("PERSISTENT DIFFERENTIAL TAXA")
    print("-" * 80)

    persistent = persistence[
        persistence["significant_timepoints"] >= 2
    ]

    if len(persistent):
        print(persistent.to_string(index=False))
    else:
        print("No taxa significant at two or more timepoints.")


print()
print("=" * 80)
print("OUTPUTS")
print("=" * 80)

for f in sorted(OUT.glob("*ancombc*.tsv")):
    print(f)
