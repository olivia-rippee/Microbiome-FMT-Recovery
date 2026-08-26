#!/usr/bin/env bash

# ============================================================
# PHASE 9 — DOWNLOAD AND VERIFY COHORT FASTQ FILES
# ============================================================
#
# Purpose:
#   Download paired-end FASTQ files for the validated
#   longitudinal FMT cohort using NCBI SRA Toolkit.
#
# Cohort:
#   52 participants
#   4 timepoints per participant
#   208 paired-end samples
#
# Input:
#   data/cohort/PRJNA298590_cohort_download_manifest.tsv
#
# Outputs:
#   data/raw/cohort/
#   data/cohort/fastq_manifest.tsv
#   data/cohort/samples.tsv
#   data/cohort/fastq_verification.tsv
#   logs/09_download_cohort_fastqs.txt
#
# Important:
#   - Restartable.
#   - Existing valid paired FASTQs are skipped.
#   - Missing samples are retrieved with prefetch.
#   - FASTQs are generated with fasterq-dump --split-files.
#   - pigz is used for compression.
#   - Final gzip integrity and mate-pair presence are validated.
#
# ============================================================

set -uo pipefail

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SOURCE_MANIFEST="data/cohort/PRJNA298590_cohort_download_manifest.tsv"

FASTQ_DIR="data/raw/cohort"
SRA_DIR="data/sra/cohort"
WORK_DIR="work/phase9"

QIIME_MANIFEST="data/cohort/fastq_manifest.tsv"
SAMPLE_METADATA="data/cohort/samples.tsv"
VERIFICATION="data/cohort/fastq_verification.tsv"

LOG_DIR="logs"
LOG_FILE="$LOG_DIR/09_download_cohort_fastqs.txt"

THREADS=4

mkdir -p "$FASTQ_DIR"
mkdir -p "$SRA_DIR"
mkdir -p "$WORK_DIR"
mkdir -p "$LOG_DIR"

exec > >(tee "$LOG_FILE") 2>&1


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

echo "============================================================"
echo "PHASE 9 — COHORT FASTQ DOWNLOAD AND VERIFICATION"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "Project:"
pwd

echo
echo "SRA Toolkit:"
prefetch --version
fasterq-dump --version

echo
echo "Compression:"
pigz --version | head -n 1


# ------------------------------------------------------------
# Validate input
# ------------------------------------------------------------

if [ ! -f "$SOURCE_MANIFEST" ]; then

    echo
    echo "ERROR: Missing cohort manifest:"
    echo "$SOURCE_MANIFEST"

    exit 1

fi


# ------------------------------------------------------------
# Cohort summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "COHORT SUMMARY"
echo "------------------------------------------------------------"

python - <<'PY'
import pandas as pd

f = "data/cohort/PRJNA298590_cohort_download_manifest.tsv"

df = pd.read_csv(
    f,
    sep="\t",
)

print(f"Samples: {len(df)}")
print(f"Participants: {df['participant'].nunique()}")

print("\nSamples per timepoint:")
print(
    df["timepoint"]
    .value_counts()
    .sort_index()
    .to_string()
)

if "run_total_bases" in df.columns:

    total = df["run_total_bases"].sum()

    print(f"\nTotal bases: {total:,}")
    print(f"Approx sequence volume: {total / 1e9:.2f} Gbases")
PY


# ------------------------------------------------------------
# Download / convert samples
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "RETRIEVING PAIRED FASTQ FILES"
echo "------------------------------------------------------------"

TOTAL=0
COMPLETE=0
FAILED=0

while IFS=$'\t' read -r SAMPLE PARTICIPANT TIMEPOINT RUN
do

    TOTAL=$((TOTAL + 1))

    R1="$FASTQ_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    SAMPLE_WORK="$WORK_DIR/$SAMPLE"

    echo
    echo "============================================================"
    echo "Sample $TOTAL: $SAMPLE"
    echo "Participant: $PARTICIPANT"
    echo "Timepoint:   $TIMEPOINT"
    echo "Run:         $RUN"
    echo "============================================================"

    # --------------------------------------------------------
    # Skip already-complete sample
    # --------------------------------------------------------

    if \
        [ -s "$R1" ] && \
        [ -s "$R2" ] && \
        gzip -t "$R1" 2>/dev/null && \
        gzip -t "$R2" 2>/dev/null
    then

        echo "Existing paired FASTQs are valid."
        echo "Skipping sample."

        COMPLETE=$((COMPLETE + 1))

        continue

    fi

    # --------------------------------------------------------
    # Remove partial/corrupt outputs
    # --------------------------------------------------------

    if [ -e "$R1" ] && ! gzip -t "$R1" 2>/dev/null; then

        echo "Removing incomplete/corrupt R1."
        rm -f "$R1"

    fi

    if [ -e "$R2" ] && ! gzip -t "$R2" 2>/dev/null; then

        echo "Removing incomplete/corrupt R2."
        rm -f "$R2"

    fi

    rm -rf "$SAMPLE_WORK"
    mkdir -p "$SAMPLE_WORK"

    SAMPLE_FAILED=0

    # --------------------------------------------------------
    # Prefetch SRA
    # --------------------------------------------------------

    echo
    echo "Prefetching $RUN..."

    prefetch \
        "$RUN" \
        --output-directory "$SRA_DIR"

    if [ $? -ne 0 ]; then

        echo "ERROR: prefetch failed for $RUN"

        SAMPLE_FAILED=1

    fi

    SRA_PATH="$SRA_DIR/$RUN/$RUN.sra"

    # Some prefetch configurations may not use .sra suffix
    if [ ! -f "$SRA_PATH" ]; then

        ALT_PATH="$SRA_DIR/$RUN/$RUN"

        if [ -f "$ALT_PATH" ]; then
            SRA_PATH="$ALT_PATH"
        fi

    fi

    if [ "$SAMPLE_FAILED" -eq 0 ] && [ ! -f "$SRA_PATH" ]; then

        echo "ERROR: Could not locate prefetched SRA file:"
        echo "$SRA_PATH"

        SAMPLE_FAILED=1

    fi

    # --------------------------------------------------------
    # Convert SRA to paired FASTQ
    # --------------------------------------------------------

    if [ "$SAMPLE_FAILED" -eq 0 ]; then

        echo
        echo "Running fasterq-dump..."

        fasterq-dump \
            "$SRA_PATH" \
            --split-files \
            --threads "$THREADS" \
            --temp "$SAMPLE_WORK" \
            --outdir "$SAMPLE_WORK"

        if [ $? -ne 0 ]; then

            echo "ERROR: fasterq-dump failed for $RUN"

            SAMPLE_FAILED=1

        fi

    fi

    RAW_R1="$SAMPLE_WORK/${RUN}_1.fastq"
    RAW_R2="$SAMPLE_WORK/${RUN}_2.fastq"

    # --------------------------------------------------------
    # Validate paired uncompressed outputs
    # --------------------------------------------------------

    if [ "$SAMPLE_FAILED" -eq 0 ]; then

        if [ ! -s "$RAW_R1" ]; then

            echo "ERROR: Missing forward FASTQ:"
            echo "$RAW_R1"

            SAMPLE_FAILED=1

        fi

        if [ ! -s "$RAW_R2" ]; then

            echo "ERROR: Missing reverse FASTQ:"
            echo "$RAW_R2"

            SAMPLE_FAILED=1

        fi

    fi

    # --------------------------------------------------------
    # Compress and move outputs
    # --------------------------------------------------------

    if [ "$SAMPLE_FAILED" -eq 0 ]; then

        echo
        echo "Compressing FASTQs with pigz..."

        pigz \
            -p "$THREADS" \
            "$RAW_R1"

        pigz \
            -p "$THREADS" \
            "$RAW_R2"

        if \
            [ ! -s "${RAW_R1}.gz" ] || \
            [ ! -s "${RAW_R2}.gz" ]
        then

            echo "ERROR: Compression failed."

            SAMPLE_FAILED=1

        fi

    fi

    if [ "$SAMPLE_FAILED" -eq 0 ]; then

        mv \
            "${RAW_R1}.gz" \
            "$R1"

        mv \
            "${RAW_R2}.gz" \
            "$R2"

    fi

    # --------------------------------------------------------
    # Final gzip validation
    # --------------------------------------------------------

    if [ "$SAMPLE_FAILED" -eq 0 ]; then

        if ! gzip -t "$R1" 2>/dev/null; then

            echo "ERROR: Final R1 gzip validation failed."

            SAMPLE_FAILED=1

        fi

        if ! gzip -t "$R2" 2>/dev/null; then

            echo "ERROR: Final R2 gzip validation failed."

            SAMPLE_FAILED=1

        fi

    fi

    # --------------------------------------------------------
    # Result
    # --------------------------------------------------------

    if [ "$SAMPLE_FAILED" -eq 1 ]; then

        echo
        echo "SAMPLE STATUS: FAILED"

        FAILED=$((FAILED + 1))

    else

        echo
        echo "SAMPLE STATUS: COMPLETE"

        COMPLETE=$((COMPLETE + 1))

        rm -rf "$SAMPLE_WORK"

    fi

done < <(
python - <<'PY'
import pandas as pd

f = "data/cohort/PRJNA298590_cohort_download_manifest.tsv"

df = pd.read_csv(
    f,
    sep="\t",
)

required = [
    "sample-id",
    "participant",
    "timepoint",
    "run_accession",
]

missing = [
    column
    for column in required
    if column not in df.columns
]

if missing:

    raise SystemExit(
        "Missing required columns: "
        + ", ".join(missing)
    )

for _, row in df.iterrows():

    print(
        "\t".join(
            [
                str(row["sample-id"]),
                str(row["participant"]),
                str(row["timepoint"]),
                str(row["run_accession"]),
            ]
        )
    )
PY
)


# ------------------------------------------------------------
# Verify entire cohort
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "VERIFYING ALL FASTQ PAIRS"
echo "------------------------------------------------------------"

python - <<'PY'
from pathlib import Path
import gzip
import pandas as pd


source = (
    "data/cohort/"
    "PRJNA298590_cohort_download_manifest.tsv"
)

fastq_dir = Path(
    "data/raw/cohort"
)

verification_file = (
    "data/cohort/fastq_verification.tsv"
)

qiime_manifest = (
    "data/cohort/fastq_manifest.tsv"
)

metadata_file = (
    "data/cohort/samples.tsv"
)


df = pd.read_csv(
    source,
    sep="\t",
)

records = []


def gzip_valid(path):

    if not path.exists():
        return False

    try:

        with gzip.open(
            path,
            "rb",
        ) as handle:

            while handle.read(
                1024 * 1024
            ):
                pass

        return True

    except Exception:

        return False


for _, row in df.iterrows():

    sample = row["sample-id"]

    r1 = (
        fastq_dir /
        f"{sample}_R1.fastq.gz"
    )

    r2 = (
        fastq_dir /
        f"{sample}_R2.fastq.gz"
    )

    r1_valid = gzip_valid(r1)
    r2_valid = gzip_valid(r2)

    records.append(
        {
            "sample-id": sample,
            "participant": row["participant"],
            "timepoint": row["timepoint"],
            "run_accession": row["run_accession"],
            "r1_exists": r1.exists(),
            "r2_exists": r2.exists(),
            "r1_size_bytes": (
                r1.stat().st_size
                if r1.exists()
                else 0
            ),
            "r2_size_bytes": (
                r2.stat().st_size
                if r2.exists()
                else 0
            ),
            "r1_gzip_valid": r1_valid,
            "r2_gzip_valid": r2_valid,
            "pair_valid": (
                r1_valid
                and r2_valid
            ),
        }
    )


verification = pd.DataFrame(
    records
)

verification.to_csv(
    verification_file,
    sep="\t",
    index=False,
)


valid_count = int(
    verification[
        "pair_valid"
    ].sum()
)

invalid_count = (
    len(verification)
    - valid_count
)

print(
    f"Samples checked: {len(verification)}"
)

print(
    f"Valid paired samples: {valid_count}"
)

print(
    f"Invalid/incomplete samples: {invalid_count}"
)


if invalid_count:

    print(
        "\nInvalid samples:"
    )

    print(
        verification[
            ~verification[
                "pair_valid"
            ]
        ][
            [
                "sample-id",
                "run_accession",
                "r1_exists",
                "r2_exists",
                "r1_gzip_valid",
                "r2_gzip_valid",
            ]
        ]
        .to_string(
            index=False
        )
    )

    raise SystemExit(
        "\nERROR: Cohort FASTQ verification failed."
    )


# ------------------------------------------------------------
# Build QIIME manifest
# ------------------------------------------------------------

project = Path.cwd()

manifest = pd.DataFrame(
    {
        "sample-id": df["sample-id"],

        "forward-absolute-filepath": [
            str(
                (
                    project
                    / fastq_dir
                    / f"{sample}_R1.fastq.gz"
                ).resolve()
            )
            for sample
            in df["sample-id"]
        ],

        "reverse-absolute-filepath": [
            str(
                (
                    project
                    / fastq_dir
                    / f"{sample}_R2.fastq.gz"
                ).resolve()
            )
            for sample
            in df["sample-id"]
        ],
    }
)

manifest.to_csv(
    qiime_manifest,
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Build sample metadata
# ------------------------------------------------------------

metadata = df[
    [
        "sample-id",
        "participant",
        "timepoint",
    ]
].copy()

days = {
    "before": 0,
    "7d": 7,
    "14d": 14,
    "30d": 30,
}

metadata[
    "days_after_fmt"
] = (
    metadata[
        "timepoint"
    ]
    .map(
        days
    )
)

metadata.to_csv(
    metadata_file,
    sep="\t",
    index=False,
)


print(
    "\nQIIME manifest:"
)

print(
    qiime_manifest
)

print(
    "\nSample metadata:"
)

print(
    metadata_file
)
PY

VERIFY_STATUS=$?


# ------------------------------------------------------------
# Disk summary
# ------------------------------------------------------------

echo
echo "------------------------------------------------------------"
echo "FASTQ STORAGE SUMMARY"
echo "------------------------------------------------------------"

du -sh "$FASTQ_DIR"

echo
echo "Compressed FASTQ files:"
find "$FASTQ_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.fastq.gz' \
    | wc -l


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "============================================================"
echo "PHASE 9 SUMMARY"
echo "============================================================"

echo
echo "Samples expected:"
echo "$TOTAL"

echo
echo "Samples completed/skipped:"
echo "$COMPLETE"

echo
echo "Samples with retrieval/conversion failures:"
echo "$FAILED"


if [ "$VERIFY_STATUS" -ne 0 ]; then

    echo
    echo "FINAL STATUS: INCOMPLETE"
    echo
    echo "One or more samples are missing or invalid."
    echo "Rerun this script to retry incomplete samples."

    echo
    echo "============================================================"

    exit 1

fi


echo
echo "FINAL STATUS: COMPLETE"

echo
echo "All cohort samples contain valid paired-end FASTQs."

echo
echo "The QIIME 2 paired-end manifest and longitudinal"
echo "sample metadata table were generated successfully."

echo
echo "============================================================"
echo "END PHASE 9"
echo "============================================================"
