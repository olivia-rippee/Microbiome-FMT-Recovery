#!/usr/bin/env python3

"""
PHASE 8 — BUILD EXPANDED LONGITUDINAL COHORT METADATA

Purpose
-------
Parse PRJNA298590 run metadata, identify recipient IDs and
FMT timepoints, validate complete longitudinal participants,
and generate clean cohort files for downstream downloading
and multi-participant analysis.

Required recipient timepoints
-----------------------------
before
7d
14d
30d

Important
---------
The source metadata uses several inconsistent naming styles,
including examples such as:

DuPont.R.056.14.days.after.FMT
DuPontR002914day
DuPontR003day14
DuPontR0047day
DuPontR005BeforeFMT

The parser identifies the timepoint first, removes the
timepoint portion from the sample name, and then extracts
the recipient ID.

This prevents compact names such as:

DuPontR002914day

from being incorrectly interpreted as participant R002914.

Donor records beginning with DuPontD or DuPont.D are not
included in the recipient longitudinal cohort. They are
retained in the unparsed/manual-review output for possible
future donor-recipient analyses.
"""

from pathlib import Path
import re
import pandas as pd


# ============================================================
# Configuration
# ============================================================

INPUT = Path(
    "data/cohort/PRJNA298590_metadata.tsv"
)

OUTPUT_DIR = Path(
    "data/cohort"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

PARSED_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_parsed_metadata.tsv"
)

INVENTORY_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_participant_inventory.tsv"
)

COMPLETE_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_complete_participants.tsv"
)

COHORT_RUNS_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_complete_cohort_runs.tsv"
)

UNPARSED_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_unparsed_runs.tsv"
)

DUPLICATE_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_duplicate_timepoints.tsv"
)

MANIFEST_OUTPUT = (
    OUTPUT_DIR /
    "PRJNA298590_cohort_download_manifest.tsv"
)

TIMEPOINTS = [
    "before",
    "7d",
    "14d",
    "30d",
]

TIMEPOINT_ORDER = {
    "before": 0,
    "7d": 1,
    "14d": 2,
    "30d": 3,
}


# ============================================================
# Helper functions
# ============================================================

def normalize_text(value):
    """
    Convert metadata text into a parser-friendly string.
    """

    if pd.isna(value):
        return ""

    return str(value).strip()


def is_donor_record(text):
    """
    Identify obvious donor samples.

    Examples:
      DuPontD00490413
      DuPont.D.012.FRE010714
    """

    text = normalize_text(text)

    return bool(
        re.search(
            r"DuPont\.?D",
            text,
            flags=re.I,
        )
    )


def parse_timepoint(text):
    """
    Identify the FMT timepoint from heterogeneous recipient names.

    Handles examples such as:

      DuPontR0097day
      DuPontR002914day
      DuPontR002930day
      DuPontR003day14
      DuPont.R.056.14.days.after.FMT
      DuPontR005BeforeFMT
    """

    raw = normalize_text(text)
    lower = raw.lower()

    # --------------------------------------------------------
    # Pre-FMT / baseline
    # --------------------------------------------------------

    before_patterns = [
        r"before[\._\-\s]*fmt",
        r"pre[\._\-\s]*fmt",
        r"baseline",
    ]

    for pattern in before_patterns:

        if re.search(
            pattern,
            lower,
        ):
            return "before"

    # --------------------------------------------------------
    # Post-FMT timepoints
    #
    # Check 30 and 14 before 7 to avoid accidental partial
    # matching.
    # --------------------------------------------------------

    patterns_30 = [
        r"30[\._\-\s]*days?",
        r"30day",
        r"day30",
        r"days?[\._\-\s]*30",
    ]

    for pattern in patterns_30:

        if re.search(
            pattern,
            lower,
        ):
            return "30d"

    patterns_14 = [
        r"14[\._\-\s]*days?",
        r"14day",
        r"day14",
        r"days?[\._\-\s]*14",
    ]

    for pattern in patterns_14:

        if re.search(
            pattern,
            lower,
        ):
            return "14d"

    patterns_7 = [
        r"7[\._\-\s]*days?",
        r"7day",
        r"day7",
        r"days?[\._\-\s]*7",
    ]

    for pattern in patterns_7:

        if re.search(
            pattern,
            lower,
        ):
            return "7d"

    return None


def strip_timepoint(text, timepoint):
    """
    Remove the timepoint portion before extracting participant ID.

    Examples:

      DuPontR00297day
        -> DuPontR0029

      DuPontR0097day
        -> DuPontR009

      DuPontR00107day
        -> DuPontR0010

      DuPontR003day14
        -> DuPontR003

      DuPont.R.056.14.days.after.FMT
        -> DuPont.R.056.
    """

    value = normalize_text(text)

    # --------------------------------------------------------
    # Remove descriptive post-FMT wording first
    # --------------------------------------------------------

    value = re.sub(
        r"after[\._\-\s]*fmt",
        "",
        value,
        flags=re.I,
    )

    # --------------------------------------------------------
    # Remove timepoint expression
    # --------------------------------------------------------

    if timepoint == "before":

        patterns = [
            r"before[\._\-\s]*fmt",
            r"pre[\._\-\s]*fmt",
            r"baseline",
        ]

    elif timepoint == "7d":

        patterns = [
            r"7[\._\-\s]*days?",
            r"7day",
            r"day7",
            r"days?[\._\-\s]*7",
        ]

    elif timepoint == "14d":

        patterns = [
            r"14[\._\-\s]*days?",
            r"14day",
            r"day14",
            r"days?[\._\-\s]*14",
        ]

    elif timepoint == "30d":

        patterns = [
            r"30[\._\-\s]*days?",
            r"30day",
            r"day30",
            r"days?[\._\-\s]*30",
        ]

    else:

        return value

    cleaned = value

    for pattern in patterns:

        cleaned = re.sub(
            pattern,
            "",
            cleaned,
            flags=re.I,
        )

    cleaned = re.sub(
        r"fmt",
        "",
        cleaned,
        flags=re.I,
    )

    return cleaned


def normalize_participant_number(number):
    """
    Normalize recipient number to R### where practical.

    Examples:

      3     -> R003
      9     -> R009
      29    -> R029
      056   -> R056
    """

    number = str(number)

    number = (
        number.lstrip("0")
        or "0"
    )

    return f"R{int(number):03d}"


def parse_participant(text, timepoint):
    """
    Extract the participant after timepoint removal.

    Handles:

      DuPont.R.056.14.days.after.FMT
      DuPontR002914day
      DuPontR003day14
      DuPontR0047day
      DuPontR005BeforeFMT
    """

    if timepoint is None:
        return None

    if is_donor_record(text):
        return None

    cleaned = strip_timepoint(
        text,
        timepoint,
    )

    patterns = [
        r"DuPont\.R\.(\d+)",
        r"DuPontR(\d+)",
        r"\bR[\._\-]?(\d+)\b",
    ]

    for pattern in patterns:

        match = re.search(
            pattern,
            cleaned,
            flags=re.I,
        )

        if match:

            return normalize_participant_number(
                match.group(1)
            )

    return None


def parse_row(row):
    """
    Parse one metadata row.

    Strategy:
      1. Look across several metadata fields.
      2. Identify a timepoint first.
      3. Remove that timepoint.
      4. Extract recipient ID.
    """

    candidate_fields = [
        "experiment_title",
        "library_name",
        "run_alias",
        "experiment_alias",
        "file_prefix",
    ]

    # --------------------------------------------------------
    # Explicit donor exclusion
    # --------------------------------------------------------

    for field in candidate_fields:

        if field not in row:
            continue

        value = normalize_text(
            row[field]
        )

        if value and is_donor_record(value):

            return (
                None,
                None,
                "donor_record",
            )

    # --------------------------------------------------------
    # Determine timepoint
    # --------------------------------------------------------

    timepoint = None

    for field in candidate_fields:

        if field not in row:
            continue

        value = normalize_text(
            row[field]
        )

        if not value:
            continue

        tp = parse_timepoint(
            value
        )

        if tp is not None:

            timepoint = tp
            break

    if timepoint is None:

        return (
            None,
            None,
            "timepoint_unparsed",
        )

    # --------------------------------------------------------
    # Determine participant
    # --------------------------------------------------------

    for field in candidate_fields:

        if field not in row:
            continue

        value = normalize_text(
            row[field]
        )

        if not value:
            continue

        participant = parse_participant(
            value,
            timepoint,
        )

        if participant is not None:

            return (
                participant,
                timepoint,
                field,
            )

    return (
        None,
        timepoint,
        "participant_unparsed",
    )


# ============================================================
# Load metadata
# ============================================================

print("=" * 72)
print("PHASE 8 — BUILD LONGITUDINAL COHORT")
print("=" * 72)

print(
    f"\nReading:\n  {INPUT}"
)

df = pd.read_csv(
    INPUT,
    sep="\t",
)

print(
    f"\nTotal SRA runs: {len(df)}"
)

if "biosample" in df.columns:

    print(
        "Unique BioSamples: "
        f"{df['biosample'].nunique()}"
    )


# ============================================================
# Restrict to paired-end sequencing
# ============================================================

if "library_layout" in df.columns:

    paired = df[
        df["library_layout"]
        .astype(str)
        .str.upper()
        .eq("PAIRED")
    ].copy()

else:

    paired = df.copy()

print(
    "Paired-end runs retained: "
    f"{len(paired)}"
)


# ============================================================
# Parse participant and timepoint
# ============================================================

parsed_values = paired.apply(
    parse_row,
    axis=1,
    result_type="expand",
)

parsed_values.columns = [
    "participant",
    "timepoint",
    "parse_source",
]

paired = pd.concat(
    [
        paired.reset_index(
            drop=True
        ),
        parsed_values.reset_index(
            drop=True
        ),
    ],
    axis=1,
)


# ============================================================
# Separate recipient and unparsed/non-recipient runs
# ============================================================

parsed = paired[
    paired["participant"].notna()
    & paired["timepoint"].notna()
].copy()

unparsed = paired[
    paired["participant"].isna()
    | paired["timepoint"].isna()
].copy()

print(
    "\nRuns with participant + recognized timepoint: "
    f"{len(parsed)}"
)

print(
    "Participants detected: "
    f"{parsed['participant'].nunique()}"
)

print(
    "Runs not assigned to the recipient cohort: "
    f"{len(unparsed)}"
)


# ============================================================
# Donor summary
# ============================================================

donor_runs = unparsed[
    unparsed["parse_source"]
    .eq("donor_record")
].copy()

print(
    "Donor/non-recipient records detected: "
    f"{len(donor_runs)}"
)


# ============================================================
# Detect duplicate participant/timepoint combinations
# ============================================================

duplicate_counts = (
    parsed
    .groupby(
        [
            "participant",
            "timepoint",
        ]
    )
    .size()
    .reset_index(
        name="n_runs"
    )
)

duplicates = duplicate_counts[
    duplicate_counts["n_runs"] > 1
].copy()

print(
    "\nDuplicate participant/timepoint combinations: "
    f"{len(duplicates)}"
)


# ============================================================
# Build participant inventory
# ============================================================

inventory = (
    parsed
    .groupby(
        [
            "participant",
            "timepoint",
        ]
    )
    .size()
    .unstack(
        fill_value=0
    )
)

for timepoint in TIMEPOINTS:

    if timepoint not in inventory.columns:

        inventory[
            timepoint
        ] = 0

inventory = inventory[
    TIMEPOINTS
]

inventory["n_timepoints"] = (
    (
        inventory[
            TIMEPOINTS
        ] > 0
    )
    .sum(
        axis=1
    )
)

inventory["exactly_one_each"] = (
    (
        inventory[
            TIMEPOINTS
        ] == 1
    )
    .all(
        axis=1
    )
)

inventory = inventory.sort_index()


# ============================================================
# Complete longitudinal cohort
# ============================================================

complete = inventory[
    inventory[
        "exactly_one_each"
    ]
].copy()

complete_ids = (
    complete.index.tolist()
)

print(
    "\n"
    + "=" * 72
)

print(
    "COMPLETE LONGITUDINAL COHORT"
)

print(
    "=" * 72
)

print(
    "\nParticipants with exactly one "
    "before + 7d + 14d + 30d run: "
    f"{len(complete)}"
)

if len(complete):

    print()

    print(
        complete[
            TIMEPOINTS
        ].to_string()
    )


# ============================================================
# Select cohort runs
# ============================================================

cohort_runs = parsed[
    parsed["participant"]
    .isin(
        complete_ids
    )
].copy()

cohort_runs[
    "timepoint_order"
] = (
    cohort_runs[
        "timepoint"
    ]
    .map(
        TIMEPOINT_ORDER
    )
)

cohort_runs = cohort_runs.sort_values(
    [
        "participant",
        "timepoint_order",
        "run_accession",
    ]
)


# ============================================================
# Validate exactly four rows per participant
# ============================================================

participant_counts = (
    cohort_runs
    .groupby(
        "participant"
    )
    .size()
)

bad_counts = participant_counts[
    participant_counts != 4
]

if len(
    bad_counts
):

    raise RuntimeError(
        "Validation failure: some complete participants "
        "do not have exactly four selected runs:\n"
        f"{bad_counts}"
    )


# ============================================================
# Construct cohort sample IDs
# ============================================================

cohort_runs[
    "sample-id"
] = (
    cohort_runs[
        "participant"
    ]
    + "_"
    + cohort_runs[
        "timepoint"
    ]
)

if cohort_runs[
    "sample-id"
].duplicated().any():

    duplicated_ids = cohort_runs.loc[
        cohort_runs[
            "sample-id"
        ].duplicated(
            keep=False
        ),
        "sample-id",
    ]

    raise RuntimeError(
        "Duplicate sample IDs detected:\n"
        + "\n".join(
            duplicated_ids.astype(
                str
            )
        )
    )


# ============================================================
# Build clean download manifest
# ============================================================

manifest_columns = [
    "sample-id",
    "participant",
    "timepoint",
    "timepoint_order",
    "run_accession",
    "biosample",
    "library_layout",
    "run_total_spots",
    "run_total_bases",
    "ena_fastq_http_1",
    "ena_fastq_http_2",
]

available_manifest_columns = [
    column
    for column in manifest_columns
    if column in cohort_runs.columns
]

manifest = cohort_runs[
    available_manifest_columns
].copy()


# ============================================================
# Validate download URLs
# ============================================================

if (
    "ena_fastq_http_1"
    in manifest.columns
    and
    "ena_fastq_http_2"
    in manifest.columns
):

    missing_urls = manifest[
        manifest[
            "ena_fastq_http_1"
        ].isna()
        |
        manifest[
            "ena_fastq_http_2"
        ].isna()
    ]

    print(
        "\nCohort rows missing paired ENA URLs: "
        f"{len(missing_urls)}"
    )

else:

    missing_urls = pd.DataFrame()


# ============================================================
# Save outputs
# ============================================================

parsed.to_csv(
    PARSED_OUTPUT,
    sep="\t",
    index=False,
)

inventory.to_csv(
    INVENTORY_OUTPUT,
    sep="\t",
)

complete.to_csv(
    COMPLETE_OUTPUT,
    sep="\t",
)

cohort_runs.to_csv(
    COHORT_RUNS_OUTPUT,
    sep="\t",
    index=False,
)

unparsed.to_csv(
    UNPARSED_OUTPUT,
    sep="\t",
    index=False,
)

duplicates.to_csv(
    DUPLICATE_OUTPUT,
    sep="\t",
    index=False,
)

manifest.to_csv(
    MANIFEST_OUTPUT,
    sep="\t",
    index=False,
)


# ============================================================
# Summary
# ============================================================

print(
    "\n"
    + "=" * 72
)

print(
    "COHORT SUMMARY"
)

print(
    "=" * 72
)

print(
    f"\nComplete participants: {len(complete)}"
)

print(
    f"Cohort runs: {len(cohort_runs)}"
)

print(
    "Expected cohort runs: "
    f"{len(complete) * 4}"
)

print(
    "Runs retained for manual review / donor analysis: "
    f"{len(unparsed)}"
)

print(
    "Donor/non-recipient records: "
    f"{len(donor_runs)}"
)

print(
    "Duplicate participant/timepoint combinations: "
    f"{len(duplicates)}"
)

if not missing_urls.empty:

    print(
        "WARNING: cohort runs with missing ENA paired FASTQ URLs: "
        f"{len(missing_urls)}"
    )

print(
    "\nSelected cohort examples:"
)

display_columns = [
    column
    for column in [
        "sample-id",
        "run_accession",
        "biosample",
    ]
    if column in cohort_runs.columns
]

print(
    cohort_runs[
        display_columns
    ]
    .head(24)
    .to_string(
        index=False
    )
)


# ============================================================
# Output files
# ============================================================

print(
    "\n"
    + "=" * 72
)

print(
    "OUTPUTS"
)

print(
    "=" * 72
)

for path in [
    PARSED_OUTPUT,
    INVENTORY_OUTPUT,
    COMPLETE_OUTPUT,
    COHORT_RUNS_OUTPUT,
    UNPARSED_OUTPUT,
    DUPLICATE_OUTPUT,
    MANIFEST_OUTPUT,
]:

    print(
        path
    )


# ============================================================
# Conclusion
# ============================================================

print(
    "\n"
    + "=" * 72
)

print(
    "PHASE 8 CONCLUSION"
)

print(
    "=" * 72
)

print(
    """
PRJNA298590 SRA metadata were parsed to identify recipient
IDs and standardized FMT timepoints.

Only paired-end sequencing runs were considered.

Compact sample names were parsed by identifying and removing
the timepoint suffix before extracting the participant ID.

Participants were retained only when exactly one sequencing
run was available for each required timepoint:

Before FMT
7 days after FMT
14 days after FMT
30 days after FMT

Donor records were kept outside the recipient longitudinal
cohort and retained for possible future donor-recipient
engraftment analyses.

Duplicate and otherwise unassigned records were preserved
for manual review rather than silently discarded.

A validated cohort download manifest was generated for the
next stage of multi-participant processing.
"""
)

print(
    "=" * 72
)

print(
    "END PHASE 8"
)

print(
    "=" * 72
)
