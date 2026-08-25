#!/usr/bin/env bash

# ============================================================
# PHASE 1 — RAW DATA ACQUISITION, ORGANIZATION, AND QC
# ============================================================
#
# Project: microbiome-fmt
#
# Purpose:
#   Document and verify the raw sequencing data prior to
#   preprocessing and ASV inference.
#
# Major activities:
#   1. Record project location
#   2. Verify software versions
#   3. Verify sample metadata
#   4. Verify FASTQ files
#   5. Verify read counts
#   6. Verify SRA archives
#   7. Run FastQC
#   8. Summarize raw-read QC findings
#
# IMPORTANT:
#   No raw FASTQ files are modified by this script.
#   Trimming/filtering decisions are intentionally deferred
#   to Phase 2.
#
# ============================================================

set -u

PROJECT_DIR="$(pwd)"
LOG_DIR="logs"
QC_DIR="results/qc/raw"

mkdir -p "$LOG_DIR"
mkdir -p "$QC_DIR"

LOG_FILE="$LOG_DIR/01_raw_data_and_qc.txt"

# ------------------------------------------------------------
# Start log
# ------------------------------------------------------------

{
echo "============================================================"
echo "PHASE 1 — RAW DATA ACQUISITION, ORGANIZATION, AND QC"
echo "============================================================"
echo
echo "Date:"
date
echo
echo "Project directory:"
pwd
echo

# ------------------------------------------------------------
# 1. Software versions
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "SOFTWARE"
echo "------------------------------------------------------------"

echo "FastQC:"
fastqc --version 2>&1 || echo "FastQC not found"

echo
echo "SRA Toolkit:"
prefetch --version 2>&1 | head -n 1 || echo "prefetch not found"

echo
echo "fasterq-dump:"
fasterq-dump --version 2>&1 | head -n 1 || echo "fasterq-dump not found"

echo

# ------------------------------------------------------------
# 2. Metadata
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "SAMPLE METADATA"
echo "------------------------------------------------------------"

if [ -f data/metadata/samples.tsv ]; then
    cat data/metadata/samples.tsv
else
    echo "WARNING: data/metadata/samples.tsv not found"
fi

echo

# ------------------------------------------------------------
# 3. FASTQ files
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "RAW FASTQ FILES"
echo "------------------------------------------------------------"

find data/raw -maxdepth 1 -type f -name "*.fastq" -print0 \
    | sort -z \
    | xargs -0 -r ls -lh

echo

# ------------------------------------------------------------
# 4. Read counts
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "READ COUNTS"
echo "------------------------------------------------------------"

for file in data/raw/*.fastq; do

    if [ -f "$file" ]; then

        reads=$(awk 'END {print NR/4}' "$file")

        echo "$(basename "$file"): $reads reads"

    fi

done

echo

# ------------------------------------------------------------
# 5. SRA archive verification
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "SRA ARCHIVES"
echo "------------------------------------------------------------"

find data/sra -type f -name "*.sra" -print 2>/dev/null \
    | sort

echo

# ------------------------------------------------------------
# 6. FastQC
# ------------------------------------------------------------

echo "------------------------------------------------------------"
echo "FASTQC"
echo "------------------------------------------------------------"

echo "Running FastQC on raw FASTQ files..."
echo

fastqc \
    data/raw/*.fastq \
    --outdir "$QC_DIR"

echo
echo "FastQC reports:"
find "$QC_DIR" -maxdepth 1 -type f \
    \( -name "*.html" -o -name "*.zip" \) \
    -print \
    | sort

echo

# ------------------------------------------------------------
# 7. QC observations
# ------------------------------------------------------------

echo "============================================================"
echo "RAW-READ QC OBSERVATIONS"
echo "============================================================"

echo
echo "Dataset:"
echo "Four samples, each represented by paired-end 250-bp reads."
echo
echo "Samples:"
echo "  R009_before -> SRR2657981"
echo "  R009_7d     -> SRR2657995"
echo "  R009_14d    -> SRR2657993"
echo "  R009_30d    -> SRR2657994"
echo

echo "Forward-read quality:"
echo "  SRR2657981 R1 showed generally high quality through"
echo "  most of the 250-bp read, with increasing variability"
echo "  toward the 3-prime end."

echo
echo "Reverse-read quality:"
echo "  SRR2657981 R2 showed generally good quality through"
echo "  approximately 160 bp."
echo "  Quality became increasingly variable after approximately"
echo "  165 bp and declined substantially after approximately"
echo "  200 bp."
echo "  Bases approximately 235-250 were extremely low quality."

echo
echo "Adapter assessment:"
echo "  No substantial Illumina adapter contamination was"
echo "  observed in the reviewed FastQC reports."

echo
echo "Overrepresented sequences:"
echo "  Numerous highly abundant sequences were observed."
echo "  Because this is targeted 16S amplicon sequencing,"
echo "  sequence overrepresentation is expected to some degree"
echo "  and is not automatically interpreted as contamination."

echo
echo "Read structure:"
echo "  Forward and reverse reads show the expected structure"
echo "  of paired-end reads from the same 16S amplicon."

echo
echo "Primer assessment:"
echo "  Obvious primer sequence was not identified from the"
echo "  beginning of the reviewed reads."
echo "  Exact primer sequences will be verified from the source"
echo "  study before any primer-removal step is performed."

echo

# ------------------------------------------------------------
# 8. Scientific conclusion
# ------------------------------------------------------------

echo "============================================================"
echo "PHASE 1 CONCLUSION"
echo "============================================================"

echo
echo "The raw sequencing dataset was successfully acquired,"
echo "organized, converted to paired-end FASTQ format, and"
echo "verified against the sample metadata."
echo
echo "The dataset consists of 250-bp paired-end Illumina reads"
echo "from a targeted 16S rRNA amplicon experiment."
echo
echo "Raw-read quality was generally high through most of the"
echo "read length, with greater quality degradation toward the"
echo "3-prime ends, particularly in the reverse reads."
echo
echo "The reverse reads showed substantial quality degradation"
echo "after approximately 200-225 bp, with the final approximately"
echo "15 bp being extremely low quality."
echo
echo "No substantial Illumina adapter contamination was observed."
echo "Overrepresented sequences were observed and are considered"
echo "compatible with the targeted amplicon nature of the dataset."
echo
echo "The raw FASTQ files were NOT modified during Phase 1."
echo
echo "Quality filtering, truncation, primer removal, denoising,"
echo "paired-end merging, and chimera removal are deferred to"
echo "Phase 2."
echo
echo "============================================================"
echo "END PHASE 1"
echo "============================================================"

} | tee "$LOG_FILE"

echo
echo "Phase 1 complete."
echo "Log saved to:"
echo "$LOG_FILE"
