# Longitudinal 16S Microbiome Analysis Following Fecal Microbiota Transplantation

## Overview

This repository contains a reproducible computational analysis of longitudinal gut microbiome data following fecal microbiota transplantation (FMT).

The project uses publicly available 16S rRNA gene amplicon sequencing data from NCBI BioProject **PRJNA298590**.

The current pilot analysis focuses on participant **R009**, sampled before FMT and at three post-FMT timepoints.

The goal of this project is to build a transparent and reproducible workflow for examining longitudinal changes in gut microbial community composition following FMT.

---

## Research Question

**How does gut microbial community composition change over time following fecal microbiota transplantation?**

The pilot analysis establishes the sequencing-processing workflow needed to subsequently examine taxonomic composition and longitudinal community changes.

---

## Pilot Dataset

The pilot dataset consists of four paired-end sequencing samples from participant R009:

| Sample | Timepoint | SRA Run |
|---|---|---|
| R009_before | Before FMT | SRR2657981 |
| R009_7d | 7 days | SRR2657995 |
| R009_14d | 14 days | SRR2657993 |
| R009_30d | 30 days | SRR2657994 |

Raw sequencing data are retrieved from NCBI SRA and are intentionally not stored in this repository.

---

## Analysis Workflow

The project is organized into reproducible analysis phases.

### Phase 1 — Data Acquisition and Quality Control

The initial workflow includes:

- identification of the appropriate NCBI SRA accessions
- SRA data acquisition
- paired-end FASTQ generation
- sample metadata organization
- FASTQ integrity verification
- read-count verification
- FastQC analysis
- inspection of forward and reverse read-quality profiles
- QIIME 2 paired-end read-quality visualization

Raw sequencing files were obtained using the NCBI SRA Toolkit and converted to paired-end FASTQ files using `fasterq-dump`.

FASTQ files were subsequently compressed using gzip for QIIME 2 import.

### Phase 2 — Read Preprocessing and ASV Inference

Paired-end reads were processed using QIIME 2 and DADA2.

Processing included:

- paired-end FASTQ import
- quality assessment
- quality filtering
- DADA2 denoising
- paired-end read merging
- chimera detection and removal
- amplicon sequence variant (ASV) inference
- generation of representative sequences
- read-retention assessment

### DADA2 Parameters

Based on FastQC results and QIIME 2 quality profiles, the following truncation lengths were selected:

| Parameter | Value |
|---|---:|
| Forward truncation | 230 bp |
| Reverse truncation | 190 bp |
| Minimum overlap | 12 bp |
| Maximum expected errors, forward | 2 |
| Maximum expected errors, reverse | 2 |

---

## Primer Handling

Expected primer motifs were investigated directly in the sequencing reads.

Primer trimming was skipped because the expected primer motifs were detected in only a very small fraction of reads, consistent with primer sequences having already been removed during upstream processing.

The decision to skip primer trimming is documented explicitly in the Phase 2 workflow rather than silently assuming primer removal.

---

## Phase 2 Results

DADA2 successfully generated an ASV feature table and representative sequences for all four longitudinal samples.

A total of **166 ASVs** were observed across the pilot dataset.

### DADA2 Read Retention

| Sample | Input Reads | Filtered | Merged | Final Non-chimeric Reads | Final Retention |
|---|---:|---:|---:|---:|---:|
| R009_before | 37,985 | 31,588 | 27,445 | 21,636 | 56.96% |
| R009_7d | 40,010 | 21,694 | 19,200 | 15,599 | 38.99% |
| R009_14d | 30,280 | 23,895 | 21,596 | 18,849 | 62.25% |
| R009_30d | 42,129 | 24,182 | 21,257 | 19,484 | 46.25% |

All four samples retained more than **15,000 non-chimeric reads** following DADA2 processing.

The 230-bp forward and 190-bp reverse truncation settings were therefore retained for downstream pilot analysis.

### Observed ASVs by Sample

| Sample | Final Reads | Observed ASVs |
|---|---:|---:|
| R009_before | 21,636 | 93 |
| R009_7d | 15,599 | 62 |
| R009_14d | 18,849 | 59 |
| R009_30d | 19,484 | 68 |

Across all four samples, **166 unique ASVs** were represented in the feature table.

These values describe sequencing output following filtering, denoising, merging, and chimera removal. Taxonomic identities have not yet been assigned.

---

## Reproducibility

The repository contains scripts documenting each major stage of the workflow:

    scripts/
    ├── 00_environment_and_data_setup.sh
    ├── 01_raw_data_and_qc.sh
    └── 02_preprocess_and_denoise.sh

### 00_environment_and_data_setup.sh

Documents:

- WSL/Ubuntu project setup
- project directory structure
- SRA Toolkit setup
- sample accession mapping
- SRA acquisition
- FASTQ generation
- FASTQ compression
- Conda installation
- QIIME 2 environment configuration
- QIIME 2 import
- quality-summary generation
- truncation decisions
- primer assessment

This script primarily serves as a reproducibility and setup record.

### 01_raw_data_and_qc.sh

Documents the raw sequencing-data verification and quality-control workflow.

### 02_preprocess_and_denoise.sh

Documents the QIIME 2 preprocessing and DADA2 ASV-inference workflow, including the final processing parameters.

---

## Software Environment

The analysis was performed under **WSL2 / Ubuntu**.

Primary software includes:

- SRA Toolkit 3.2.1
- FastQC 0.12.1
- Miniconda
- Python 3.12.13
- QIIME 2 / Rachis 2026.7
- DADA2

The QIIME 2 analysis environment used for the current workflow was:

    rachis-qiime2-2026.7

The environment includes QIIME 2 plugins for DADA2, feature-table analysis, taxonomic classification, diversity analysis, and related microbiome workflows.

---

## Repository Structure

    microbiome-fmt/
    ├── README.md
    ├── .gitignore
    │
    ├── data/
    │   └── metadata/
    │       ├── samples.tsv
    │       └── fastq_manifest.tsv
    │
    ├── scripts/
    │   ├── 00_environment_and_data_setup.sh
    │   ├── 01_raw_data_and_qc.sh
    │   └── 02_preprocess_and_denoise.sh
    │
    ├── logs/
    │   ├── 01_raw_data_and_qc.txt
    │   └── 02_preprocess_and_denoise.txt
    │
    └── results/
        └── phase2/
            ├── denoising-stats-export/
            │   └── stats.tsv
            ├── feature-frequencies-export/
            │   └── metadata.tsv
            └── sample-frequencies-export/
                └── metadata.tsv

---

## Data Availability

Raw sequencing files are intentionally excluded from this repository.

The raw data are publicly available through the NCBI Sequence Read Archive and can be retrieved using the accession information documented in this repository.

**NCBI BioProject:** PRJNA298590

Pilot SRA runs:

- SRR2657981
- SRR2657995
- SRR2657993
- SRR2657994

Large raw FASTQ files, SRA archives, generated FastQC reports, and QIIME 2 artifacts are excluded from version control through `.gitignore`.

This keeps the repository focused on reproducible code, metadata, processing decisions, logs, and compact analysis results.

---

## Current Project Status

### Completed

**Phase 1 — Data acquisition and quality control**

- SRA acquisition
- FASTQ generation
- metadata creation
- raw-read verification
- FastQC
- QIIME 2 quality assessment

**Phase 2 — Read preprocessing and ASV inference**

- paired-end QIIME 2 import
- DADA2 filtering
- denoising
- paired-end merging
- chimera removal
- ASV generation
- representative-sequence generation
- read-retention assessment

### Next

**Phase 3 — Taxonomic Classification**

The next stage will assign bacterial taxonomy to the inferred ASVs.

Subsequent analysis will examine longitudinal changes in microbial community composition across:

**Before FMT → 7 days → 14 days → 30 days**

Future stages may include:

- taxonomic classification
- taxonomic composition visualization
- alpha diversity
- beta diversity
- longitudinal community analysis
- visualization of microbiome recovery following FMT

---

## Project Notes

This repository is being developed as a reproducible bioinformatics workflow.

Processing decisions, including quality-based truncation and primer handling, are documented explicitly so that the analysis can be reviewed, reproduced, and modified as the project develops.
