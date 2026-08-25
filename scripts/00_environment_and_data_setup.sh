#!/usr/bin/env bash

# ============================================================
# PROJECT ENVIRONMENT AND DATA SETUP HISTORY
# ============================================================
#
# Project:
#   microbiome-fmt
#
# Purpose:
#   Document the commands and environment setup used to prepare
#   the project for reproducible microbiome analysis.
#
# This script is primarily a reproducibility/reference record.
# Some installation and download commands are commented out to
# prevent accidental reinstallation or repeated downloads.
#
# ============================================================


# ------------------------------------------------------------
# 1. WINDOWS / WSL ENVIRONMENT
# ------------------------------------------------------------
#
# WSL2 was enabled on Windows.
#
# Ubuntu was installed with:
#
#   wsl --install -d Ubuntu
#
# Linux username:
#
#   orippee
#
# Project Windows path:
#
#   C:\Users\ORippee\OneDrive - Kemin Industries\microbiome-fmt
#
# Equivalent WSL path:
#
#   /mnt/c/Users/ORippee/OneDrive - Kemin Industries/microbiome-fmt
#
# Enter project:
#
# cd "/mnt/c/Users/ORippee/OneDrive - Kemin Industries/microbiome-fmt"


# ------------------------------------------------------------
# 2. PROJECT DIRECTORY STRUCTURE
# ------------------------------------------------------------

mkdir -p data/raw
mkdir -p data/sra
mkdir -p data/metadata

mkdir -p scripts
mkdir -p logs

mkdir -p results/qc/raw
mkdir -p results/qiime2
mkdir -p results/phase2


# ------------------------------------------------------------
# 3. SRA TOOLKIT
# ------------------------------------------------------------
#
# Ubuntu package indexes were updated:
#
# sudo apt update
#
# SRA Toolkit was installed:
#
# sudo apt install sra-toolkit
#
# Version verified:
#
# prefetch --version
#
# Observed version:
#
#   prefetch : 3.2.1


# ------------------------------------------------------------
# 4. SAMPLE IDENTIFICATION
# ------------------------------------------------------------
#
# BioProject:
#
#   PRJNA298590
#
# Study:
#
#   Fecal microbiota transplantation for recurrent
#   C. difficile infection
#
# Sequencing:
#
#   16S V4 amplicon
#   Illumina MiSeq
#   paired-end
#   250 bp reads
#
# Selected pilot participant:
#
#   R009
#
# Samples:
#
#   R009_before -> SAMN04161322 -> SRR2657981
#   R009_7d     -> SAMN04161321 -> SRR2657995
#   R009_14d    -> SAMN04161319 -> SRR2657993
#   R009_30d    -> SAMN04161320 -> SRR2657994


# ------------------------------------------------------------
# 5. SRA DOWNLOAD
# ------------------------------------------------------------
#
# Downloads were performed using prefetch.
#
# Commands used:
#
# prefetch SRR2657981
# prefetch SAMN04161321
# prefetch SAMN04161319
# prefetch SAMN04161320
#
# BioSample accessions were resolved by SRA Toolkit to:
#
# SAMN04161321 -> SRR2657995
# SAMN04161319 -> SRR2657993
# SAMN04161320 -> SRR2657994


# ------------------------------------------------------------
# 6. SRA -> FASTQ CONVERSION
# ------------------------------------------------------------
#
# Paired FASTQ files were generated using fasterq-dump.
#
# Commands:
#
# fasterq-dump SRR2657981 --split-files --outdir data/raw
# fasterq-dump SRR2657995 --split-files --outdir data/raw
# fasterq-dump SRR2657993 --split-files --outdir data/raw
# fasterq-dump SRR2657994 --split-files --outdir data/raw
#
# Output FASTQ pairs:
#
# data/raw/SRR2657981_1.fastq
# data/raw/SRR2657981_2.fastq
#
# data/raw/SRR2657995_1.fastq
# data/raw/SRR2657995_2.fastq
#
# data/raw/SRR2657993_1.fastq
# data/raw/SRR2657993_2.fastq
#
# data/raw/SRR2657994_1.fastq
# data/raw/SRR2657994_2.fastq


# ------------------------------------------------------------
# 7. SRA ARCHIVE ORGANIZATION
# ------------------------------------------------------------
#
# SRA archives were moved under:
#
# data/sra/
#
# Example:
#
# data/sra/SRR2657981/SRR2657981.sra
#
# Verification:
#
# find data/sra -name "*.sra" -type f


# ------------------------------------------------------------
# 8. SAMPLE METADATA
# ------------------------------------------------------------
#
# Created:
#
# data/metadata/samples.tsv
#
# Current contents:
#
# sample-id       timepoint       srr
# R009_before     before          SRR2657981
# R009_7d         7d              SRR2657995
# R009_14d        14d             SRR2657993
# R009_30d        30d             SRR2657994
#
# Tab-delimited format was verified using:
#
# cat -A data/metadata/samples.tsv
#
# ^I characters confirmed tab separators.


# ------------------------------------------------------------
# 9. RAW FASTQ VERIFICATION
# ------------------------------------------------------------
#
# FASTQ read counts were verified.
#
# Observed counts:
#
# SRR2657981:
#   37,985 reads per direction
#
# SRR2657993:
#   30,280 reads per direction
#
# SRR2657994:
#   42,129 reads per direction
#
# SRR2657995:
#   40,010 reads per direction
#
# Representative FASTQ content was inspected using:
#
# head -n 8 data/raw/SRR2657981_1.fastq
#
# and:
#
# head -n 2 data/raw/SRR2657981_2.fastq


# ------------------------------------------------------------
# 10. FASTQC
# ------------------------------------------------------------
#
# FastQC version:
#
#   FastQC v0.12.1
#
# Version check:
#
# fastqc --version
#
# Raw-read reports were generated under:
#
# results/qc/raw/
#
# FastQC reports were inspected for:
#
# - per-base quality
# - per-sequence quality
# - GC distribution
# - adapter contamination
# - overrepresented sequences
#
# Key QC observations are documented in:
#
# logs/01_raw_data_and_qc.txt


# ------------------------------------------------------------
# 11. CONDA INSTALLATION
# ------------------------------------------------------------
#
# The system QIIME 2 package used Python 3.14 and was
# incompatible with QIIME 2 2024.5.
#
# Miniconda was therefore installed to create an isolated,
# compatible QIIME 2 environment.
#
# Installer:
#
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#
# Installation:
#
# bash Miniconda3-latest-Linux-x86_64.sh
#
# Shell activation:
#
# eval "$(/home/orippee/miniconda3/bin/conda shell.bash hook)"
#
# conda init bash
#
# source ~/.bashrc
#
# Conda version observed:
#
#   conda 26.5.3


# ------------------------------------------------------------
# 12. QIIME 2 ENVIRONMENT
# ------------------------------------------------------------
#
# QIIME 2 / Rachis environment:
#
#   rachis-qiime2-2026.7
#
# Activation:
#
# conda activate rachis-qiime2-2026.7
#
# Observed environment:
#
# Python: 3.12.13
# rachis: 2026.7.0
# q2cli: 2026.7.0
#
# Installed plugins include:
#
# - dada2
# - cutadapt
# - demux
# - diversity
# - feature-classifier
# - feature-table
# - taxa
# - longitudinal
#
# Environment check:
#
# qiime info


# ------------------------------------------------------------
# 13. FASTQ COMPRESSION
# ------------------------------------------------------------
#
# QIIME 2's paired-end manifest import expected gzipped FASTQ.
#
# FASTQ files were compressed with:
#
# gzip data/raw/*.fastq
#
# This is lossless compression.
#
# Output files now end in:
#
# .fastq.gz


# ------------------------------------------------------------
# 14. QIIME 2 FASTQ MANIFEST
# ------------------------------------------------------------
#
# Created:
#
# data/metadata/fastq_manifest.tsv
#
# Columns:
#
# sample-id
# forward-absolute-filepath
# reverse-absolute-filepath
#
# Manifest tab formatting was verified with:
#
# cat -A data/metadata/fastq_manifest.tsv


# ------------------------------------------------------------
# 15. QIIME 2 IMPORT
# ------------------------------------------------------------
#
# Imported the paired-end FASTQs using:
#
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path data/metadata/fastq_manifest.tsv \
#   --output-path results/qiime2/demux-paired-end.qza \
#   --input-format PairedEndFastqManifestPhred33V2
#
# Artifact verification:
#
# qiime tools peek results/qiime2/demux-paired-end.qza
#
# Observed type:
#
# SampleData[PairedEndSequencesWithQuality]


# ------------------------------------------------------------
# 16. QIIME 2 QUALITY SUMMARY
# ------------------------------------------------------------
#
# Generated:
#
# qiime demux summarize \
#   --i-data results/qiime2/demux-paired-end.qza \
#   --o-visualization results/qiime2/demux-paired-end.qzv
#
# Visualization:
#
# results/qiime2/demux-paired-end.qzv
#
# This was inspected in QIIME 2 View.


# ------------------------------------------------------------
# 17. TRUNCATION DECISION
# ------------------------------------------------------------
#
# Based on FastQC and QIIME 2 quality plots:
#
# Forward reads:
#   generally high quality across most of the read
#
# Reverse reads:
#   increasing variability after ~160 bp
#   substantial quality degradation after ~200 bp
#   very poor quality near ~235-250 bp
#
# Initial DADA2 truncation parameters selected:
#
#   forward: 230 bp
#   reverse: 190 bp


# ------------------------------------------------------------
# 18. PRIMER ASSESSMENT
# ------------------------------------------------------------
#
# Standard V4 primer motifs were checked against the reads.
#
# Commands:
#
# zgrep -c "GTGCCAG" data/raw/SRR2657981_1.fastq.gz
#
# Observed:
#
#   19 matches
#
# zgrep -c "GGACTAC" data/raw/SRR2657981_2.fastq.gz
#
# Observed:
#
#   59 matches
#
# Relative to tens of thousands of reads, these counts are
# very small.
#
# Interpretation:
#
# Primer trimming will be skipped in Phase 2 because the
# motifs are present in only a very small fraction of reads,
# consistent with primer removal having occurred upstream.


# ------------------------------------------------------------
# 19. CURRENT ANALYSIS STATE
# ------------------------------------------------------------
#
# Phase 1:
#
# COMPLETE
#
# Phase 2:
#
# Script created:
#
# scripts/02_preprocess_and_denoise.sh
#
# Current parameters:
#
# TRUNC_LEN_F=230
# TRUNC_LEN_R=190
#
# Primer trimming:
#
# skipped
#
# Next major action:
#
# Run Phase 2 DADA2 preprocessing and inspect denoising /
# read-retention statistics.


# ============================================================
# END ENVIRONMENT AND DATA SETUP HISTORY
# ============================================================
