#!/usr/bin/env bash

# Example configuration for the Illumina–ONT small-variant concordance workflow.
# Copy this file to config/config.sh and replace the placeholder paths.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

REFERENCE_FASTA="/path/to/GRCh38_no_alt_analysis_set.fna"
PANEL_BED="$PROJECT_DIR/resources/autoinflammatory.bed"

MIN_ONT_DEPTH=20
MAIN_CONTIG_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

METADATA_DIR="$PROJECT_DIR/metadata"
SCRIPT_DIR="$PROJECT_DIR/scripts"
WORK_DIR="$PROJECT_DIR/work"
RESULTS_DIR="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

SAMPLE_PAIRS="$METADATA_DIR/sample_pairs.tsv"
