#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR_REAL="$(
    cd "$(dirname "${BASH_SOURCE[0]}")" &&
    pwd
)"

PROJECT_DIR_REAL="$(
    cd "$SCRIPT_DIR_REAL/.." &&
    pwd
)"

CONFIG_FILE="$PROJECT_DIR_REAL/config/config.sh"

# shellcheck source=/dev/null
source "$CONFIG_FILE"

mkdir -p \
    "$WORK_DIR" \
    "$RESULTS_DIR" \
    "$LOG_DIR"

echo "Step 1/4: validating inputs"

python3 "$SCRIPT_DIR_REAL/validate_inputs.py" \
    --reference "$REFERENCE_FASTA" \
    --bed "$PANEL_BED" \
    --sample-pairs "$SAMPLE_PAIRS" \
    --output "$RESULTS_DIR/input_validation.tsv"

echo
echo "Step 2/4: preparing normalized panel-restricted VCFs"

"$SCRIPT_DIR_REAL/prepare_vcfs.sh" \
    "$CONFIG_FILE"

echo
echo "Step 3/4: comparing small-variant records"

tail -n +2 "$SAMPLE_PAIRS" |
while IFS=$'\t' read -r \
    sample_key \
    comparison_id \
    illumina_sample_id \
    ont_sample_id \
    illumina_vcf \
    ont_vcf \
    ont_alignment
do
    [[ -z "$sample_key" ]] && continue

    pair_work_dir="$WORK_DIR/$comparison_id"
    pair_results_dir="$RESULTS_DIR/$comparison_id"

    mkdir -p "$pair_results_dir"

    python3 "$SCRIPT_DIR_REAL/compare_small_variants.py" \
        --sample-key "$sample_key" \
        --comparison-id "$comparison_id" \
        --illumina-vcf \
            "$pair_work_dir/illumina.prepared.vcf.gz" \
        --ont-vcf \
            "$pair_work_dir/ont.prepared.vcf.gz" \
        --ont-alignment "$ont_alignment" \
        --reference "$REFERENCE_FASTA" \
        --minimum-ont-depth "$MIN_ONT_DEPTH" \
        --output-dir "$pair_results_dir"
done

echo
echo "Step 4/4: combining comparison summaries"

python3 "$SCRIPT_DIR_REAL/combine_summaries.py" \
    --results-dir "$RESULTS_DIR" \
    --output \
        "$RESULTS_DIR/concordance_summary_all_samples.tsv"

echo
echo "Analysis completed."
echo "Combined summary:"
echo "$RESULTS_DIR/concordance_summary_all_samples.tsv"
