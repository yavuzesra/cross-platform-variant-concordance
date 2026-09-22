#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/config.sh" >&2
    exit 1
fi

CONFIG_FILE="$1"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_FILE"

mkdir -p \
    "$WORK_DIR" \
    "$RESULTS_DIR" \
    "$LOG_DIR"

MAIN_CONTIGS="$WORK_DIR/main_contigs.txt"

cat > "$MAIN_CONTIGS" <<'EOF'
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
EOF

MAIN_CONTIG_LIST="$(
    paste -sd, "$MAIN_CONTIGS"
)"

prepare_one_vcf() {
    local input_vcf="$1"
    local output_vcf="$2"
    local log_file="$3"
    local label="$4"

    mkdir -p "$(dirname "$output_vcf")"
    mkdir -p "$(dirname "$log_file")"

    {
        echo "label=$label"
        echo "input_vcf=$input_vcf"
        echo "output_vcf=$output_vcf"
        echo "reference_fasta=$REFERENCE_FASTA"
        echo "panel_bed=$PANEL_BED"
        echo "filter_policy=FILTER equals PASS"
        echo "normalization=bcftools norm -m -any"
        echo "reference_mismatch_policy=exclude"
        echo "main_contigs=$MAIN_CONTIG_LIST"
        echo
        echo "bcftools_version:"
        bcftools --version | head -n 1
        echo
        echo "normalization_messages:"
    } > "$log_file"

    bcftools view \
        -f PASS \
        -t "$MAIN_CONTIG_LIST" \
        -T "$PANEL_BED" \
        -Ou \
        "$input_vcf" \
    | bcftools norm \
        -f "$REFERENCE_FASTA" \
        -m -any \
        -c x \
        -Oz \
        -o "$output_vcf" \
        2>> "$log_file"

    bcftools index \
        --force \
        --tbi \
        "$output_vcf"

    {
        echo
        echo "prepared_record_count:"
        bcftools view \
            -H \
            "$output_vcf" \
        | wc -l
    } >> "$log_file"
}

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
    if [[ -z "$sample_key" ]]; then
        continue
    fi

    pair_work_dir="$WORK_DIR/$comparison_id"
    pair_log_dir="$LOG_DIR/$comparison_id"

    mkdir -p \
        "$pair_work_dir" \
        "$pair_log_dir"

    illumina_output="$pair_work_dir/illumina.prepared.vcf.gz"
    ont_output="$pair_work_dir/ont.prepared.vcf.gz"

    echo "Preparing $comparison_id"

    prepare_one_vcf \
        "$illumina_vcf" \
        "$illumina_output" \
        "$pair_log_dir/illumina_preparation.log" \
        "${comparison_id}:Illumina"

    prepare_one_vcf \
        "$ont_vcf" \
        "$ont_output" \
        "$pair_log_dir/ont_preparation.log" \
        "${comparison_id}:ONT"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$sample_key" \
        "$comparison_id" \
        "$illumina_output" \
        "$ont_output" \
        "$ont_alignment" \
        > "$pair_work_dir/prepared_inputs.tsv"
done

echo "VCF preparation completed."
