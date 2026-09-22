#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_key",
    "comparison_id",
    "illumina_sample_id",
    "ont_sample_id",
    "illumina_vcf",
    "ont_vcf",
    "ont_alignment",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate inputs for the Illumina–ONT small-variant "
            "concordance analysis."
        )
    )
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--sample-pairs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )


def single_vcf_sample(vcf_path: Path) -> tuple[str | None, str]:
    result = run_command(
        ["bcftools", "query", "-l", str(vcf_path)]
    )

    if result.returncode != 0:
        return None, result.stderr.strip()

    samples = [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip()
    ]

    if len(samples) != 1:
        return None, (
            f"Expected exactly one sample, found {len(samples)}: "
            f"{','.join(samples)}"
        )

    return samples[0], ""


def check_alignment(alignment_path: Path) -> tuple[bool, str]:
    result = run_command(
        ["samtools", "quickcheck", "-v", str(alignment_path)]
    )

    if result.returncode == 0:
        return True, ""

    message = result.stdout.strip() or result.stderr.strip()
    return False, message


def check_required_tools() -> list[str]:
    required_tools = [
        "bcftools",
        "samtools",
        "bgzip",
        "tabix",
        "python3",
    ]

    return [
        tool
        for tool in required_tools
        if shutil.which(tool) is None
    ]


def main() -> int:
    args = parse_arguments()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    audit_rows: list[dict[str, str]] = []

    missing_tools = check_required_tools()

    for tool in missing_tools:
        errors.append(f"Required tool not found: {tool}")

    common_inputs = [
        ("reference_fasta", args.reference),
        ("reference_fai", Path(f"{args.reference}.fai")),
        ("panel_bed", args.bed),
        ("sample_pairs", args.sample_pairs),
    ]

    for label, path in common_inputs:
        if not path.is_file():
            errors.append(f"Missing {label}: {path}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    with args.sample_pairs.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            print(
                "ERROR: sample_pairs.tsv has no header.",
                file=sys.stderr,
            )
            return 1

        missing_columns = REQUIRED_COLUMNS - set(reader.fieldnames)

        if missing_columns:
            print(
                "ERROR: Missing metadata columns: "
                + ", ".join(sorted(missing_columns)),
                file=sys.stderr,
            )
            return 1

        rows = list(reader)

    if not rows:
        print(
            "ERROR: sample_pairs.tsv contains no sample rows.",
            file=sys.stderr,
        )
        return 1

    seen_keys: set[str] = set()
    seen_comparisons: set[str] = set()

    for row_number, row in enumerate(rows, start=2):
        sample_key = row["sample_key"].strip()
        comparison_id = row["comparison_id"].strip()
        illumina_sample_id = row["illumina_sample_id"].strip()
        ont_sample_id = row["ont_sample_id"].strip()

        illumina_vcf = Path(row["illumina_vcf"].strip())
        ont_vcf = Path(row["ont_vcf"].strip())
        ont_alignment = Path(row["ont_alignment"].strip())

        row_errors: list[str] = []

        if not sample_key:
            row_errors.append("sample_key is empty")
        elif not re.fullmatch(r"\d+", sample_key):
            row_errors.append(
                f"sample_key is not numeric: {sample_key}"
            )

        if sample_key in seen_keys:
            row_errors.append(
                f"duplicate sample_key: {sample_key}"
            )
        seen_keys.add(sample_key)

        if comparison_id in seen_comparisons:
            row_errors.append(
                f"duplicate comparison_id: {comparison_id}"
            )
        seen_comparisons.add(comparison_id)

        expected_comparison_id = (
            f"{illumina_sample_id}_vs_{ont_sample_id.replace('_', '')}"
        )

        if comparison_id != expected_comparison_id:
            row_errors.append(
                "comparison_id does not follow the expected form: "
                f"expected {expected_comparison_id}, "
                f"observed {comparison_id}"
            )

        identifiers_to_check = {
            "illumina_sample_id": illumina_sample_id,
            "ont_sample_id": ont_sample_id,
            "comparison_id": comparison_id,
            "illumina_vcf_path": str(illumina_vcf),
            "ont_vcf_path": str(ont_vcf),
            "ont_alignment_path": str(ont_alignment),
        }

        for label, value in identifiers_to_check.items():
            if sample_key not in value:
                row_errors.append(
                    f"{label} does not contain sample key "
                    f"{sample_key}: {value}"
                )

        for label, path in [
            ("illumina_vcf", illumina_vcf),
            ("ont_vcf", ont_vcf),
            ("ont_alignment", ont_alignment),
        ]:
            if not path.is_file():
                row_errors.append(
                    f"missing {label}: {path}"
                )

        illumina_vcf_sample = ""
        ont_vcf_sample = ""
        alignment_status = "NOT_TESTED"

        if illumina_vcf.is_file():
            sample_name, message = single_vcf_sample(illumina_vcf)

            if sample_name is None:
                row_errors.append(
                    f"could not validate Illumina VCF sample: "
                    f"{message}"
                )
            else:
                illumina_vcf_sample = sample_name

                if sample_name != illumina_sample_id:
                    row_errors.append(
                        "Illumina VCF sample mismatch: "
                        f"metadata={illumina_sample_id}, "
                        f"VCF={sample_name}"
                    )

        if ont_vcf.is_file():
            sample_name, message = single_vcf_sample(ont_vcf)

            if sample_name is None:
                row_errors.append(
                    f"could not validate ONT VCF sample: "
                    f"{message}"
                )
            else:
                ont_vcf_sample = sample_name

                if sample_name != ont_sample_id:
                    row_errors.append(
                        "ONT VCF sample mismatch: "
                        f"metadata={ont_sample_id}, "
                        f"VCF={sample_name}"
                    )

        if ont_alignment.is_file():
            alignment_ok, message = check_alignment(ont_alignment)

            if alignment_ok:
                alignment_status = "PASS"
            else:
                alignment_status = "FAIL"
                row_errors.append(
                    f"ONT alignment quickcheck failed: {message}"
                )

        audit_rows.append(
            {
                "sample_key": sample_key,
                "comparison_id": comparison_id,
                "illumina_sample_id_metadata": illumina_sample_id,
                "illumina_sample_id_vcf": illumina_vcf_sample,
                "ont_sample_id_metadata": ont_sample_id,
                "ont_sample_id_vcf": ont_vcf_sample,
                "ont_alignment_quickcheck": alignment_status,
                "validation_status": (
                    "PASS" if not row_errors else "FAIL"
                ),
                "validation_messages": "; ".join(row_errors),
                "illumina_vcf": str(illumina_vcf),
                "ont_vcf": str(ont_vcf),
                "ont_alignment": str(ont_alignment),
            }
        )

        for error in row_errors:
            errors.append(
                f"Metadata row {row_number} "
                f"({comparison_id}): {error}"
            )

    fieldnames = [
        "sample_key",
        "comparison_id",
        "illumina_sample_id_metadata",
        "illumina_sample_id_vcf",
        "ont_sample_id_metadata",
        "ont_sample_id_vcf",
        "ont_alignment_quickcheck",
        "validation_status",
        "validation_messages",
        "illumina_vcf",
        "ont_vcf",
        "ont_alignment",
    ]

    with args.output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(audit_rows)

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)

        print(
            f"Input validation failed. Audit: {args.output}",
            file=sys.stderr,
        )
        return 1

    print(f"Input validation passed: {len(audit_rows)} pairs")
    print(f"Audit written to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
