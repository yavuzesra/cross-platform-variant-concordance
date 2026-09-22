#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: str
    filt: str
    gt: str
    dp: str
    ad: str
    af: str
    gq: str

    @property
    def key(self) -> tuple[str, int, str, str]:
        return self.chrom, self.pos, self.ref, self.alt


OUTPUT_COLUMNS = [
    "sample_key",
    "comparison_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_type",
    "comparison_category",
    "illumina_gt",
    "illumina_qual",
    "illumina_filter",
    "illumina_dp",
    "illumina_ad",
    "illumina_af",
    "illumina_gq",
    "ont_gt",
    "ont_qual",
    "ont_filter",
    "ont_dp",
    "ont_ad",
    "ont_af",
    "ont_gq",
    "ont_local_alignment_depth",
    "ont_depth_category",
    "illumina_callability",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare normalized, panel-restricted Illumina and ONT "
            "small-variant VCF files."
        )
    )
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--illumina-vcf", required=True, type=Path)
    parser.add_argument("--ont-vcf", required=True, type=Path)
    parser.add_argument("--ont-alignment", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--minimum-ont-depth", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(command)
            + "\n\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr
        )

    return result


def clean(value: str) -> str:
    value = value.strip()

    if value in {"", "."}:
        return ""

    return value


def normalize_genotype(genotype: str) -> str:
    genotype = clean(genotype)

    if not genotype:
        return ""

    alleles = genotype.replace("|", "/").split("/")

    if len(alleles) == 2 and all(
        allele not in {"", "."} for allele in alleles
    ):
        alleles = sorted(alleles)

    return "/".join(alleles)


def classify_variant(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if len(ref) < len(alt):
        return "insertion"

    if len(ref) > len(alt):
        return "deletion"

    return "complex_substitution"


def query_vcf(vcf_path: Path) -> dict[tuple[str, int, str, str], Variant]:
    query_format = (
        "%CHROM\\t"
        "%POS\\t"
        "%REF\\t"
        "%ALT\\t"
        "%QUAL\\t"
        "%FILTER"
        "[\\t%GT\\t%DP\\t%AD\\t%AF\\t%GQ]"
        "\\n"
    )

    result = run_command(
        [
            "bcftools",
            "query",
            "-f",
            query_format,
            str(vcf_path),
        ]
    )

    variants: dict[tuple[str, int, str, str], Variant] = {}

    for line_number, line in enumerate(
        result.stdout.splitlines(),
        start=1,
    ):
        if not line.strip():
            continue

        fields = line.split("\t")

        if len(fields) < 11:
            raise RuntimeError(
                f"Unexpected bcftools output at line {line_number}: "
                f"{line!r}"
            )

        chrom, pos, ref, alt, qual, filt, gt, dp, ad, af, gq = (
            fields[:11]
        )

        variant = Variant(
            chrom=chrom,
            pos=int(pos),
            ref=ref,
            alt=alt,
            qual=clean(qual),
            filt=clean(filt),
            gt=normalize_genotype(gt),
            dp=clean(dp),
            ad=clean(ad),
            af=clean(af),
            gq=clean(gq),
        )

        if variant.key in variants:
            raise RuntimeError(
                "Duplicate normalized variant key found in "
                f"{vcf_path}: {variant.key}"
            )

        variants[variant.key] = variant

    return variants


def calculate_ont_depth(
    variants: list[Variant],
    alignment: Path,
    reference: Path,
) -> dict[tuple[str, int], int]:
    if not variants:
        return {}

    positions = sorted(
        {(variant.chrom, variant.pos) for variant in variants}
    )

    with tempfile.TemporaryDirectory(
        prefix="ont_depth_"
    ) as temp_directory:
        temporary_bed = Path(temp_directory) / "positions.bed"

        with temporary_bed.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            for chrom, pos in positions:
                handle.write(f"{chrom}\t{pos - 1}\t{pos}\n")

        command = [
            "samtools",
            "depth",
            "-aa",
            "-G",
            "3844",
            "-q",
            "0",
            "-Q",
            "0",
            "-b",
            str(temporary_bed),
            "--reference",
            str(reference),
            str(alignment),
        ]

        result = run_command(command)

    depth_by_position = {
        position: 0
        for position in positions
    }

    for line in result.stdout.splitlines():
        if not line.strip():
            continue

        chrom, pos, depth = line.split("\t")[:3]
        depth_by_position[(chrom, int(pos))] = int(depth)

    return depth_by_position


def blank_platform_fields(prefix: str) -> dict[str, str]:
    return {
        f"{prefix}_gt": "",
        f"{prefix}_qual": "",
        f"{prefix}_filter": "",
        f"{prefix}_dp": "",
        f"{prefix}_ad": "",
        f"{prefix}_af": "",
        f"{prefix}_gq": "",
    }


def platform_fields(
    prefix: str,
    variant: Variant,
) -> dict[str, str]:
    return {
        f"{prefix}_gt": variant.gt,
        f"{prefix}_qual": variant.qual,
        f"{prefix}_filter": variant.filt,
        f"{prefix}_dp": variant.dp,
        f"{prefix}_ad": variant.ad,
        f"{prefix}_af": variant.af,
        f"{prefix}_gq": variant.gq,
    }


def write_tsv(
    output_path: Path,
    rows: list[dict[str, str]],
) -> None:
    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    illumina_variants = query_vcf(args.illumina_vcf)
    ont_variants = query_vcf(args.ont_vcf)

    illumina_keys = set(illumina_variants)
    ont_keys = set(ont_variants)

    shared_keys = illumina_keys & ont_keys
    illumina_only_keys = illumina_keys - ont_keys
    ont_only_keys = ont_keys - illumina_keys

    illumina_only_variants = [
        illumina_variants[key]
        for key in sorted(illumina_only_keys)
    ]

    ont_depth = calculate_ont_depth(
        illumina_only_variants,
        args.ont_alignment,
        args.reference,
    )

    rows: list[dict[str, str]] = []

    for key in sorted(illumina_keys | ont_keys):
        illumina = illumina_variants.get(key)
        ont = ont_variants.get(key)

        chrom, pos, ref, alt = key

        row: dict[str, str] = {
            "sample_key": args.sample_key,
            "comparison_id": args.comparison_id,
            "chrom": chrom,
            "pos": str(pos),
            "ref": ref,
            "alt": alt,
            "variant_type": classify_variant(ref, alt),
            "ont_local_alignment_depth": "",
            "ont_depth_category": "",
            "illumina_callability": "",
        }

        if illumina is not None and ont is not None:
            if illumina.gt == ont.gt:
                category = "concordant_variant_and_genotype"
            else:
                category = "concordant_variant_discordant_genotype"

            row["comparison_category"] = category
            row.update(platform_fields("illumina", illumina))
            row.update(platform_fields("ont", ont))
            row["illumina_callability"] = "variant_record_present"

        elif illumina is not None:
            depth = ont_depth.get((chrom, pos), 0)

            if depth >= args.minimum_ont_depth:
                depth_category = (
                    f"ONT_depth_ge_{args.minimum_ont_depth}"
                )
            else:
                depth_category = (
                    f"ONT_depth_lt_{args.minimum_ont_depth}"
                )

            row["comparison_category"] = "illumina_only_record"
            row.update(platform_fields("illumina", illumina))
            row.update(blank_platform_fields("ont"))
            row["ont_local_alignment_depth"] = str(depth)
            row["ont_depth_category"] = depth_category
            row["illumina_callability"] = "variant_record_present"

        else:
            assert ont is not None

            row["comparison_category"] = "ont_only_record"
            row.update(blank_platform_fields("illumina"))
            row.update(platform_fields("ont", ont))
            row["illumina_callability"] = "unknown"

        rows.append(row)

    output_groups = {
        "variant_concordance_all.tsv": rows,
        "concordant_variants.tsv": [
            row
            for row in rows
            if row["comparison_category"]
            == "concordant_variant_and_genotype"
        ],
        "genotype_discordant_variants.tsv": [
            row
            for row in rows
            if row["comparison_category"]
            == "concordant_variant_discordant_genotype"
        ],
        "illumina_only_records.tsv": [
            row
            for row in rows
            if row["comparison_category"]
            == "illumina_only_record"
        ],
        "ont_only_records.tsv": [
            row
            for row in rows
            if row["comparison_category"]
            == "ont_only_record"
        ],
    }

    for filename, selected_rows in output_groups.items():
        write_tsv(
            args.output_dir / filename,
            selected_rows,
        )

    category_counts: dict[str, int] = {}

    for row in rows:
        category = row["comparison_category"]
        category_counts[category] = (
            category_counts.get(category, 0) + 1
        )

    concordant_variant_count = len(shared_keys)
    genotype_concordant_count = category_counts.get(
        "concordant_variant_and_genotype",
        0,
    )
    genotype_discordant_count = category_counts.get(
        "concordant_variant_discordant_genotype",
        0,
    )

    illumina_only_depth_ge = sum(
        1
        for row in rows
        if row["comparison_category"] == "illumina_only_record"
        and row["ont_depth_category"]
        == f"ONT_depth_ge_{args.minimum_ont_depth}"
    )

    illumina_only_depth_lt = sum(
        1
        for row in rows
        if row["comparison_category"] == "illumina_only_record"
        and row["ont_depth_category"]
        == f"ONT_depth_lt_{args.minimum_ont_depth}"
    )

    summary_fields = [
        "sample_key",
        "comparison_id",
        "illumina_prepared_records",
        "ont_prepared_records",
        "shared_variant_records",
        "concordant_variant_and_genotype",
        "concordant_variant_discordant_genotype",
        "illumina_only_records",
        f"illumina_only_ont_depth_ge_{args.minimum_ont_depth}",
        f"illumina_only_ont_depth_lt_{args.minimum_ont_depth}",
        "ont_only_records",
        "illumina_callability_for_ont_only",
        "record_overlap_fraction_illumina",
        "record_overlap_fraction_ont",
        "genotype_concordance_among_shared",
    ]

    def fraction(numerator: int, denominator: int) -> str:
        if denominator == 0:
            return ""

        return f"{numerator / denominator:.6f}"

    summary_row = {
        "sample_key": args.sample_key,
        "comparison_id": args.comparison_id,
        "illumina_prepared_records": str(len(illumina_keys)),
        "ont_prepared_records": str(len(ont_keys)),
        "shared_variant_records": str(concordant_variant_count),
        "concordant_variant_and_genotype": str(
            genotype_concordant_count
        ),
        "concordant_variant_discordant_genotype": str(
            genotype_discordant_count
        ),
        "illumina_only_records": str(len(illumina_only_keys)),
        f"illumina_only_ont_depth_ge_{args.minimum_ont_depth}": str(
            illumina_only_depth_ge
        ),
        f"illumina_only_ont_depth_lt_{args.minimum_ont_depth}": str(
            illumina_only_depth_lt
        ),
        "ont_only_records": str(len(ont_only_keys)),
        "illumina_callability_for_ont_only": "unknown",
        "record_overlap_fraction_illumina": fraction(
            concordant_variant_count,
            len(illumina_keys),
        ),
        "record_overlap_fraction_ont": fraction(
            concordant_variant_count,
            len(ont_keys),
        ),
        "genotype_concordance_among_shared": fraction(
            genotype_concordant_count,
            concordant_variant_count,
        ),
    }

    with (
        args.output_dir / "summary.tsv"
    ).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=summary_fields,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(summary_row)

    metadata_lines = [
        "Analysis type: descriptive cross-platform VCF-record concordance",
        "Truth set used: no",
        (
            "Comparison key: normalized CHROM, POS, REF and ALT"
        ),
        (
            "Genotype comparison: phase separators ignored; "
            "diploid allele order normalized"
        ),
        (
            "Comparison territory: variants overlapping "
            "autoinflammatory.bed"
        ),
        (
            "Illumina callability for ONT-only records: unknown "
            "because Illumina BAM, gVCF and callable-region files "
            "were unavailable"
        ),
        (
            "ONT local depth for Illumina-only records: calculated "
            "from the supplied ONT alignment"
        ),
        (
            "ONT depth excluded flags: 4, 256, 512, 1024 and 2048 "
            "(combined mask 3844)"
        ),
        (
            f"Operational ONT depth threshold: "
            f"{args.minimum_ont_depth}"
        ),
        f"Illumina prepared VCF: {args.illumina_vcf}",
        f"ONT prepared VCF: {args.ont_vcf}",
        f"ONT alignment: {args.ont_alignment}",
        f"Normalization reference: {args.reference}",
    ]

    (
        args.output_dir / "analysis_metadata.txt"
    ).write_text(
        "\n".join(metadata_lines) + "\n",
        encoding="utf-8",
    )

    print(
        f"{args.comparison_id}: "
        f"{len(illumina_keys)} Illumina, "
        f"{len(ont_keys)} ONT, "
        f"{concordant_variant_count} shared"
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
