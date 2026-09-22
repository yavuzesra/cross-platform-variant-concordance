#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine per-comparison summary TSV files."
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()

    summary_files = sorted(
        path
        for path in args.results_dir.glob("*/summary.tsv")
        if path.is_file()
    )

    if not summary_files:
        print(
            "ERROR: No per-comparison summary.tsv files found.",
            file=sys.stderr,
        )
        return 1

    combined_rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None

    for summary_file in summary_files:
        with summary_file.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            if reader.fieldnames is None:
                raise RuntimeError(
                    f"Missing header: {summary_file}"
                )

            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise RuntimeError(
                    f"Header mismatch: {summary_file}"
                )

            combined_rows.extend(reader)

    assert fieldnames is not None

    args.output.parent.mkdir(parents=True, exist_ok=True)

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
        writer.writerows(combined_rows)

    print(
        f"Combined {len(combined_rows)} comparison summaries: "
        f"{args.output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
