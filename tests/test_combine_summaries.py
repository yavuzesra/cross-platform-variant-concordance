from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def test_combine_summaries(tmp_path: Path) -> None:
    results = tmp_path / "results"
    for name, value in [("pair_a", "10"), ("pair_b", "20")]:
        folder = results / name
        folder.mkdir(parents=True)
        with (folder / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["comparison_id", "shared_variant_records"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"comparison_id": name, "shared_variant_records": value})

    output = tmp_path / "combined.tsv"
    script = Path(__file__).resolve().parents[1] / "scripts" / "combine_summaries.py"
    result = subprocess.run(
        [sys.executable, str(script), "--results-dir", str(results), "--output", str(output)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    text = output.read_text(encoding="utf-8")
    assert "pair_a" in text
    assert "pair_b" in text
