from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "compare_small_variants.py"
spec = importlib.util.spec_from_file_location("compare_small_variants", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)


def test_clean_missing_values() -> None:
    assert module.clean(".") == ""
    assert module.clean("  ") == ""
    assert module.clean("PASS") == "PASS"


def test_normalize_genotype_ignores_phasing_and_order() -> None:
    assert module.normalize_genotype("0|1") == "0/1"
    assert module.normalize_genotype("1|0") == "0/1"
    assert module.normalize_genotype("1/1") == "1/1"


def test_classify_variant() -> None:
    assert module.classify_variant("A", "G") == "SNV"
    assert module.classify_variant("A", "AT") == "insertion"
    assert module.classify_variant("AT", "A") == "deletion"
    assert module.classify_variant("AT", "GC") == "complex_substitution"
