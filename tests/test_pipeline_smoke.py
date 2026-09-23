"""Offline smoke tests for the leakage and terminology controls."""

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from condition_drug_pairs import ConditionDrugPairBuilder  # noqa: E402
from scalable_graph_builder import _fit_weights  # noqa: E402


def test_clinical_modifiers_and_conjunctions_are_preserved() -> None:
    builder = object.__new__(ConditionDrugPairBuilder)
    builder.include_placebo = False
    assert builder.normalize_for_match("metastatic breast cancer") == "metastatic breast cancer"
    assert builder.normalize_items(["hand, foot and mouth disease"]) == [
        "hand, foot and mouth disease"
    ]


def test_graph_weight_fit_marks_unlabelled_non_edges() -> None:
    edges = {
        (f"drug_{i}", f"disease_{j}")
        for i in range(5)
        for j in range(5)
        if i == j or (i + j) % 2 == 0
    }
    payload, training = _fit_weights(
        edges,
        target_pairs=set(),
        random_seed=7,
        max_training_positives=20,
    )
    assert set(training["LabelSemantics"]) == {"observed_edge", "unlabelled_non_edge"}
    assert "unlabelled" in payload["negative_sampling_semantics"]
    assert 0.0 <= payload["test_average_precision"] <= 1.0
    assert isinstance(training, pd.DataFrame)
