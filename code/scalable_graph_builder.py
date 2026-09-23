"""Scalable, leakage-controlled graph features for full-registry reruns.

The legacy graph builder materialises every drug x disease non-edge and a dense
Katz inverse.  That is appropriate only for a small sampled graph.  This module
keeps the full observed MeSH graph, computes features only for requested target
pairs, and re-fits graph weights from deterministic held-out edges/non-edges.
"""

from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


FEATURE_NAMES = [
    "GraphDistanceToIndication",
    "RandomWalkScore",
    "StructuralLikelihood",
    "PreferentialAttachment",
    "KatzSimilarity",
]


def _norm(value: Any) -> str:
    return str(value or "").strip().casefold()


def _adjacency(
    edges: Iterable[tuple[str, str]],
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    drug_to_diseases: dict[str, set[str]] = defaultdict(set)
    disease_to_drugs: dict[str, set[str]] = defaultdict(set)
    for drug, disease in edges:
        drug_to_diseases[drug].add(disease)
        disease_to_drugs[disease].add(drug)
    return dict(drug_to_diseases), dict(disease_to_drugs)


def _pair_features(
    drug: str,
    disease: str,
    drug_to_diseases: dict[str, set[str]],
    disease_to_drugs: dict[str, set[str]],
    *,
    katz_beta: float = 0.005,
) -> dict[str, float]:
    source_diseases = drug_to_diseases.get(drug, set())
    target_drugs = disease_to_drugs.get(disease, set())
    drug_degree = len(source_diseases)
    disease_degree = len(target_drugs)

    path_count = 0
    random_walk = 0.0
    if drug_degree and disease_degree:
        if len(source_diseases) <= len(target_drugs):
            for middle_disease in source_diseases:
                middle_drugs = disease_to_drugs.get(middle_disease, set())
                connecting = middle_drugs.intersection(target_drugs)
                path_count += len(connecting)
                for middle_drug in connecting:
                    denominator = (
                        drug_degree
                        * max(1, len(middle_drugs))
                        * max(1, len(drug_to_diseases.get(middle_drug, set())))
                    )
                    random_walk += 1.0 / denominator
        else:
            for middle_drug in target_drugs:
                middle_diseases = drug_to_diseases.get(middle_drug, set())
                connecting = middle_diseases.intersection(source_diseases)
                path_count += len(connecting)
                middle_drug_degree = max(1, len(middle_diseases))
                for middle_disease in connecting:
                    denominator = (
                        drug_degree
                        * max(1, len(disease_to_drugs.get(middle_disease, set())))
                        * middle_drug_degree
                    )
                    random_walk += 1.0 / denominator

    n_drugs = max(2, len(drug_to_diseases))
    n_diseases = max(2, len(disease_to_drugs))
    drug_centrality = drug_degree / (n_drugs + n_diseases - 1)
    disease_centrality = disease_degree / (n_drugs + n_diseases - 1)
    return {
        # In a bipartite graph, a withheld pair's first alternative route is
        # length 3.  Legacy distance-to-another-indication is then 1/(1+2).
        "GraphDistanceToIndication": (1.0 / 3.0) if path_count else 0.0,
        "RandomWalkScore": random_walk,
        "StructuralLikelihood": (1.0 + drug_centrality) * (1.0 + disease_centrality),
        "PreferentialAttachment": drug_centrality * disease_centrality,
        "KatzSimilarity": (katz_beta**3) * path_count,
        "AlternativePathCountLength3": float(path_count),
        "DrugDegree": float(drug_degree),
        "DiseaseDegree": float(disease_degree),
    }


def _sample_non_edges(
    drugs: list[str],
    diseases: list[str],
    observed: set[tuple[str, str]],
    count: int,
    rng: random.Random,
) -> list[tuple[str, str]]:
    output: set[tuple[str, str]] = set()
    attempts = 0
    max_attempts = max(10_000, count * 50)
    while len(output) < count and attempts < max_attempts:
        pair = (rng.choice(drugs), rng.choice(diseases))
        attempts += 1
        if pair not in observed:
            output.add(pair)
    if len(output) != count:
        raise RuntimeError(f"Could sample only {len(output):,}/{count:,} graph non-edges")
    return sorted(output)


def _fit_weights(
    edges: set[tuple[str, str]],
    *,
    target_pairs: set[tuple[str, str]],
    random_seed: int,
    max_training_positives: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    rng = random.Random(random_seed)
    eligible_edges = sorted(edges - target_pairs)
    rng.shuffle(eligible_edges)
    holdout_count = min(max_training_positives, max(1, int(0.20 * len(eligible_edges))))
    positive_examples = eligible_edges[:holdout_count]
    base_edges = set(eligible_edges[holdout_count:])
    drug_to_diseases, disease_to_drugs = _adjacency(base_edges)

    drugs = sorted(drug_to_diseases)
    diseases = sorted(disease_to_drugs)
    negative_examples = _sample_non_edges(
        drugs, diseases, edges, len(positive_examples), rng
    )

    rows = []
    for label, pairs in ((1, positive_examples), (0, negative_examples)):
        for drug, disease in pairs:
            features = _pair_features(
                drug, disease, drug_to_diseases, disease_to_drugs
            )
            rows.append({"Drug": drug, "Disease": disease, "Label": label, **features})
    dataset = pd.DataFrame(rows)
    x = dataset[FEATURE_NAMES].astype(float)
    y = dataset["Label"].astype(int)
    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=0.30,
        random_state=random_seed,
        stratify=y,
    )
    scaler = StandardScaler()
    x_train_scaled = scaler.fit_transform(x_train)
    x_test_scaled = scaler.transform(x_test)
    model = LogisticRegression(
        max_iter=2000,
        class_weight="balanced",
        random_state=random_seed,
    )
    model.fit(x_train_scaled, y_train)
    probability = model.predict_proba(x_test_scaled)[:, 1]
    auc = float(roc_auc_score(y_test, probability))

    coef_scaled = model.coef_[0]
    coef_raw = coef_scaled / scaler.scale_
    intercept_raw = float(
        model.intercept_[0] - np.sum(coef_scaled * scaler.mean_ / scaler.scale_)
    )
    weights = {name: float(value) for name, value in zip(FEATURE_NAMES, coef_raw)}
    payload = {
        "method": "balanced logistic regression on deterministic 20% held-out registered edges and equal sampled non-edges",
        "leakage_control": "held-out positive edges are absent from the graph used to compute their features; requested target pairs are excluded from fitting",
        "random_seed": random_seed,
        "base_graph_edges": len(base_edges),
        "training_positive_examples": len(positive_examples),
        "training_negative_examples": len(negative_examples),
        "test_roc_auc": auc,
        "feature_names": FEATURE_NAMES,
        "scaled_coefficients": {
            name: float(value) for name, value in zip(FEATURE_NAMES, coef_scaled)
        },
        "scaler_mean": {
            name: float(value) for name, value in zip(FEATURE_NAMES, scaler.mean_)
        },
        "scaler_scale": {
            name: float(value) for name, value in zip(FEATURE_NAMES, scaler.scale_)
        },
        "raw_feature_weights": weights,
        "raw_feature_intercept": intercept_raw,
    }
    return payload, dataset


def build_targeted_full_graph(
    *,
    matched_path: Path,
    graph_dir: Path,
    target_pairs: list[tuple[str, str]],
    random_seed: int = 42,
    max_training_positives: int = 20_000,
) -> tuple[pd.DataFrame, Path, Path, Path]:
    graph_dir.mkdir(parents=True, exist_ok=True)
    records = json.loads(matched_path.read_text(encoding="utf-8"))
    edges: set[tuple[str, str]] = set()
    pair_ncts: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in records:
        pair = (_norm(row.get("intervention")), _norm(row.get("condition")))
        if not all(pair) or "placebo" in pair[0]:
            continue
        edges.add(pair)
        nct_id = str(row.get("nct_id") or "").strip()
        if nct_id:
            pair_ncts[pair].add(nct_id)

    targets = {(_norm(drug), _norm(disease)) for drug, disease in target_pairs}
    weights_payload, training_df = _fit_weights(
        edges,
        target_pairs=targets,
        random_seed=random_seed,
        max_training_positives=max_training_positives,
    )
    weights_path = graph_dir / "updated_graph_weights.json"
    weights_path.write_text(json.dumps(weights_payload, indent=2), encoding="utf-8")
    training_df.to_csv(graph_dir / "graph_weight_training_sample.csv", index=False)

    withheld_edges = edges - targets
    drug_to_diseases, disease_to_drugs = _adjacency(withheld_edges)
    known_rows = []
    unknown_rows = []
    intercept = float(weights_payload["raw_feature_intercept"])
    weights = weights_payload["raw_feature_weights"]
    for pair in sorted(targets):
        features = _pair_features(
            pair[0], pair[1], drug_to_diseases, disease_to_drugs
        )
        score = intercept + sum(weights[name] * features[name] for name in FEATURE_NAMES)
        probability = 1.0 / (1.0 + math.exp(-max(-700.0, min(700.0, score))))
        row = {
            "Drug": pair[0],
            "Disease": pair[1],
            **{name: features[name] for name in FEATURE_NAMES},
            "AlternativePathCountLength3": int(features["AlternativePathCountLength3"]),
            "DrugDegree": int(features["DrugDegree"]),
            "DiseaseDegree": int(features["DiseaseDegree"]),
            "GraphLogit": score,
            "GraphProbability": probability,
            "TargetEdgeWithheld": True,
            "CanonicalTrialCount": len(pair_ncts.get(pair, set())),
            "Label": int(pair in edges),
        }
        (known_rows if pair in edges else unknown_rows).append(row)

    known_path = graph_dir / "graph_features_known.csv"
    unknown_path = graph_dir / "graph_features_unknown.csv"
    columns = [
        "Drug",
        "Disease",
        *FEATURE_NAMES,
        "AlternativePathCountLength3",
        "DrugDegree",
        "DiseaseDegree",
        "GraphLogit",
        "GraphProbability",
        "TargetEdgeWithheld",
        "CanonicalTrialCount",
        "Label",
    ]
    pd.DataFrame(known_rows, columns=columns).to_csv(known_path, index=False)
    pd.DataFrame(unknown_rows, columns=columns).to_csv(unknown_path, index=False)

    graph = nx.Graph()
    graph.add_nodes_from(drug_to_diseases, bipartite="drug")
    graph.add_nodes_from(disease_to_drugs, bipartite="disease")
    graph.add_edges_from(withheld_edges)
    nx.write_graphml(graph, graph_dir / "bipartite_target_edge_withheld.graphml")

    audit = pd.DataFrame(
        [
            {"metric": "canonical_observed_edges", "value": len(edges)},
            {"metric": "withheld_graph_edges", "value": len(withheld_edges)},
            {"metric": "unique_drugs", "value": len(drug_to_diseases)},
            {"metric": "unique_diseases", "value": len(disease_to_drugs)},
            {"metric": "weight_training_roc_auc", "value": weights_payload["test_roc_auc"]},
            {"metric": "target_pairs_scored", "value": len(targets)},
            {"metric": "target_edges_withheld", "value": len(edges.intersection(targets))},
        ]
    )
    return audit, known_path, unknown_path, weights_path

