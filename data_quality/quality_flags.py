"""Rule-based quality flags and coverage tiers."""

from __future__ import annotations

from typing import Optional


def coverage_tier(
    *,
    has_mapping: bool,
    has_graph: bool,
    has_literature: bool,
    has_safety: bool,
    has_bayesian: bool,
) -> str:
    if has_mapping and has_graph and has_literature and has_safety and has_bayesian:
        return "full_bayesian_audit"
    if has_graph and has_literature and has_bayesian:
        return "bayesian_without_safety"
    if has_graph and has_literature:
        return "literature_and_graph"
    if has_graph:
        return "graph_only"
    if has_literature:
        return "literature_only"
    if has_mapping:
        return "matched_pairs_only"
    return "insufficient_coverage"


def evidence_domain_support(
    *,
    literature_count: int,
    therapeutic_count: int,
    adverse_count: int,
    irrelevant_count: int,
    structural_score: Optional[float],
) -> str:
    literature_supported = (
        literature_count >= 5
        and literature_count > 0
        and irrelevant_count / literature_count < 0.70
        and (therapeutic_count + adverse_count) > 0
    )
    network_supported = (structural_score or 0.0) >= 0.60

    if literature_supported and network_supported:
        return "literature_and_network"
    if literature_supported:
        return "literature_only"
    if network_supported:
        return "network_only"
    return "neither"


def quality_flag(
    *,
    readiness_score: float,
    entity_mapping_score: float,
    literature_count: int,
    adverse_rate: float,
    irrelevant_rate: float,
    safety_gamma: Optional[float],
    structural_score: Optional[float],
    credible_interval_width: Optional[float],
) -> str:
    gamma = 0.5 if safety_gamma is None else safety_gamma
    uncertainty_width = 1.0 if credible_interval_width is None else credible_interval_width
    structural = 0.0 if structural_score is None else structural_score

    if entity_mapping_score < 0.60:
        return "Terminology uncertainty"
    if literature_count < 3 and structural >= 0.60:
        return "Sparse literature; structurally plausible"
    if literature_count < 3:
        return "Insufficient evidence"
    if adverse_rate >= 0.20 and gamma >= 0.50:
        return "Safety-conflicted evidence"
    if gamma >= 0.75:
        return "Safety-concerning"
    if irrelevant_rate >= 0.70:
        return "Literature noise dominated"
    if adverse_rate >= 0.20:
        return "Literature-conflicted"
    if readiness_score >= 75.0 and uncertainty_width <= 0.25:
        return "High evidence quality"
    if readiness_score >= 50.0:
        return "Moderate evidence quality"
    return "Low evidence quality"
