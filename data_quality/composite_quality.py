"""Composite evidence-readiness scoring rules.

The score measures auditability and evidence readiness, not clinical efficacy.
Weights are intentionally explicit so manuscript methods can cite the exact
formula used to generate the quality index.
"""

from __future__ import annotations

from typing import Mapping, Optional


COMPOSITE_QUALITY_WEIGHTS: Mapping[str, float] = {
    "entity_mapping": 0.15,
    "literature_completeness": 0.15,
    "semantic_relevance": 0.15,
    "adverse_cleanliness": 0.10,
    "safety_alignment": 0.15,
    "structural_consistency": 0.15,
    "posterior_certainty": 0.15,
}


def clamp01(value: Optional[float], default: float = 0.0) -> float:
    try:
        if value is None:
            value = default
        return max(0.0, min(1.0, float(value)))
    except Exception:
        return max(0.0, min(1.0, float(default)))


def evidence_readiness_score(
    *,
    entity_mapping: float,
    literature_completeness: float,
    semantic_relevance: float,
    adverse_cleanliness: float,
    safety_alignment: float,
    structural_consistency: Optional[float],
    posterior_certainty: float,
) -> float:
    """Return a 0-100 evidence-readiness score."""
    structural = 0.0 if structural_consistency is None else structural_consistency
    components = {
        "entity_mapping": clamp01(entity_mapping),
        "literature_completeness": clamp01(literature_completeness),
        "semantic_relevance": clamp01(semantic_relevance),
        "adverse_cleanliness": clamp01(adverse_cleanliness),
        "safety_alignment": clamp01(safety_alignment),
        "structural_consistency": clamp01(structural),
        "posterior_certainty": clamp01(posterior_certainty),
    }
    weighted = sum(components[key] * COMPOSITE_QUALITY_WEIGHTS[key] for key in components)
    return round(100.0 * clamp01(weighted), 3)

