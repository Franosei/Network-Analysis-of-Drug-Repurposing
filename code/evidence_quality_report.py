"""
Offline data-quality and evidence-quality report generator.

This module is intentionally additive: it reads the artifacts already produced by
the clinical-trial, graph, literature, safety, and Bayesian steps, then exports
auditable quality tables without re-running PubMed, OpenAI, or graph building.

Primary outputs:
  - pair_level_evidence_quality.csv
  - source_level_data_quality_audit.csv
  - summary_dashboard.csv
  - clinical_trial_data_readiness.csv
  - terminology_standardisation_quality.csv
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from data_quality.composite_quality import evidence_readiness_score
from data_quality.quality_flags import evidence_domain_support as rule_evidence_domain_support
from data_quality.quality_flags import quality_flag as rule_quality_flag

try:
    from scipy.integrate import trapezoid
    from scipy.special import rel_entr
    from scipy.stats import beta as beta_dist
except Exception:  # pragma: no cover - fallback for minimal environments
    trapezoid = None
    rel_entr = None
    beta_dist = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]

FEATURE_COLUMNS = [
    "GraphDistanceToIndication",
    "RandomWalkScore",
    "StructuralLikelihood",
    "PreferentialAttachment",
    "KatzSimilarity",
]

STRUCTURAL_WEIGHTS = {
    "GraphDistanceToIndication": 0.25,
    "RandomWalkScore": 0.25,
    "KatzSimilarity": 0.20,
    "StructuralLikelihood": 0.20,
    "PreferentialAttachment": 0.10,
}

VALID_LITERATURE_CATEGORIES = {"therapeutic", "adverse", "irrelevant"}


@dataclass(frozen=True)
class RunRecord:
    drug: str
    disease: str
    timestamp: str
    path: Path
    components: Dict[str, Any]


@dataclass(frozen=True)
class LiteratureFile:
    path: Path
    timestamp: str
    flat_name: str
    source_kind: str


def safe_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        if value is None or value == "":
            return default
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return default
        return out
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except Exception:
        return default


def safe_div(num: Any, den: Any, default: float = 0.0) -> float:
    n = safe_float(num, 0.0) or 0.0
    d = safe_float(den, 0.0) or 0.0
    if d == 0:
        return default
    return n / d


def clamp01(value: Any, default: float = 0.0) -> float:
    v = safe_float(value, default)
    if v is None:
        v = default
    return max(0.0, min(1.0, float(v)))


def round_or_none(value: Any, digits: int = 6) -> Optional[float]:
    v = safe_float(value, None)
    if v is None:
        return None
    return round(v, digits)


def norm_entity(text: Any) -> str:
    text = str(text or "").strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


def flat_key(text: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm_entity(text))


def pair_key(drug: Any, disease: Any) -> Tuple[str, str]:
    return norm_entity(drug), norm_entity(disease)


def read_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default


def is_placebo_like(text: Any) -> bool:
    return "placebo" in norm_entity(text)


def first_present(row: Dict[str, Any], keys: Iterable[str], default: Any = None) -> Any:
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return default


def timestamp_from_path(path: Path) -> str:
    match = re.search(r"(\d{8}_\d{6})", path.name)
    if match:
        return match.group(1)
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y%m%d_%H%M%S")


def count_input_trials(data_dir: Path) -> Optional[int]:
    if not data_dir.exists():
        return None

    total = 0
    found = False
    for path in data_dir.glob("*.json"):
        payload = read_json(path, None)
        if isinstance(payload, list):
            total += len(payload)
            found = True
        elif isinstance(payload, dict):
            studies = payload.get("studies") or payload.get("trials") or payload.get("data")
            if isinstance(studies, list):
                total += len(studies)
                found = True

    return total if found else None


def classify_unmatched_reason(entry: Dict[str, Any]) -> str:
    reason = norm_entity(entry.get("reason"))
    intervention = first_present(entry, ["raw_intervention", "intervention"])

    if "condition" in reason:
        return "unmatched_condition"
    if "intervention" in reason or "drug" in reason:
        return "unmatched_drug"
    if intervention in (None, ""):
        return "unmatched_condition"
    return "unmatched_unknown"


def canonical_mapping_method(method: Any, score: Any = None) -> str:
    raw = norm_entity(method)
    confidence = safe_float(score, None)

    if raw in {"", "none", "null"}:
        return "mapped_legacy_no_provenance"
    if raw in {"exact", "exact_match", "mesh_exact"}:
        return "exact_match"
    if raw in {"unmapped", "no_match", "empty", "no_tokens", "no_candidates", "token_score_failed"}:
        return "unmapped"
    if raw in {"fuzzy", "fuzzy_match", "token_score", "token-guided", "token_guided"}:
        if confidence is not None and confidence < 0.80:
            return "fuzzy_low_confidence"
        return "fuzzy_high_confidence"
    return raw.replace(" ", "_")


def mapping_confidence(method: str) -> float:
    return {
        "exact_match": 1.00,
        "fuzzy_high_confidence": 0.85,
        "fuzzy_low_confidence": 0.55,
        "mapped_legacy_no_provenance": 0.75,
        "unmapped": 0.00,
        "not_attempted": 0.00,
    }.get(method, 0.50)


def extract_mapping_method(row: Dict[str, Any], entity: str) -> str:
    method_keys = [
        f"{entity}_mapping_method",
        f"{entity}_match_method",
        f"{entity}_method",
        f"{entity}_mesh_method",
        f"{entity}_mapping",
    ]
    score_keys = [
        f"{entity}_mapping_score",
        f"{entity}_match_score",
        f"{entity}_score",
        f"{entity}_confidence",
    ]
    method = first_present(row, method_keys)
    score = first_present(row, score_keys)

    if method is None:
        debug = row.get("debug")
        if isinstance(debug, dict):
            method = debug.get("method")
            score = debug.get("score")

    return canonical_mapping_method(method, score)


def build_clinical_readiness(
    matched_rows: List[Dict[str, Any]],
    unmatched_rows: List[Dict[str, Any]],
    data_dir: Path,
) -> Dict[str, Any]:
    total_trials = count_input_trials(data_dir)
    matched_count = len(matched_rows)
    unmatched_count = len(unmatched_rows)
    pair_attempt_count = matched_count + unmatched_count

    unmatched_classes = Counter(classify_unmatched_reason(row) for row in unmatched_rows)
    unmatched_condition_count = unmatched_classes.get("unmatched_condition", 0)
    unresolved_drug_count = unmatched_classes.get("unmatched_drug", 0)

    matched_pairs = [
        pair_key(row.get("intervention"), row.get("condition"))
        for row in matched_rows
        if row.get("intervention") and row.get("condition")
    ]
    non_placebo_pairs = [p for p in matched_pairs if not is_placebo_like(p[0])]
    unique_non_placebo_pairs = set(non_placebo_pairs)
    placebo_count = sum(1 for drug, _ in matched_pairs if is_placebo_like(drug))
    duplicate_pair_count = max(len(non_placebo_pairs) - len(unique_non_placebo_pairs), 0)

    condition_resolved_rows = matched_count + unresolved_drug_count
    drug_resolution_denominator = matched_count + unresolved_drug_count

    return {
        "total_extracted_trial_records": total_trials,
        "drug_condition_rows": pair_attempt_count,
        "matched_drug_condition_rows": matched_count,
        "successfully_normalised_conditions": condition_resolved_rows,
        "successfully_normalised_drugs": matched_count,
        "unmatched_condition_count": unmatched_condition_count,
        "unmatched_condition_rate": safe_div(unmatched_condition_count, pair_attempt_count),
        "unresolved_drug_count": unresolved_drug_count,
        "unresolved_drug_rate": safe_div(unresolved_drug_count, drug_resolution_denominator),
        "placebo_exclusion_count": placebo_count,
        "duplicate_pair_count": duplicate_pair_count,
        "final_graph_ready_pair_count": len(unique_non_placebo_pairs),
        "legacy_note": (
            "total_extracted_trial_records is blank when raw ClinicalTrials.gov JSON files "
            "are not present; matched/unmatched rows are counted from processed_data."
        ),
    }


def build_mapping_quality_summary(
    matched_rows: List[Dict[str, Any]],
    unmatched_rows: List[Dict[str, Any]],
) -> pd.DataFrame:
    counters: Counter[Tuple[str, str]] = Counter()

    for row in matched_rows:
        condition_method = extract_mapping_method(row, "condition")
        drug_method = extract_mapping_method(row, "drug")
        if drug_method == "mapped_legacy_no_provenance":
            drug_method = extract_mapping_method(row, "intervention")
        counters[("condition", condition_method)] += 1
        counters[("drug", drug_method)] += 1

    for row in unmatched_rows:
        unmatched_class = classify_unmatched_reason(row)
        if unmatched_class == "unmatched_condition":
            counters[("condition", "unmapped")] += 1
            counters[("drug", "not_attempted")] += 1
        elif unmatched_class == "unmatched_drug":
            counters[("condition", "mapped_legacy_no_provenance")] += 1
            counters[("drug", "unmapped")] += 1
        else:
            counters[("condition", "unmapped")] += 1
            counters[("drug", "unmapped")] += 1

    rows: List[Dict[str, Any]] = []
    totals = Counter()
    for (entity_type, _), count in counters.items():
        totals[entity_type] += count

    for (entity_type, method), count in sorted(counters.items()):
        denominator = totals[entity_type]
        rows.append(
            {
                "entity_type": entity_type,
                "mapping_method": method,
                "count": count,
                "proportion": round(safe_div(count, denominator), 6),
                "nominal_confidence": mapping_confidence(method),
            }
        )

    return pd.DataFrame(rows)


def mapping_success_rate(mapping_df: pd.DataFrame) -> float:
    if mapping_df.empty:
        return 0.0
    total = float(mapping_df["count"].sum())
    unresolved = float(
        mapping_df.loc[
            mapping_df["mapping_method"].isin(["unmapped", "not_attempted"]),
            "count",
        ].sum()
    )
    return safe_div(total - unresolved, total)


def unresolved_entity_rate(mapping_df: pd.DataFrame) -> float:
    if mapping_df.empty:
        return 0.0
    total = float(mapping_df["count"].sum())
    unresolved = float(
        mapping_df.loc[
            mapping_df["mapping_method"].isin(["unmapped", "not_attempted"]),
            "count",
        ].sum()
    )
    return safe_div(unresolved, total)


def load_graph_features(known_path: Path, unknown_path: Path) -> pd.DataFrame:
    frames = []
    for path in [known_path, unknown_path]:
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()

    graph_df = pd.concat(frames, ignore_index=True)
    graph_df["Drug"] = graph_df["Drug"].astype(str).map(norm_entity)
    graph_df["Disease"] = graph_df["Disease"].astype(str).map(norm_entity)
    graph_df["pair_key"] = list(zip(graph_df["Drug"], graph_df["Disease"]))

    for col in FEATURE_COLUMNS:
        if col not in graph_df.columns:
            graph_df[col] = 0.0
        graph_df[col] = pd.to_numeric(graph_df[col], errors="coerce").fillna(0.0)
        graph_df[f"{col}_percentile"] = graph_df[col].rank(method="average", pct=True)

    graph_df["structural_consistency_score"] = 0.0
    for col, weight in STRUCTURAL_WEIGHTS.items():
        graph_df["structural_consistency_score"] += graph_df[f"{col}_percentile"] * weight

    graph_df["structural_consistency_score"] = graph_df["structural_consistency_score"].clip(0.0, 1.0)
    return graph_df


def load_latest_runs(runs_dir: Path) -> Dict[Tuple[str, str], RunRecord]:
    latest: Dict[Tuple[str, str], RunRecord] = {}
    if not runs_dir.exists():
        return latest

    for path in runs_dir.glob("run_*.json"):
        payload = read_json(path, {})
        if not isinstance(payload, dict):
            continue

        drug = str(payload.get("drug") or "").strip()
        disease = str(payload.get("disease") or "").strip()
        components = payload.get("components")
        if not drug or not disease or not isinstance(components, dict):
            continue

        timestamp = str(payload.get("timestamp") or timestamp_from_path(path))
        record = RunRecord(
            drug=drug,
            disease=disease,
            timestamp=timestamp,
            path=path,
            components=components,
        )
        key = pair_key(drug, disease)
        if key not in latest or record.timestamp > latest[key].timestamp:
            latest[key] = record

    return latest


def index_literature_files(literature_dir: Path) -> List[LiteratureFile]:
    out: List[LiteratureFile] = []
    if not literature_dir.exists():
        return out

    for path in literature_dir.glob("classified_*.json"):
        if path.name.startswith("classified_pubmed_"):
            source_kind = "pubmed"
            stem = path.stem.replace("classified_pubmed_", "", 1)
        elif path.name.startswith("classified_articles_"):
            source_kind = "article_archive"
            stem = path.stem.replace("classified_articles_", "", 1)
        else:
            continue

        timestamp = timestamp_from_path(path)
        stem_without_ts = re.sub(r"_\d{8}_\d{6}$", "", stem)
        out.append(
            LiteratureFile(
                path=path,
                timestamp=timestamp,
                flat_name=flat_key(stem_without_ts),
                source_kind=source_kind,
            )
        )

    return out


def find_literature_file(
    drug: str,
    disease: str,
    literature_files: List[LiteratureFile],
) -> Optional[LiteratureFile]:
    drug_flat = flat_key(drug)
    disease_flat = flat_key(disease)
    matches = [
        item
        for item in literature_files
        if drug_flat and disease_flat and drug_flat in item.flat_name and disease_flat in item.flat_name
    ]
    if not matches:
        return None

    matches.sort(key=lambda item: (item.source_kind == "pubmed", item.timestamp), reverse=True)
    return matches[0]


def literature_metrics(
    record: RunRecord,
    literature_files: List[LiteratureFile],
) -> Dict[str, Any]:
    file_record = find_literature_file(record.drug, record.disease, literature_files)
    components = record.components
    m_from_run = safe_int(components.get("M"), 0)

    if file_record is None:
        return {
            "literature_file": None,
            "literature_source_kind": None,
            "records_retrieved": m_from_run,
            "records_usable": m_from_run,
            "records_with_abstracts": m_from_run if m_from_run else 0,
            "records_with_intro": None,
            "records_with_conclusion": None,
            "records_excluded_missing_text": 0,
            "records_finally_classified": m_from_run,
            "retrieval_usability_rate": 1.0 if m_from_run else 0.0,
        }

    records = read_json(file_record.path, [])
    if not isinstance(records, list):
        records = []

    usable = len(records)
    with_abstract = sum(1 for row in records if str(row.get("abstract", "")).strip())
    with_intro = sum(1 for row in records if str(row.get("introduction", "")).strip())
    with_conclusion = sum(1 for row in records if str(row.get("conclusion", "")).strip())
    classified = sum(1 for row in records if norm_entity(row.get("category")) in VALID_LITERATURE_CATEGORIES)
    retrieved = max(m_from_run, usable)

    return {
        "literature_file": str(file_record.path.relative_to(PROJECT_ROOT)),
        "literature_source_kind": file_record.source_kind,
        "records_retrieved": retrieved,
        "records_usable": usable,
        "records_with_abstracts": with_abstract,
        "records_with_intro": with_intro,
        "records_with_conclusion": with_conclusion,
        "records_excluded_missing_text": max(retrieved - usable, 0),
        "records_finally_classified": classified,
        "retrieval_usability_rate": safe_div(usable, retrieved),
    }


def beta_metrics(components: Dict[str, Any]) -> Dict[str, Any]:
    prior_a = safe_float(components.get("prior_a"), None)
    prior_b = safe_float(components.get("prior_b"), None)
    post_a = safe_float(components.get("post_a"), None)
    post_b = safe_float(components.get("post_b"), None)

    if post_a is None or post_b is None or post_a <= 0 or post_b <= 0:
        return {
            "posterior_variance": None,
            "posterior_ci_low": None,
            "posterior_ci_high": None,
            "credible_interval_width": None,
            "kl_divergence": None,
            "prior_mean": None,
            "posterior_mean_integrated": None,
            "mean_shift": None,
        }

    total = post_a + post_b
    posterior_variance = (post_a * post_b) / ((total**2) * (total + 1.0))

    if beta_dist is None:
        mean = post_a / total
        sd = math.sqrt(posterior_variance)
        low = clamp01(mean - 1.96 * sd)
        high = clamp01(mean + 1.96 * sd)
        return {
            "posterior_variance": posterior_variance,
            "posterior_ci_low": low,
            "posterior_ci_high": high,
            "credible_interval_width": high - low,
            "kl_divergence": None,
            "prior_mean": safe_div(prior_a, (prior_a or 0.0) + (prior_b or 0.0), None),
            "posterior_mean_integrated": mean,
            "mean_shift": None,
        }

    low = float(beta_dist.ppf(0.025, post_a, post_b))
    high = float(beta_dist.ppf(0.975, post_a, post_b))

    kl = None
    prior_mean = None
    posterior_mean_integrated = post_a / total
    mean_shift = None
    if (
        prior_a is not None
        and prior_b is not None
        and prior_a > 0
        and prior_b > 0
        and trapezoid is not None
        and rel_entr is not None
    ):
        x = np.linspace(0.001, 0.999, 1000)
        prior_pdf = np.clip(beta_dist.pdf(x, prior_a, prior_b), 1e-12, None)
        post_pdf = np.clip(beta_dist.pdf(x, post_a, post_b), 1e-12, None)
        prior_pdf = prior_pdf / trapezoid(prior_pdf, x)
        post_pdf = post_pdf / trapezoid(post_pdf, x)
        kl = float(trapezoid(rel_entr(post_pdf, prior_pdf), x))
        prior_mean = float(trapezoid(x * prior_pdf, x))
        posterior_mean_integrated = float(trapezoid(x * post_pdf, x))
        mean_shift = posterior_mean_integrated - prior_mean

    return {
        "posterior_variance": posterior_variance,
        "posterior_ci_low": low,
        "posterior_ci_high": high,
        "credible_interval_width": high - low,
        "kl_divergence": kl,
        "prior_mean": prior_mean,
        "posterior_mean_integrated": posterior_mean_integrated,
        "mean_shift": mean_shift,
    }


def literature_completeness_score(m_articles: int, target_count: int = 30) -> float:
    if m_articles <= 0:
        return 0.0
    return clamp01(math.log1p(m_articles) / math.log1p(target_count))


def evidence_pattern(m_articles: int, therapeutic: int, adverse: int, irrelevant: int, gamma: float) -> str:
    if m_articles <= 0:
        return "no_literature_retrieved"

    therapeutic_rate = safe_div(therapeutic, m_articles)
    adverse_rate = safe_div(adverse, m_articles)
    irrelevant_rate = safe_div(irrelevant, m_articles)

    if m_articles < 3 and irrelevant_rate < 0.50:
        return "sparse_but_clean"
    if m_articles < 3:
        return "sparse_evidence"
    if irrelevant_rate >= 0.70:
        return "dominated_by_irrelevant_articles"
    if therapeutic_rate >= 0.25 and adverse_rate >= 0.20:
        return "abundant_but_conflicting"
    if adverse_rate >= 0.20 or gamma >= 0.70:
        return "safety_concerning"
    if therapeutic_rate >= 0.35 and irrelevant_rate < 0.50:
        return "therapeutically_relevant_literature"
    return "mixed_or_limited_signal"


def evidence_domain_support(
    m_articles: int,
    therapeutic: int,
    adverse: int,
    irrelevant: int,
    structural_score: Optional[float],
) -> str:
    return rule_evidence_domain_support(
        literature_count=m_articles,
        therapeutic_count=therapeutic,
        adverse_count=adverse,
        irrelevant_count=irrelevant,
        structural_score=structural_score,
    )


def quality_flag(
    readiness_score: float,
    entity_score: float,
    m_articles: int,
    adverse_rate: float,
    irrelevant_rate: float,
    gamma: float,
    structural_score: Optional[float],
    ci_width: Optional[float],
) -> str:
    return rule_quality_flag(
        readiness_score=readiness_score,
        entity_mapping_score=entity_score,
        literature_count=m_articles,
        adverse_rate=adverse_rate,
        irrelevant_rate=irrelevant_rate,
        safety_gamma=gamma,
        structural_score=structural_score,
        credible_interval_width=ci_width,
    )


def pair_mapping_quality(
    drug: str,
    disease: str,
    matched_rows: List[Dict[str, Any]],
    graph_df: pd.DataFrame,
) -> Tuple[float, str, str, str]:
    drug_norm, disease_norm = pair_key(drug, disease)
    matched_pairs = {
        pair_key(row.get("intervention"), row.get("condition"))
        for row in matched_rows
        if row.get("intervention") and row.get("condition")
    }

    if (drug_norm, disease_norm) in matched_pairs:
        return (
            0.75,
            "mapped_legacy_no_provenance",
            "mapped_legacy_no_provenance",
            "pair_seen_in_processed_trials",
        )

    if not graph_df.empty:
        graph_drugs = set(graph_df["Drug"].dropna().unique())
        graph_diseases = set(graph_df["Disease"].dropna().unique())
        if drug_norm in graph_drugs and disease_norm in graph_diseases:
            return (
                0.75,
                "mapped_legacy_no_provenance",
                "mapped_legacy_no_provenance",
                "entities_seen_in_graph_artifacts",
            )

    return (0.50, "not_captured", "not_captured", "not_found_in_local_mapping_artifacts")


def composite_readiness_score(
    entity_score: float,
    lit_completeness: float,
    semantic_relevance: float,
    adverse_cleanliness: float,
    safety_quality: float,
    structural_score: Optional[float],
    uncertainty_quality: float,
) -> float:
    return evidence_readiness_score(
        entity_mapping=entity_score,
        literature_completeness=lit_completeness,
        semantic_relevance=semantic_relevance,
        adverse_cleanliness=adverse_cleanliness,
        safety_alignment=safety_quality,
        structural_consistency=structural_score,
        posterior_certainty=uncertainty_quality,
    )


def build_pair_level_table(
    runs: Dict[Tuple[str, str], RunRecord],
    graph_df: pd.DataFrame,
    matched_rows: List[Dict[str, Any]],
    literature_files: List[LiteratureFile],
) -> pd.DataFrame:
    if not runs:
        return pd.DataFrame()

    graph_lookup: Dict[Tuple[str, str], Dict[str, Any]] = {}
    if not graph_df.empty:
        graph_lookup = (
            graph_df.drop_duplicates("pair_key", keep="first")
            .set_index("pair_key")
            .to_dict(orient="index")
        )

    rows: List[Dict[str, Any]] = []
    for key, run in sorted(runs.items(), key=lambda item: (item[1].drug.lower(), item[1].disease.lower())):
        components = run.components
        counts = components.get("counts") if isinstance(components.get("counts"), dict) else {}
        therapeutic = safe_int(counts.get("therapeutic"), 0)
        adverse = safe_int(counts.get("adverse"), 0)
        irrelevant = safe_int(counts.get("irrelevant"), 0)
        m_articles = safe_int(components.get("M"), therapeutic + adverse + irrelevant)
        if m_articles == 0:
            m_articles = therapeutic + adverse + irrelevant

        gamma_value = safe_float(components.get("gamma"), None)
        gamma_for_scoring = clamp01(gamma_value, 0.5) if gamma_value is not None else 0.5

        lit_metrics = literature_metrics(run, literature_files)
        uncertainty = beta_metrics(components)

        graph_row = graph_lookup.get(key, {})
        structural_score = safe_float(graph_row.get("structural_consistency_score"), None)
        entity_score, drug_method, disease_method, terminology_flag = pair_mapping_quality(
            run.drug,
            run.disease,
            matched_rows,
            graph_df,
        )

        therapeutic_rate = safe_div(therapeutic, m_articles)
        adverse_rate = safe_div(adverse, m_articles)
        irrelevant_rate = safe_div(irrelevant, m_articles)
        lit_score = literature_completeness_score(m_articles)
        semantic_relevance = clamp01(1.0 - irrelevant_rate)
        adverse_cleanliness = clamp01(1.0 - adverse_rate)
        safety_quality = clamp01(1.0 - gamma_for_scoring)
        ci_width = uncertainty.get("credible_interval_width")
        uncertainty_quality = clamp01(1.0 - (ci_width if ci_width is not None else 1.0))

        readiness = composite_readiness_score(
            entity_score=entity_score,
            lit_completeness=lit_score,
            semantic_relevance=semantic_relevance,
            adverse_cleanliness=adverse_cleanliness,
            safety_quality=safety_quality,
            structural_score=structural_score,
            uncertainty_quality=uncertainty_quality,
        )

        matching_effects = components.get("matching_effects")
        if not isinstance(matching_effects, list):
            matching_effects = []

        rows.append(
            {
                "drug": run.drug,
                "disease": run.disease,
                "run_timestamp": run.timestamp,
                "run_log": str(run.path.relative_to(PROJECT_ROOT)),
                "drug_mapping_method": drug_method,
                "disease_mapping_method": disease_method,
                "entity_mapping_quality_score": entity_score,
                "terminology_quality_flag": terminology_flag,
                "records_retrieved": lit_metrics["records_retrieved"],
                "records_usable": lit_metrics["records_usable"],
                "records_with_abstracts": lit_metrics["records_with_abstracts"],
                "records_with_intro": lit_metrics["records_with_intro"],
                "records_with_conclusion": lit_metrics["records_with_conclusion"],
                "records_excluded_missing_text": lit_metrics["records_excluded_missing_text"],
                "records_finally_classified": lit_metrics["records_finally_classified"],
                "retrieval_usability_rate": round_or_none(lit_metrics["retrieval_usability_rate"]),
                "literature_file": lit_metrics["literature_file"],
                "therapeutic_count": therapeutic,
                "adverse_count": adverse,
                "irrelevant_count": irrelevant,
                "literature_completeness_count": m_articles,
                "therapeutic_relevance_ratio": round_or_none(therapeutic_rate),
                "adverse_evidence_burden": round_or_none(adverse_rate),
                "irrelevant_retrieval_noise_rate": round_or_none(irrelevant_rate),
                "literature_completeness_score": round_or_none(lit_score),
                "literature_evidence_pattern": evidence_pattern(
                    m_articles,
                    therapeutic,
                    adverse,
                    irrelevant,
                    gamma_for_scoring,
                ),
                "gamma_safety_overlap": round_or_none(gamma_value),
                "safety_conflict_score": round_or_none(gamma_value),
                "safety_overlap_term_count": len(matching_effects) if matching_effects else None,
                "top_safety_overlap_terms": "; ".join(map(str, matching_effects[:10])) if matching_effects else None,
                "GraphDistanceToIndication": round_or_none(graph_row.get("GraphDistanceToIndication")),
                "RandomWalkScore": round_or_none(graph_row.get("RandomWalkScore")),
                "StructuralLikelihood": round_or_none(graph_row.get("StructuralLikelihood")),
                "PreferentialAttachment": round_or_none(graph_row.get("PreferentialAttachment")),
                "KatzSimilarity": round_or_none(graph_row.get("KatzSimilarity")),
                "structural_consistency_score": round_or_none(structural_score),
                "evidence_domain_support": evidence_domain_support(
                    m_articles,
                    therapeutic,
                    adverse,
                    irrelevant,
                    structural_score,
                ),
                "p_raw": round_or_none(components.get("p_raw")),
                "p_penalised": round_or_none(components.get("p_penalised")),
                "p_final": round_or_none(components.get("p_final")),
                "cM": round_or_none(components.get("cM")),
                "p_likelihood": round_or_none(components.get("p_likelihood")),
                "posterior_mean": round_or_none(components.get("post_mean")),
                "posterior_variance": round_or_none(uncertainty.get("posterior_variance"), 8),
                "posterior_ci_low": round_or_none(uncertainty.get("posterior_ci_low")),
                "posterior_ci_high": round_or_none(uncertainty.get("posterior_ci_high")),
                "credible_interval_width": round_or_none(ci_width),
                "kl_divergence": round_or_none(uncertainty.get("kl_divergence")),
                "mean_shift": round_or_none(uncertainty.get("mean_shift")),
                "evidence_readiness_score": readiness,
                "quality_flag": quality_flag(
                    readiness,
                    entity_score,
                    m_articles,
                    adverse_rate,
                    irrelevant_rate,
                    gamma_for_scoring,
                    structural_score,
                    ci_width,
                ),
            }
        )

    return pd.DataFrame(rows)


def add_audit_row(
    rows: List[Dict[str, Any]],
    source: str,
    dimension: str,
    metric: str,
    value: Any,
    denominator: Any = None,
    rate: Any = None,
    note: str = "",
) -> None:
    rows.append(
        {
            "source": source,
            "quality_dimension": dimension,
            "metric": metric,
            "value": value,
            "denominator": denominator,
            "rate": round_or_none(rate),
            "note": note,
        }
    )


def build_source_level_audit(
    clinical: Dict[str, Any],
    mapping_df: pd.DataFrame,
    pair_df: pd.DataFrame,
    graph_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    for metric in [
        "total_extracted_trial_records",
        "drug_condition_rows",
        "matched_drug_condition_rows",
        "successfully_normalised_conditions",
        "successfully_normalised_drugs",
        "unmatched_condition_count",
        "unresolved_drug_count",
        "placebo_exclusion_count",
        "duplicate_pair_count",
        "final_graph_ready_pair_count",
    ]:
        add_audit_row(rows, "ClinicalTrials.gov", "data_readiness", metric, clinical.get(metric))

    add_audit_row(
        rows,
        "ClinicalTrials.gov",
        "data_readiness",
        "unmatched_condition_rate",
        clinical.get("unmatched_condition_count"),
        clinical.get("drug_condition_rows"),
        clinical.get("unmatched_condition_rate"),
    )
    add_audit_row(
        rows,
        "ClinicalTrials.gov",
        "data_readiness",
        "unresolved_drug_rate",
        clinical.get("unresolved_drug_count"),
        clinical.get("matched_drug_condition_rows", 0) + clinical.get("unresolved_drug_count", 0),
        clinical.get("unresolved_drug_rate"),
    )

    for _, row in mapping_df.iterrows():
        add_audit_row(
            rows,
            "MeSH",
            "terminology_standardisation",
            f"{row['entity_type']}:{row['mapping_method']}",
            int(row["count"]),
            None,
            row["proportion"],
            f"nominal_confidence={row['nominal_confidence']}",
        )

    if not pair_df.empty:
        add_audit_row(rows, "PubMed/PMC", "retrieval_quality", "audited_pairs", len(pair_df))
        add_audit_row(rows, "PubMed/PMC", "retrieval_quality", "records_retrieved", int(pair_df["records_retrieved"].sum()))
        add_audit_row(rows, "PubMed/PMC", "retrieval_quality", "records_usable", int(pair_df["records_usable"].sum()))
        add_audit_row(
            rows,
            "PubMed/PMC",
            "retrieval_quality",
            "usable_literature_rate",
            int(pair_df["records_usable"].sum()),
            int(pair_df["records_retrieved"].sum()),
            safe_div(pair_df["records_usable"].sum(), pair_df["records_retrieved"].sum()),
        )
        add_audit_row(
            rows,
            "PubMed/PMC",
            "evidence_quality",
            "average_irrelevant_retrieval_noise_rate",
            pair_df["irrelevant_retrieval_noise_rate"].mean(),
        )
        add_audit_row(
            rows,
            "PubMed/PMC",
            "evidence_quality",
            "average_therapeutic_relevance_ratio",
            pair_df["therapeutic_relevance_ratio"].mean(),
        )

        gamma_series = pd.to_numeric(pair_df["gamma_safety_overlap"], errors="coerce")
        add_audit_row(rows, "openFDA/FAERS", "safety_alignment", "average_safety_conflict_score", gamma_series.mean())
        add_audit_row(rows, "openFDA/FAERS", "safety_alignment", "pairs_gamma_ge_0_50", int((gamma_series >= 0.50).sum()), len(pair_df), safe_div((gamma_series >= 0.50).sum(), len(pair_df)))
        add_audit_row(rows, "openFDA/FAERS", "safety_alignment", "pairs_missing_gamma", int(gamma_series.isna().sum()), len(pair_df), safe_div(gamma_series.isna().sum(), len(pair_df)))

        add_audit_row(rows, "Bayesian posterior", "uncertainty_quality", "average_posterior_mean", pair_df["posterior_mean"].mean())
        add_audit_row(rows, "Bayesian posterior", "uncertainty_quality", "average_posterior_variance", pair_df["posterior_variance"].mean())
        add_audit_row(rows, "Bayesian posterior", "uncertainty_quality", "average_credible_interval_width", pair_df["credible_interval_width"].mean())
        add_audit_row(rows, "Bayesian posterior", "uncertainty_quality", "average_kl_information_gain", pair_df["kl_divergence"].mean())

    if not graph_df.empty:
        add_audit_row(rows, "Graph", "structural_consistency", "feature_rows", len(graph_df))
        add_audit_row(rows, "Graph", "structural_consistency", "known_pair_rows", int((graph_df["Label"] == 1).sum()) if "Label" in graph_df else None)
        add_audit_row(rows, "Graph", "structural_consistency", "unknown_pair_rows", int((graph_df["Label"] == 0).sum()) if "Label" in graph_df else None)
        add_audit_row(rows, "Graph", "structural_consistency", "unique_drugs", graph_df["Drug"].nunique())
        add_audit_row(rows, "Graph", "structural_consistency", "unique_diseases", graph_df["Disease"].nunique())
        add_audit_row(rows, "Graph", "structural_consistency", "average_structural_consistency_score", graph_df["structural_consistency_score"].mean())

    return pd.DataFrame(rows)


def build_summary_dashboard(
    clinical: Dict[str, Any],
    mapping_df: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    def add(metric: str, value: Any, note: str = "") -> None:
        rows.append({"metric": metric, "value": value, "note": note})

    add("report_generated_at", datetime.now().isoformat(timespec="seconds"))
    add("pair_count", 0 if pair_df.empty else len(pair_df))
    add("clinical_drug_condition_rows", clinical.get("drug_condition_rows"))
    add("final_graph_ready_pair_count", clinical.get("final_graph_ready_pair_count"))
    add("mapping_success_rate", round(mapping_success_rate(mapping_df), 6))
    add("unresolved_entity_rate", round(unresolved_entity_rate(mapping_df), 6))

    if not pair_df.empty:
        add("usable_literature_rate", round(safe_div(pair_df["records_usable"].sum(), pair_df["records_retrieved"].sum()), 6))
        add("average_irrelevant_retrieval_rate", round(pair_df["irrelevant_retrieval_noise_rate"].mean(), 6))
        add("average_safety_overlap_score", round(pair_df["gamma_safety_overlap"].mean(), 6))
        add("average_posterior_uncertainty", round(pair_df["credible_interval_width"].mean(), 6), "Mean 95% credible-interval width.")
        add("average_evidence_readiness_score", round(pair_df["evidence_readiness_score"].mean(), 3))
        add("median_evidence_readiness_score", round(pair_df["evidence_readiness_score"].median(), 3))

        for label, count in pair_df["quality_flag"].value_counts().sort_index().items():
            add(f"quality_category_count:{label}", int(count))
        for support, count in pair_df["evidence_domain_support"].value_counts().sort_index().items():
            add(f"evidence_domain_support_count:{support}", int(count))

    add("legacy_mapping_provenance_note", "Current matched pair artifacts do not store exact vs fuzzy match provenance for successful mappings.")
    return pd.DataFrame(rows)


def write_markdown_report(
    output_dir: Path,
    clinical: Dict[str, Any],
    pair_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> None:
    lines = [
        "# Data Quality and Evidence Quality Report",
        "",
        "This report was generated offline from existing pipeline artifacts. It measures evidence readiness and data quality, not clinical efficacy.",
        "",
        "## Outputs",
        "",
        "- `pair_level_evidence_quality.csv`: one row per audited drug-disease run log.",
        "- `source_level_data_quality_audit.csv`: ClinicalTrials.gov, MeSH, PubMed/PMC, safety, graph, and posterior quality metrics.",
        "- `summary_dashboard.csv`: high-level counts, rates, and quality category totals.",
        "- `clinical_trial_data_readiness.csv`: extraction/normalisation readiness metrics.",
        "- `terminology_standardisation_quality.csv`: mapping-method proportions and nominal confidence.",
        "",
        "## Clinical-Trial Readiness",
        "",
        f"- Drug-condition rows audited: {clinical.get('drug_condition_rows')}",
        f"- Matched rows: {clinical.get('matched_drug_condition_rows')}",
        f"- Unmatched condition rate: {round_or_none(clinical.get('unmatched_condition_rate'))}",
        f"- Unresolved drug rate: {round_or_none(clinical.get('unresolved_drug_rate'))}",
        f"- Final graph-ready unique non-placebo pairs: {clinical.get('final_graph_ready_pair_count')}",
        "",
        "## Pair-Level Evidence",
        "",
    ]

    if pair_df.empty:
        lines.append("No Bayesian run logs were found for pair-level reporting.")
    else:
        top = pair_df.sort_values("evidence_readiness_score", ascending=False).head(10)
        lines.append(f"- Audited pair count: {len(pair_df)}")
        lines.append(f"- Average evidence readiness score: {pair_df['evidence_readiness_score'].mean():.3f}")
        lines.append(f"- Average posterior uncertainty width: {pair_df['credible_interval_width'].mean():.6f}")
        lines.append("")
        lines.append("Top evidence-readiness pairs:")
        for _, row in top.iterrows():
            lines.append(
                f"- {row['drug']} -> {row['disease']}: "
                f"{row['evidence_readiness_score']:.3f} ({row['quality_flag']})"
            )

    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- Evidence readiness scores combine mapping reliability, literature completeness, semantic relevance, adverse burden, safety overlap, structural consistency, and posterior uncertainty.",
            "- Successful MeSH mappings in the current legacy processed pair file do not include exact/fuzzy provenance, so they are labelled `mapped_legacy_no_provenance`.",
            "- Safety overlap terms are included when present in run logs; older logs may only contain gamma.",
            "",
            "## Dashboard Snapshot",
            "",
        ]
    )
    for _, row in summary_df.iterrows():
        lines.append(f"- {row['metric']}: {row['value']}")

    (output_dir / "data_quality_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_reports(args: argparse.Namespace) -> Dict[str, Path]:
    matched_path = PROJECT_ROOT / args.matched
    unmatched_path = PROJECT_ROOT / args.unmatched
    known_graph_path = PROJECT_ROOT / args.known_graph
    unknown_graph_path = PROJECT_ROOT / args.unknown_graph
    runs_dir = PROJECT_ROOT / args.runs_dir
    literature_dir = PROJECT_ROOT / args.literature_dir
    data_dir = PROJECT_ROOT / args.data_dir
    output_dir = PROJECT_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    matched_rows = read_json(matched_path, [])
    unmatched_rows = read_json(unmatched_path, [])
    if not isinstance(matched_rows, list):
        matched_rows = []
    if not isinstance(unmatched_rows, list):
        unmatched_rows = []

    clinical = build_clinical_readiness(matched_rows, unmatched_rows, data_dir)
    mapping_df = build_mapping_quality_summary(matched_rows, unmatched_rows)
    graph_df = load_graph_features(known_graph_path, unknown_graph_path)
    runs = load_latest_runs(runs_dir)
    literature_files = index_literature_files(literature_dir)
    pair_df = build_pair_level_table(runs, graph_df, matched_rows, literature_files)
    source_df = build_source_level_audit(clinical, mapping_df, pair_df, graph_df)
    summary_df = build_summary_dashboard(clinical, mapping_df, pair_df)

    clinical_df = pd.DataFrame([clinical])
    category_df = (
        pair_df["quality_flag"].value_counts().rename_axis("quality_flag").reset_index(name="count")
        if not pair_df.empty
        else pd.DataFrame(columns=["quality_flag", "count"])
    )

    paths = {
        "pair_level": output_dir / "pair_level_evidence_quality.csv",
        "source_level": output_dir / "source_level_data_quality_audit.csv",
        "summary": output_dir / "summary_dashboard.csv",
        "clinical_readiness": output_dir / "clinical_trial_data_readiness.csv",
        "terminology": output_dir / "terminology_standardisation_quality.csv",
        "quality_counts": output_dir / "quality_category_counts.csv",
        "markdown": output_dir / "data_quality_report.md",
    }

    pair_df.to_csv(paths["pair_level"], index=False)
    source_df.to_csv(paths["source_level"], index=False)
    summary_df.to_csv(paths["summary"], index=False)
    clinical_df.to_csv(paths["clinical_readiness"], index=False)
    mapping_df.to_csv(paths["terminology"], index=False)
    category_df.to_csv(paths["quality_counts"], index=False)
    write_markdown_report(output_dir, clinical, pair_df, summary_df)

    return paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate offline data-quality and evidence-quality tables from existing pipeline artifacts."
    )
    parser.add_argument("--matched", default="processed_data/condition_drug_pairs.json")
    parser.add_argument("--unmatched", default="processed_data/unmatched_pairs.json")
    parser.add_argument("--known-graph", default="graph/graph_features_known.csv")
    parser.add_argument("--unknown-graph", default="graph/graph_features_unknown.csv")
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument("--literature-dir", default="literatures")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--output-dir", default="reports/data_quality")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = generate_reports(args)
    print("Generated data-quality outputs:")
    for label, path in paths.items():
        print(f"- {label}: {path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
