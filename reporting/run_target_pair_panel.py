"""Run a focused target-pair evidence-readiness panel.

This script reuses the archived full-run artifacts when a requested pair already
has a Bayesian run log. Missing pairs are refreshed through the same PubMed,
LLM, safety, graph, and posterior calculations used by the main pipeline.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = PROJECT_ROOT / "code"
for path in (PROJECT_ROOT, CODE_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from data_quality.quality_flags import coverage_tier  # noqa: E402
from evidence_quality_report import (  # noqa: E402
    beta_metrics,
    clamp01,
    composite_readiness_score,
    evidence_domain_support,
    literature_completeness_score,
    load_graph_features,
    load_latest_runs,
    norm_entity,
    pair_key,
    pair_mapping_quality,
    quality_flag,
    read_json,
    round_or_none,
    safe_div,
    safe_float,
    safe_int,
)


TARGET_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("Dexamethasone", "COVID-19"),
    ("Remdesivir", "COVID-19"),
    ("Molnupiravir", "COVID-19"),
    ("Baricitinib", "COVID-19"),
    ("BIO101", "COVID-19"),
    ("Valproic acid", "pancreatic cancer"),
    ("Simvastatin", "pancreatic cancer"),
    ("Disulfiram", "SARS-CoV-2"),
    ("Carfilzomib", "SARS-CoV-2"),
    ("Amodiaquine", "SARS-CoV-2"),
    ("Sildenafil", "erectile dysfunction"),
    ("Thalidomide", "multiple myeloma"),
)

FEATURE_COLUMNS: Tuple[str, ...] = (
    "GraphDistanceToIndication",
    "RandomWalkScore",
    "StructuralLikelihood",
    "PreferentialAttachment",
    "KatzSimilarity",
)


@dataclass(frozen=True)
class PanelRecord:
    drug: str
    disease: str
    timestamp: str
    path: Path
    components: Dict[str, Any]
    source: str


def safe_run_slug(text: str) -> str:
    return re.sub(r"[^a-z0-9_,.-]+", "_", text.strip().lower()).strip("_")


def flat_key(text: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm_entity(text))


def disease_flat_variants(text: Any) -> List[str]:
    base = flat_key(text)
    variants = [base]
    raw = norm_entity(text)
    for part in re.split(r"\s*(?:/|\bor\b|\band\b)\s*", raw):
        if part:
            variants.append(flat_key(part))
    covid_aliases = {"covid19", "sarscov2", "sarscov2infection", "coronavirusdisease2019"}
    if base in covid_aliases:
        variants.extend(sorted(covid_aliases))
    if base == "pancreaticcancer":
        variants.append("pancreaticneoplasms")
    if "hairloss" in base or "alopecia" in base:
        variants.extend(["alopecia", "androgeneticalopecia", "alopeciaandrogenetic"])
    if "leukaemia" in base:
        variants.append(base.replace("leukaemia", "leukemia"))
    if "myelodysplasticsyndromes" in base:
        variants.append("myelodysplasticsyndromes")
    if "acutemyeloid" in base:
        variants.extend(["leukemiamyeloidacute", "acutemyeloidleukemia", "acutemyeloidleukaemia"])
    if "pancreaticductaladenocarcinoma" in base or "pancreaticadenocarcinoma" in base or "pancreatic" in base:
        variants.extend(["pancreaticneoplasms", "pancreaticcancer", "pancreaticductaladenocarcinoma"])
    return list(dict.fromkeys(item for item in variants if item))


def rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def load_run_config(base_run_dir: Path) -> Dict[str, Any]:
    return read_json(base_run_dir / "run_config.json", {})


def get_run_params(run_config: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    search = run_config.get("search_parameters") if isinstance(run_config.get("search_parameters"), dict) else {}
    llm = run_config.get("llm") if isinstance(run_config.get("llm"), dict) else {}
    bayes = run_config.get("bayesian") if isinstance(run_config.get("bayesian"), dict) else {}
    safety = (
        run_config.get("safety_penalty_settings")
        if isinstance(run_config.get("safety_penalty_settings"), dict)
        else {}
    )

    return {
        "pubmed_max_articles": args.pubmed_max_articles
        if args.pubmed_max_articles is not None
        else int(search.get("pubmed_max_articles", 40)),
        "pubmed_filter_level": args.pubmed_filter_level or str(search.get("pubmed_filter_level", "high")),
        "pubmed_years_back": args.pubmed_years_back
        if args.pubmed_years_back is not None
        else int(search.get("pubmed_years_back", 10)),
        "force_pubmed_refresh": args.force_pubmed_refresh,
        "llm_model": args.llm_model or str(llm.get("model", "gpt-4o-mini")),
        "llm_batch_size": args.llm_batch_size if args.llm_batch_size is not None else int(llm.get("batch_size", 5)),
        "llm_delay_s": args.llm_delay_s if args.llm_delay_s is not None else float(llm.get("delay_s", 2.0)),
        "llm_max_retries": args.llm_max_retries
        if args.llm_max_retries is not None
        else int(llm.get("max_retries", 2)),
        "cmax": args.cmax if args.cmax is not None else float(bayes.get("cmax", 200.0)),
        "tau": args.tau if args.tau is not None else float(bayes.get("tau", 25.0)),
        "likelihood_strength": args.likelihood_strength
        if args.likelihood_strength is not None
        else float(bayes.get("likelihood_strength", 50.0)),
        "likelihood_intercept": args.likelihood_intercept
        if args.likelihood_intercept is not None
        else float(bayes.get("likelihood_intercept", 0.0)),
        "weights": bayes.get("weights") if isinstance(bayes.get("weights"), dict) else {},
        "safety_penalty_scale": float(safety.get("penalty_scale", 0.5)),
    }


def graph_lookup(graph_df: pd.DataFrame) -> Dict[Tuple[str, str], Dict[str, Any]]:
    if graph_df.empty:
        return {}
    return (
        graph_df.drop_duplicates("pair_key", keep="first")
        .set_index("pair_key")
        .to_dict(orient="index")
    )


def load_matched_rows(base_run_dir: Path) -> List[Dict[str, Any]]:
    payload = read_json(base_run_dir / "processed_data" / "condition_drug_pairs.json", [])
    return payload if isinstance(payload, list) else []


def resolve_canonical_pair(
    drug: str,
    disease: str,
    matched_rows: List[Dict[str, Any]],
) -> Tuple[str, str, str]:
    """Resolve requested labels to the processed ClinicalTrials.gov terms."""
    drug_norm, disease_norm = pair_key(drug, disease)
    exact_matches = [
        row
        for row in matched_rows
        if norm_entity(row.get("intervention")) == drug_norm and norm_entity(row.get("condition")) == disease_norm
    ]
    if exact_matches:
        row = exact_matches[0]
        return str(row.get("intervention") or drug), str(row.get("condition") or disease), "exact_processed_pair"

    drug_flat = flat_key(drug)
    disease_flats = disease_flat_variants(disease)
    approximate: List[Tuple[str, str]] = []
    drug_candidates: List[str] = []
    disease_candidates: List[str] = []
    for row in matched_rows:
        intervention_text = " | ".join(
            str(row.get(key, "")) for key in ["intervention", "raw_intervention", "intervention_cleaned_term"]
        )
        condition_text = " | ".join(
            str(row.get(key, "")) for key in ["condition", "raw_condition", "condition_cleaned_term"]
        )
        intervention_flat = flat_key(intervention_text)
        condition_flat = flat_key(condition_text)
        if drug_flat in intervention_flat and any(disease_flat in condition_flat for disease_flat in disease_flats):
            canonical_drug = str(row.get("intervention") or drug).strip()
            canonical_disease = str(row.get("condition") or disease).strip()
            if canonical_drug and canonical_disease:
                approximate.append((canonical_drug, canonical_disease))
        if drug_flat in intervention_flat:
            canonical_drug = str(row.get("intervention") or "").strip()
            if canonical_drug:
                drug_candidates.append(canonical_drug)
        if any(disease_flat in condition_flat for disease_flat in disease_flats):
            canonical_disease = str(row.get("condition") or "").strip()
            if canonical_disease:
                disease_candidates.append(canonical_disease)

    if approximate:
        (canonical_drug, canonical_disease), _ = Counter(approximate).most_common(1)[0]
        return canonical_drug, canonical_disease, "approximate_processed_pair"

    if drug_candidates and disease_candidates:
        canonical_drug, _ = Counter(drug_candidates).most_common(1)[0]
        canonical_disease, _ = Counter(disease_candidates).most_common(1)[0]
        return canonical_drug, canonical_disease, "entity_level_processed_terms"

    return drug, disease, "literal_request"


def matching_run_record(
    drug: str,
    disease: str,
    available_runs: Dict[Tuple[str, str], Any],
    source_lookup: Dict[Tuple[str, str], str],
) -> Optional[PanelRecord]:
    key = pair_key(drug, disease)
    record = available_runs.get(key)
    if record is None:
        return None
    return PanelRecord(
        drug=drug,
        disease=disease,
        timestamp=record.timestamp,
        path=record.path,
        components=record.components,
        source=source_lookup.get(key, "existing_run"),
    )


def copy_archived_run(record: PanelRecord, output_runs_dir: Path) -> Path:
    output_runs_dir.mkdir(parents=True, exist_ok=True)
    destination = output_runs_dir / record.path.name
    if destination.resolve() != record.path.resolve():
        shutil.copy2(record.path, destination)
    return destination


def sigmoid(value: float) -> float:
    if value >= 0:
        exp_value = math.exp(-value)
        return 1.0 / (1.0 + exp_value)
    exp_value = math.exp(value)
    return exp_value / (1.0 + exp_value)


def beta_params_from_prob(probability: float, strength: float) -> Tuple[float, float]:
    p = min(max(float(probability), 1e-6), 1.0 - 1e-6)
    strength = max(float(strength), 0.0)
    return 1.0 + strength * p, 1.0 + strength * (1.0 - p)


def components_with_graph(
    components: Dict[str, Any],
    graph_row: Dict[str, Any],
    params: Dict[str, Any],
) -> Dict[str, Any]:
    if not graph_row:
        return dict(components)

    out = dict(components)
    graph_score = float(params["likelihood_intercept"])
    feature_values: Dict[str, float] = {}
    weights = params.get("weights") if isinstance(params.get("weights"), dict) else {}
    for feature in FEATURE_COLUMNS:
        value = float(graph_row.get(feature, 0.0) or 0.0)
        feature_values[feature] = value
        graph_score += float(weights.get(feature, 0.0)) * value

    p_likelihood = min(max(sigmoid(graph_score), 1e-6), 1.0 - 1e-6)
    like_a, like_b = beta_params_from_prob(p_likelihood, float(params["likelihood_strength"]))

    prior_a = safe_float(out.get("prior_a"), None)
    prior_b = safe_float(out.get("prior_b"), None)
    if prior_a is None or prior_b is None:
        return out

    post_a = prior_a + (like_a - 1.0)
    post_b = prior_b + (like_b - 1.0)
    out.update(
        {
            "graph_score": graph_score,
            "p_likelihood": p_likelihood,
            "likelihood_strength": float(params["likelihood_strength"]),
            "like_a": like_a,
            "like_b": like_b,
            "feature_values": feature_values,
            "post_a": post_a,
            "post_b": post_b,
            "post_mean": post_a / (post_a + post_b),
        }
    )
    return out


def make_classifier(params: Dict[str, Any], output_dir: Path) -> Any:
    from pubmed_utils import LLMClassifier, LLMConfig, PubMedSearchConfig

    search_cfg = PubMedSearchConfig(
        max_results=int(params["pubmed_max_articles"]),
        years_back=int(params["pubmed_years_back"]),
        cache_dir=str(output_dir / "cache"),
    )
    llm_cfg = LLMConfig(
        model=str(params["llm_model"]),
        delay_s=float(params["llm_delay_s"]),
        batch_size=int(params["llm_batch_size"]),
        max_retries=int(params["llm_max_retries"]),
    )
    return LLMClassifier(search_cfg=search_cfg, llm_cfg=llm_cfg)


def compute_fresh_record(
    drug: str,
    disease: str,
    classifier: Any,
    params: Dict[str, Any],
    graph_row: Dict[str, Any],
    output_dir: Path,
    timestamp: str,
) -> PanelRecord:
    from bayesian_predictor import beta_params_from_prob, clamp01 as bayes_clamp01
    from bayesian_predictor import concentration_c, sigmoid

    literature_dir = output_dir / "literatures"
    runs_dir = output_dir / "runs"
    literature_dir.mkdir(parents=True, exist_ok=True)
    runs_dir.mkdir(parents=True, exist_ok=True)

    prior_result = classifier.build_semantic_prior(
        drug=drug,
        disease=disease,
        max_articles=int(params["pubmed_max_articles"]),
        filter_level=str(params["pubmed_filter_level"]),
        save_dir=str(literature_dir),
        use_cache=not bool(params["force_pubmed_refresh"]),
    )

    counts = dict(prior_result.get("raw_counts", {}))
    m_articles = int(prior_result.get("total_articles", 0))
    p_raw = bayes_clamp01(prior_result.get("prior", 0.5))
    p_pen = bayes_clamp01(prior_result.get("penalised_prior", 0.5))
    p_final = bayes_clamp01(prior_result.get("enhanced_prior", 0.5))
    gamma = prior_result.get("gamma")
    gamma = bayes_clamp01(gamma) if gamma is not None else None

    c_m = concentration_c(m_articles, cmax=float(params["cmax"]), tau=float(params["tau"]))
    prior_a, prior_b = beta_params_from_prob(p_final, c_m)

    graph_score = float(params["likelihood_intercept"])
    feature_values: Dict[str, float] = {}
    weights = params.get("weights") if isinstance(params.get("weights"), dict) else {}
    if graph_row:
        for feature in FEATURE_COLUMNS:
            value = float(graph_row.get(feature, 0.0) or 0.0)
            feature_values[feature] = value
            graph_score += float(weights.get(feature, 0.0)) * value

    p_like = bayes_clamp01(sigmoid(graph_score))
    like_a, like_b = beta_params_from_prob(p_like, float(params["likelihood_strength"]))
    post_a = prior_a + (like_a - 1.0)
    post_b = prior_b + (like_b - 1.0)
    post_mean = post_a / (post_a + post_b)

    components = {
        "p_raw": p_raw,
        "p_penalised": p_pen,
        "p_final": p_final,
        "gamma": gamma,
        "counts": counts,
        "M": m_articles,
        "cM": c_m,
        "prior_a": prior_a,
        "prior_b": prior_b,
        "graph_score": graph_score,
        "p_likelihood": p_like,
        "likelihood_strength": float(params["likelihood_strength"]),
        "like_a": like_a,
        "like_b": like_b,
        "feature_values": feature_values,
        "post_a": post_a,
        "post_b": post_b,
        "post_mean": float(post_mean),
        "matching_effects": prior_result.get("matching_effects", []),
        "safety_relation": prior_result.get("safety_relation"),
        "safety_penalty_scale": prior_result.get("safety_penalty_scale", params["safety_penalty_scale"]),
    }
    payload = {
        "timestamp": timestamp,
        "drug": drug,
        "disease": disease,
        "case_type": "target_pair_panel",
        "config": {
            "pubmed_max_articles": int(params["pubmed_max_articles"]),
            "pubmed_filter_level": str(params["pubmed_filter_level"]),
            "pubmed_years_back": int(params["pubmed_years_back"]),
            "pubmed_use_cache": not bool(params["force_pubmed_refresh"]),
            "llm_model": str(params["llm_model"]),
            "llm_batch_size": int(params["llm_batch_size"]),
            "llm_delay_s": float(params["llm_delay_s"]),
            "llm_max_retries": int(params["llm_max_retries"]),
            "cmax": float(params["cmax"]),
            "tau": float(params["tau"]),
            "likelihood_strength": float(params["likelihood_strength"]),
            "likelihood_intercept": float(params["likelihood_intercept"]),
            "weights": weights,
        },
        "components": components,
    }
    run_path = runs_dir / f"run_{safe_run_slug(drug)}_{safe_run_slug(disease)}_{timestamp}.json"
    run_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return PanelRecord(drug=drug, disease=disease, timestamp=timestamp, path=run_path, components=components, source="fresh_target_run")


def format_literature_composition(therapeutic: int, adverse: int, irrelevant: int, m_articles: int) -> str:
    if m_articles <= 0:
        return "M=0; no classified literature retrieved"
    return (
        f"M={m_articles}; therapeutic={therapeutic} ({safe_div(therapeutic, m_articles):.2f}), "
        f"adverse={adverse} ({safe_div(adverse, m_articles):.2f}), "
        f"irrelevant={irrelevant} ({safe_div(irrelevant, m_articles):.2f})"
    )


def format_safety_signal(gamma: Optional[float], matching_effects: Any) -> str:
    if gamma is None:
        return "gamma=NA; no safety-overlap score returned"
    terms = matching_effects if isinstance(matching_effects, list) else []
    if terms:
        shown_terms = ", ".join(str(item) for item in terms[:4])
        return f"gamma={gamma:.3f}; overlap terms={shown_terms}"
    return f"gamma={gamma:.3f}; no named overlap terms"


def format_graph_support(
    graph_row: Dict[str, Any],
    components: Dict[str, Any],
    requested_drug: str,
    requested_disease: str,
    canonical_drug: str,
    canonical_disease: str,
) -> str:
    p_like = safe_float(components.get("p_likelihood"), None)
    if not graph_row:
        p_text = "NA" if p_like is None else f"{p_like:.3f}"
        return f"no exact graph row; neutral likelihood p={p_text}"
    structural = safe_float(graph_row.get("structural_consistency_score"), None)
    label = safe_int(graph_row.get("Label"), 0)
    label_text = "known" if label == 1 else "candidate"
    requested = pair_key(requested_drug, requested_disease)
    canonical = pair_key(canonical_drug, canonical_disease)
    canonical_text = ""
    if requested != canonical:
        canonical_text = f"; canonical={canonical_drug}-{canonical_disease}"
    return (
        f"{label_text}{canonical_text}; structural={structural:.3f}; "
        f"p_likelihood={safe_float(components.get('p_likelihood'), 0.5):.3f}"
    )


def format_posterior(components: Dict[str, Any], uncertainty: Dict[str, Any]) -> str:
    mean = safe_float(components.get("post_mean"), None)
    low = safe_float(uncertainty.get("posterior_ci_low"), None)
    high = safe_float(uncertainty.get("posterior_ci_high"), None)
    if mean is None:
        return "posterior unavailable"
    if low is None or high is None:
        return f"mean={mean:.3f}; CI=NA"
    return f"mean={mean:.3f}; 95% CrI={low:.3f}-{high:.3f}"


def interpret_flag(flag: str, readiness: float, domain_support: str, posterior_mean: Optional[float]) -> str:
    if flag == "High evidence quality":
        return "Evidence-ready for clinical-context review; uncertainty is well constrained."
    if flag == "Moderate evidence quality":
        return f"Moderate audit readiness with {domain_support.replace('_', ' ')} support."
    if flag == "Safety-conflicted evidence":
        return "Therapeutic signal is offset by adverse/safety overlap; review safety evidence first."
    if flag == "Safety-concerning":
        return "Safety overlap is high enough to limit prioritisation."
    if flag == "Literature noise dominated":
        return "Retrieved literature is dominated by irrelevant evidence; search/mapping needs review."
    if flag == "Terminology uncertainty":
        return "Local terminology or graph mapping is incomplete; interpret posterior cautiously."
    if flag == "Insufficient evidence":
        return "Evidence base is too sparse for a reliable readiness judgement."
    if flag == "Sparse literature; structurally plausible":
        return "Sparse literature, but network structure provides supporting context."
    if readiness < 50:
        return "Low audit readiness; additional evidence or mapping validation is needed."
    if posterior_mean is not None and posterior_mean >= 0.60:
        return "Posterior is elevated, but readiness depends on supporting evidence quality."
    return "Mixed evidence; use as a triage result rather than a clinical conclusion."


def build_panel_row(
    record: PanelRecord,
    graph_row: Dict[str, Any],
    matched_rows: List[Dict[str, Any]],
    graph_df: pd.DataFrame,
    canonical_drug: str,
    canonical_disease: str,
    canonical_resolution: str,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    components = components_with_graph(record.components, graph_row, params)
    counts = components.get("counts") if isinstance(components.get("counts"), dict) else {}
    therapeutic = safe_int(counts.get("therapeutic"), 0)
    adverse = safe_int(counts.get("adverse"), 0)
    irrelevant = safe_int(counts.get("irrelevant"), 0)
    m_articles = safe_int(components.get("M"), therapeutic + adverse + irrelevant)
    if m_articles == 0:
        m_articles = therapeutic + adverse + irrelevant

    gamma_value = safe_float(components.get("gamma"), None)
    gamma_for_scoring = clamp01(gamma_value, 0.5) if gamma_value is not None else 0.5
    uncertainty = beta_metrics(components)
    ci_width = safe_float(uncertainty.get("credible_interval_width"), None)
    structural_score = safe_float(graph_row.get("structural_consistency_score"), None) if graph_row else None

    entity_score, drug_method, disease_method, terminology_flag = pair_mapping_quality(
        canonical_drug,
        canonical_disease,
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
    domain_support = evidence_domain_support(m_articles, therapeutic, adverse, irrelevant, structural_score)
    flag = quality_flag(
        readiness,
        entity_score,
        m_articles,
        adverse_rate,
        irrelevant_rate,
        gamma_for_scoring,
        structural_score,
        ci_width,
    )
    tier = coverage_tier(
        has_mapping=entity_score >= 0.60,
        has_graph=bool(graph_row),
        has_literature=m_articles > 0,
        has_safety=gamma_value is not None,
        has_bayesian=safe_float(components.get("post_mean"), None) is not None,
    )
    matching_effects = components.get("matching_effects")
    posterior_mean = safe_float(components.get("post_mean"), None)

    return {
        "drug_disease_pair": f"{record.drug}-{record.disease}",
        "coverage_tier": tier,
        "literature_composition": format_literature_composition(therapeutic, adverse, irrelevant, m_articles),
        "safety_signal": format_safety_signal(gamma_value, matching_effects),
        "graph_support": format_graph_support(
            graph_row,
            components,
            record.drug,
            record.disease,
            canonical_drug,
            canonical_disease,
        ),
        "bayesian_posterior": format_posterior(components, uncertainty),
        "evidence_readiness_score": round(readiness, 3),
        "quality_flag": flag,
        "interpretation": interpret_flag(flag, readiness, domain_support, posterior_mean),
        "drug": record.drug,
        "disease": record.disease,
        "record_source": record.source,
        "canonical_drug": canonical_drug,
        "canonical_disease": canonical_disease,
        "canonical_resolution": canonical_resolution,
        "run_log": rel_path(record.path),
        "run_timestamp": record.timestamp,
        "drug_mapping_method": drug_method,
        "disease_mapping_method": disease_method,
        "entity_mapping_quality_score": round_or_none(entity_score),
        "terminology_quality_flag": terminology_flag,
        "therapeutic_count": therapeutic,
        "adverse_count": adverse,
        "irrelevant_count": irrelevant,
        "literature_count": m_articles,
        "therapeutic_rate": round_or_none(therapeutic_rate),
        "adverse_rate": round_or_none(adverse_rate),
        "irrelevant_rate": round_or_none(irrelevant_rate),
        "gamma_safety_overlap": round_or_none(gamma_value),
        "GraphDistanceToIndication": round_or_none(graph_row.get("GraphDistanceToIndication")) if graph_row else None,
        "RandomWalkScore": round_or_none(graph_row.get("RandomWalkScore")) if graph_row else None,
        "StructuralLikelihood": round_or_none(graph_row.get("StructuralLikelihood")) if graph_row else None,
        "PreferentialAttachment": round_or_none(graph_row.get("PreferentialAttachment")) if graph_row else None,
        "KatzSimilarity": round_or_none(graph_row.get("KatzSimilarity")) if graph_row else None,
        "structural_consistency_score": round_or_none(structural_score),
        "evidence_domain_support": domain_support,
        "p_raw": round_or_none(components.get("p_raw")),
        "p_penalised": round_or_none(components.get("p_penalised")),
        "p_final": round_or_none(components.get("p_final")),
        "cM": round_or_none(components.get("cM")),
        "graph_score": round_or_none(components.get("graph_score")),
        "p_likelihood": round_or_none(components.get("p_likelihood")),
        "posterior_mean": round_or_none(posterior_mean),
        "posterior_ci_low": round_or_none(uncertainty.get("posterior_ci_low")),
        "posterior_ci_high": round_or_none(uncertainty.get("posterior_ci_high")),
        "credible_interval_width": round_or_none(ci_width),
        "kl_divergence": round_or_none(uncertainty.get("kl_divergence")),
        "mean_shift": round_or_none(uncertainty.get("mean_shift")),
    }


def parse_pair_line(line: str) -> Optional[Tuple[str, str]]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    for sep in ["\t", ",", "|"]:
        if sep in line:
            left, right = line.split(sep, 1)
            return left.strip(), right.strip()
    for sep in ["--", " - ", "->"]:
        if sep in line:
            left, right = line.split(sep, 1)
            return left.strip(), right.strip()
    raise ValueError(f"Could not parse pair line: {line}")


def load_target_pairs(path: Optional[Path]) -> List[Tuple[str, str]]:
    if path is None:
        return list(TARGET_PAIRS)
    pairs: List[Tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        parsed = parse_pair_line(line)
        if parsed is not None:
            pairs.append(parsed)
    return pairs


def markdown_cell(value: Any) -> str:
    text = "" if value is None else str(value)
    text = text.replace("\n", " ").replace("|", "\\|")
    return text


def dataframe_to_markdown(df: pd.DataFrame, columns: Iterable[str]) -> str:
    selected = list(columns)
    header = "| " + " | ".join(selected) + " |"
    separator = "| " + " | ".join("---" for _ in selected) + " |"
    body = []
    for _, row in df.loc[:, selected].iterrows():
        body.append("| " + " | ".join(markdown_cell(row[col]) for col in selected) + " |")
    return "\n".join([header, separator, *body])


def write_markdown_table(df: pd.DataFrame, path: Path, columns: Iterable[str]) -> None:
    lines = [
        f"# Target Pair Evidence-Readiness Panel ({date.today().isoformat()})",
        "",
        dataframe_to_markdown(df, columns),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run-dir", type=Path, default=PROJECT_ROOT / "outputs" / "20260610_bayesian")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / "target_pair_panel_20260617")
    parser.add_argument("--pairs-file", type=Path, default=None)
    parser.add_argument("--no-reuse-archived", action="store_true")
    parser.add_argument("--force-pubmed-refresh", action="store_true")
    parser.add_argument("--pubmed-max-articles", type=int, default=None)
    parser.add_argument("--pubmed-filter-level", type=str, default=None)
    parser.add_argument("--pubmed-years-back", type=int, default=None)
    parser.add_argument("--llm-model", type=str, default=None)
    parser.add_argument("--llm-batch-size", type=int, default=None)
    parser.add_argument("--llm-delay-s", type=float, default=None)
    parser.add_argument("--llm-max-retries", type=int, default=None)
    parser.add_argument("--cmax", type=float, default=None)
    parser.add_argument("--tau", type=float, default=None)
    parser.add_argument("--likelihood-strength", type=float, default=None)
    parser.add_argument("--likelihood-intercept", type=float, default=None)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    base_run_dir = args.base_run_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "runs").mkdir(parents=True, exist_ok=True)

    target_pairs = load_target_pairs(args.pairs_file)
    run_config = load_run_config(base_run_dir)
    params = get_run_params(run_config, args)

    graph_df = load_graph_features(
        base_run_dir / "graph" / "graph_features_known.csv",
        base_run_dir / "graph" / "graph_features_unknown.csv",
    )
    graph_rows = graph_lookup(graph_df)
    matched_rows = load_matched_rows(base_run_dir)
    output_runs = load_latest_runs(output_dir / "runs")
    archived_runs = load_latest_runs(base_run_dir / "runs") if not args.no_reuse_archived else {}
    available_runs: Dict[Tuple[str, str], Any] = dict(archived_runs)
    available_runs.update(output_runs)
    source_lookup: Dict[Tuple[str, str], str] = {key: "archived_full_run" for key in archived_runs}
    source_lookup.update({key: "existing_target_run" for key in output_runs})

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    classifier = None
    rows: List[Dict[str, Any]] = []
    run_records: List[PanelRecord] = []

    for drug, disease in target_pairs:
        drug = drug.strip()
        disease = disease.strip()
        canonical_drug, canonical_disease, canonical_resolution = resolve_canonical_pair(
            drug,
            disease,
            matched_rows,
        )
        key = pair_key(canonical_drug, canonical_disease)
        graph_row = graph_rows.get(key, {})
        record = matching_run_record(drug, disease, available_runs, source_lookup)
        if record is not None:
            copied_path = copy_archived_run(record, output_dir / "runs")
            record = PanelRecord(
                drug=drug,
                disease=disease,
                timestamp=record.timestamp,
                path=copied_path,
                components=record.components,
                source=record.source,
            )
        else:
            if classifier is None:
                classifier = make_classifier(params, output_dir)
            record = compute_fresh_record(
                drug=drug,
                disease=disease,
                classifier=classifier,
                params=params,
                graph_row=graph_row,
                output_dir=output_dir,
                timestamp=timestamp,
            )
        run_records.append(record)
        rows.append(
            build_panel_row(
                record,
                graph_row,
                matched_rows,
                graph_df,
                canonical_drug,
                canonical_disease,
                canonical_resolution,
                params,
            )
        )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "target_pair_framework_summary.csv"
    md_path = output_dir / "target_pair_framework_summary.md"
    audit_path = output_dir / "target_pair_run_manifest.json"
    display_columns = [
        "drug_disease_pair",
        "coverage_tier",
        "literature_composition",
        "safety_signal",
        "graph_support",
        "bayesian_posterior",
        "evidence_readiness_score",
        "quality_flag",
        "interpretation",
    ]

    df.to_csv(csv_path, index=False)
    write_markdown_table(df, md_path, display_columns)
    audit = {
        "timestamp": timestamp,
        "base_run_dir": rel_path(base_run_dir),
        "output_dir": rel_path(output_dir),
        "parameters": params,
        "target_pairs": [{"drug": drug, "disease": disease} for drug, disease in target_pairs],
        "run_logs": [
            {
                "drug": record.drug,
                "disease": record.disease,
                "source": record.source,
                "run_log": rel_path(record.path),
                "timestamp": record.timestamp,
            }
            for record in run_records
        ],
    }
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"Wrote CSV: {rel_path(csv_path)}")
    print(f"Wrote Markdown: {rel_path(md_path)}")
    print()
    print(dataframe_to_markdown(df, display_columns))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
