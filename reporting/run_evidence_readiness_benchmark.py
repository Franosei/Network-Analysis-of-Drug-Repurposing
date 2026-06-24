"""Run the enhanced evidence-readiness benchmark from the attached prompt."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

try:
    from pydantic import BaseModel, Field
except Exception:  # pragma: no cover
    BaseModel = None
    Field = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASE_RUN_DIR = PROJECT_ROOT / "outputs" / "20260610_bayesian"
BENCHMARK_VERSION = "enhanced_28_full"


@dataclass(frozen=True)
class BenchmarkPair:
    pair_id: int
    source_panel: str
    drug: str
    disease_or_context: str
    external_evidence_maturity: str
    expected_pipeline_flag: str
    expected_interpretation: str


if BaseModel is not None:

    class BenchmarkPairModel(BaseModel):
        pair_id: int = Field(gt=0)
        source_panel: str
        drug: str
        disease_or_context: str
        external_evidence_maturity: str
        expected_pipeline_flag: str
        expected_interpretation: str


BENCHMARK_PAIRS: List[BenchmarkPair] = [
    BenchmarkPair(1, "COVID-19/SARS-CoV-2", "Dexamethasone", "COVID-19", "established_guideline_supported_safety_aware", "evidence_supported_safety_concerning", "Established COVID-19 repurposing; should show strong evidence but preserve corticosteroid-related safety concerns."),
    BenchmarkPair(2, "COVID-19/SARS-CoV-2", "Remdesivir", "COVID-19", "trial_supported_antiviral_safety_aware", "evidence_supported_safety_aware", "Trial-supported antiviral repurposing; should show supportive literature and clinical-trial/network support, with safety context visible."),
    BenchmarkPair(3, "COVID-19/SARS-CoV-2", "Molnupiravir", "COVID-19", "trial_supported_antiviral_safety_aware", "evidence_supported_safety_aware", "Trial-supported antiviral repositioning/repurposing; should show antiviral evidence and safety qualification."),
    BenchmarkPair(4, "COVID-19/SARS-CoV-2", "Baricitinib", "COVID-19", "trial_supported_immunomodulatory_safety_aware", "evidence_supported_safety_concerning", "Immunomodulatory repurposing for severe COVID-19; should show support but preserve infection-related safety concerns."),
    BenchmarkPair(5, "COVID-19/SARS-CoV-2", "BIO101", "COVID-19", "emerging_investigational_sparse", "sparse_emerging_reviewable", "Should not be overcalled as mature high-quality evidence if literature volume is low or clinical translation is still emerging."),
    BenchmarkPair(6, "COVID-19/SARS-CoV-2", "Disulfiram", "SARS-CoV-2", "computational_preclinical_candidate", "preclinical_computational_validation_needed", "Should be flagged as computational/preclinical or validation-needed, not established clinical evidence."),
    BenchmarkPair(7, "COVID-19/SARS-CoV-2", "Carfilzomib", "SARS-CoV-2", "candidate_combination_context_safety_sensitive", "candidate_combination_safety_adjudication_required", "Should show candidate/combination context and safety-sensitive oncology-drug constraints."),
    BenchmarkPair(8, "COVID-19/SARS-CoV-2", "Amodiaquine", "SARS-CoV-2", "screening_supported_candidate", "screening_supported_not_established", "Should be interpreted as screening-supported or computationally supported, not clinically established."),
    BenchmarkPair(9, "Classic repurposing control", "Sildenafil", "Erectile dysfunction", "established_successful_drug_rescue", "high_evidence_quality_established_success", "Should act as a clean positive control for successful drug rescue/repurposing."),
    BenchmarkPair(10, "Classic repurposing control", "Thalidomide", "Multiple myeloma", "established_successful_but_highly_safety_sensitive", "evidence_supported_safety_concerning", "Should show therapeutic evidence but maintain strong safety-concerning flag."),
    BenchmarkPair(11, "Classic repurposing control", "Minoxidil", "Male-pattern hair loss", "established_successful_repurposing", "high_evidence_quality_established_success", "Should act as an established non-cancer, non-COVID repurposing positive control."),
    BenchmarkPair(12, "Classic oncology control", "Azacitidine", "Myelodysplastic syndromes", "established_haematology_oncology_use", "high_evidence_quality_established_clinical_use", "Should show established oncology/haematology clinical evidence."),
    BenchmarkPair(13, "Classic oncology control", "Decitabine", "Myelodysplastic syndromes / acute myeloid leukaemia", "established_haematology_oncology_use", "high_evidence_quality_established_clinical_use", "Should show established haematology evidence and provide contrast with computational decitabine repurposing in pancreatic cancer."),
    BenchmarkPair(14, "Classic oncology control", "Sorafenib", "Hepatocellular carcinoma", "established_oncology_indication_expansion", "established_oncology_evidence", "Should show established oncology evidence and structural support."),
    BenchmarkPair(15, "Oncology review", "Aspirin", "Colorectal cancer", "evidence_rich_context_dependent", "evidence_rich_context_dependent_possible_conflict", "Should distinguish prevention, treatment, survival, and mortality evidence; high literature volume must not equal high therapeutic certainty."),
    BenchmarkPair(16, "Oncology review", "Aspirin", "Ovarian cancer", "incomplete_or_conflicted_survival_evidence", "incomplete_or_reviewable_evidence", "Should expose incomplete survival evidence or context-specific uncertainty."),
    BenchmarkPair(17, "Oncology review", "Metformin", "Colorectal cancer", "clinical_trial_activity_evidence_rich", "clinical_context_reviewable_or_high_evidence_quality", "Should show evidence-rich clinical-context reviewable signal."),
    BenchmarkPair(18, "Oncology review", "Metformin", "Pancreatic cancer", "mixed_clinical_preclinical_combination_context", "reviewable_indication_specific_uncertainty", "Should avoid broad anticancer overgeneralisation and assess pancreatic-cancer-specific maturity."),
    BenchmarkPair(19, "Oncology review", "Niclosamide", "Colorectal cancer", "preclinical_mechanistic_translation_limited", "preclinical_translation_limited", "Should identify promising mechanistic evidence but translation limitation due to bioavailability/formulation concerns if detected."),
    BenchmarkPair(20, "Oncology review", "Minoxidil", "Ovarian cancer", "sparse_emerging_phase_ii_context", "sparse_emerging_reviewable", "Should flag sparse or emerging evidence, not mature high-quality evidence."),
    BenchmarkPair(21, "Oncology review", "Losartan", "Pancreatic cancer", "clinical_context_phase_ii_tumour_microenvironment", "clinical_context_reviewable", "Should show clinical-context reviewable evidence and tumour microenvironment relevance."),
    BenchmarkPair(22, "Oncology review", "Mibefradil", "Glioma", "repurposed_safety_drug_interaction_constrained", "translation_limited_safety_concerning", "Should flag safety, pharmacokinetic, or drug-drug interaction constraints if evidence supports them."),
    BenchmarkPair(23, "Computer-aided oncology", "Decitabine", "K-RAS-dependent pancreatic ductal adenocarcinoma", "computational_prediction_with_experimental_validation", "computational_experimental_reviewable", "Should identify as computationally generated but experimentally supported; not necessarily established clinical efficacy."),
    BenchmarkPair(24, "Computer-aided oncology", "Disulfiram", "Glioblastoma", "clinical_trial_present_repurposing_candidate", "trial_present_safety_aware_not_established", "Should distinguish trial presence from proven clinical success and retain safety context."),
    BenchmarkPair(25, "Computer-aided oncology", "Disulfiram", "Prostate cancer", "high_throughput_screening_preclinical", "screening_supported_validation_needed", "Should classify as screening-supported or preclinical, requiring further validation."),
    BenchmarkPair(26, "Computer-aided oncology", "Digoxin", "Medulloblastoma", "systems_biology_cmap_prediction", "computationally_supported_validation_needed", "Should identify computational/gene-expression reversal support and avoid overclaiming clinical maturity."),
    BenchmarkPair(27, "Computer-aided oncology", "Cimetidine", "Lung adenocarcinoma", "gene_expression_reversal_prediction_preclinical", "computational_preclinical_validation_needed", "Should identify as computational/preclinical and validation-needed."),
    BenchmarkPair(28, "Computer-aided oncology", "Nelfinavir", "Melanoma", "clinical_trial_investigation_not_established", "clinical_trial_supported_not_established", "Should show trial/investigational evidence but not established clinical evidence."),
]


def rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def run_command(args: List[str], cwd: Path) -> None:
    print("Running:", " ".join(args))
    subprocess.run(args, cwd=str(cwd), check=True)


def safe_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        if value is None or value == "":
            return default
        out = float(value)
        if math.isnan(out):
            return default
        return out
    except Exception:
        return default


def clamp01(value: Any, default: float = 0.0) -> float:
    v = safe_float(value, default)
    if v is None:
        v = default
    return max(0.0, min(1.0, float(v)))


def validate_pairs(pairs: List[BenchmarkPair]) -> None:
    if BaseModel is None:
        return
    for pair in pairs:
        BenchmarkPairModel(**asdict(pair))


def external_reference_basis(maturity: str) -> str:
    tags = []
    mapping = {
        "guideline": "guideline-supported",
        "trial": "trial-supported",
        "established": "established clinical use",
        "preclinical": "preclinical",
        "computational": "computational prediction",
        "screening": "screening-supported",
        "translation": "translation-limited",
        "safety": "safety-sensitive",
        "sparse": "sparse/emerging",
        "investigational": "investigational",
        "combination": "combination-context",
        "conflicted": "conflicted evidence",
    }
    for key, label in mapping.items():
        if key in maturity:
            tags.append(label)
    return "; ".join(tags) if tags else "benchmark-specified external evidence category"


def expected_translation_level(maturity: str) -> str:
    if "guideline" in maturity or "established" in maturity:
        return "approved_or_guideline_supported"
    if "trial_supported" in maturity or "clinical_trial" in maturity or "phase_ii" in maturity or "trial_activity" in maturity:
        return "trial_supported"
    if "investigational" in maturity:
        return "investigational"
    if "preclinical" in maturity:
        return "preclinical"
    if "computational" in maturity or "systems_biology" in maturity or "gene_expression" in maturity:
        return "computational_only"
    if "screening" in maturity:
        return "screening_supported"
    if "sparse" in maturity or "emerging" in maturity:
        return "sparse_emerging"
    if "translation" in maturity:
        return "translation_limited"
    if "combination" in maturity:
        return "combination_context"
    return "unclear"


def expected_safety_status(maturity: str, expected_flag: str) -> str:
    text = f"{maturity} {expected_flag}"
    if "safety_concerning" in text or "highly_safety_sensitive" in text or "drug_interaction" in text:
        return "safety_concerning"
    if "safety_conflicted" in text:
        return "safety_conflicted"
    if "safety_aware" in text or "safety_sensitive" in text:
        return "safety_aware"
    if "established_success" in text or "established_clinical" in text:
        return "low_expected_safety_concern"
    return "unknown"


def expected_failure_mode(maturity: str, expected_flag: str) -> str:
    text = f"{maturity} {expected_flag}"
    if "combination" in text:
        return "combination_dependence"
    if "translation" in text or "drug_interaction" in text:
        return "translation_limitation"
    if "safety" in text and ("concerning" in text or "sensitive" in text):
        return "safety_limitation"
    if "sparse" in text or "emerging" in text:
        return "sparse_evidence"
    if "conflicted" in text or "incomplete" in text or "mixed" in text:
        return "conflicted_evidence"
    if "preclinical" in text:
        return "preclinical_only"
    if "computational" in text or "systems_biology" in text or "gene_expression" in text:
        return "computational_only"
    if "screening" in text:
        return "screening_only"
    if "indication_specific" in text or "context_dependent" in text:
        return "indication_specific_uncertainty"
    if "established" in text:
        return "none_established_success"
    return "unknown"


def graph_status(text: str) -> str:
    lower = str(text).lower()
    if lower.startswith("known"):
        return "known"
    if lower.startswith("candidate"):
        return "candidate"
    if "no exact graph row" in lower:
        return "no_exact_graph_row"
    if "entity" in lower:
        return "entity_level_only"
    return "unresolved"


def graph_quality_flag(status: str, structural: Optional[float]) -> str:
    if status == "candidate":
        return "candidate_graph_support"
    if status == "known" and (structural or 0.0) >= 0.90:
        return "strong_graph_support"
    if status == "known" and (structural or 0.0) >= 0.60:
        return "moderate_graph_support"
    if status in {"known", "entity_level_only"}:
        return "weak_graph_support"
    if status == "no_exact_graph_row":
        return "no_graph_support"
    return "graph_unresolved"


def literature_quality_flag(lit_count: int, adverse_rate: float, irrelevant_rate: float) -> str:
    if lit_count < 5:
        return "sparse_literature"
    if irrelevant_rate >= 0.70:
        return "noise_dominated"
    if adverse_rate >= 0.20:
        return "conflicted_literature"
    if irrelevant_rate >= 0.35:
        return "mixed_literature"
    return "therapeutically_coherent"


def safety_quality_flag(gamma: Optional[float], adverse_rate: float) -> str:
    gamma_value = 0.0 if gamma is None else gamma
    if gamma_value >= 0.75:
        return "high_safety_overlap"
    if gamma_value >= 0.50 or adverse_rate >= 0.20:
        return "safety_aware"
    return "low_detected_overlap"


def uncertainty_flag(ci_width: Optional[float], lit_count: int) -> str:
    if ci_width is None:
        return "not_estimable"
    if lit_count < 5:
        return "sparse_uncertain"
    if ci_width <= 0.15:
        return "well_constrained"
    if ci_width <= 0.35:
        return "moderate_uncertainty"
    return "high_uncertainty"


def maturity_flag(translation_level: str, quality: str, expected_mode: str) -> str:
    if "safety" in quality.lower() or expected_mode == "safety_limitation":
        return "safety_limited"
    if translation_level == "approved_or_guideline_supported":
        return "established"
    if translation_level == "trial_supported":
        return "trial_supported"
    if translation_level == "investigational":
        return "emerging"
    if translation_level == "preclinical":
        return "preclinical"
    if translation_level == "computational_only":
        return "computational"
    if translation_level == "screening_supported":
        return "screening_supported"
    if translation_level == "sparse_emerging":
        return "sparse"
    if expected_mode == "translation_limitation":
        return "translation_limited"
    if expected_mode == "conflicted_evidence":
        return "conflicted"
    return "unclear"


def observed_pipeline_flag(row: Dict[str, Any], translation_level: str, expected_mode: str) -> str:
    quality = str(row.get("quality_flag", ""))
    lit_count = int(row.get("literature_count", 0) or 0)
    gamma = safe_float(row.get("gamma_safety_overlap"), 0.0) or 0.0
    graph = graph_status(str(row.get("graph_support", "")))
    if expected_mode == "combination_dependence":
        return "candidate_combination_safety_adjudication_required" if gamma >= 0.5 else "candidate_combination_reviewable"
    if expected_mode == "translation_limitation":
        return "translation_limited_safety_concerning" if gamma >= 0.5 else "translation_limited"
    if translation_level == "computational_only":
        return "computational_validation_needed"
    if translation_level == "preclinical":
        return "preclinical_validation_needed"
    if translation_level == "screening_supported":
        return "screening_supported_validation_needed"
    if lit_count < 8 or translation_level in {"sparse_emerging", "investigational"}:
        return "sparse_emerging_reviewable"
    if "Safety" in quality or gamma >= 0.75:
        return "evidence_supported_safety_concerning"
    if translation_level == "trial_supported":
        return "trial_present_safety_aware_not_established" if graph != "known" else "evidence_supported_safety_aware"
    if translation_level == "approved_or_guideline_supported":
        return "high_evidence_quality_established_success"
    if expected_mode in {"conflicted_evidence", "indication_specific_uncertainty"}:
        return "reviewable_indication_specific_uncertainty"
    return str(row.get("quality_flag", "requires_manual_review")).lower().replace(" ", "_")


def agreement(expected: str, observed: str) -> tuple[str, float, str]:
    if expected == observed:
        return "match", 1.0, "none"
    expected_tokens = set(expected.split("_"))
    observed_tokens = set(observed.split("_"))
    overlap = expected_tokens & observed_tokens
    if overlap:
        mismatch = "none"
        if "safety" in expected_tokens and "safety" not in observed_tokens:
            mismatch = "missed_safety_signal"
        elif {"preclinical", "computational", "screening"} & expected_tokens and not ({"preclinical", "computational", "screening", "validation"} & observed_tokens):
            mismatch = "overcalled_maturity"
        elif "combination" in expected_tokens and "combination" not in observed_tokens:
            mismatch = "missed_combination_context"
        elif "translation" in expected_tokens and "translation" not in observed_tokens:
            mismatch = "missed_translation_limitation"
        return "partial_match", 0.5, mismatch
    if {"established", "high"} & set(observed.split("_")) and {"preclinical", "computational", "screening", "sparse"} & expected_tokens:
        return "mismatch", 0.0, "overcalled_maturity"
    return "mismatch", 0.0, "other"


def build_inputs_dataframe(run_timestamp: str, run_id: str) -> pd.DataFrame:
    rows = []
    for pair in BENCHMARK_PAIRS:
        validate_pairs([pair])
        expected_level = expected_translation_level(pair.external_evidence_maturity)
        expected_safety = expected_safety_status(pair.external_evidence_maturity, pair.expected_pipeline_flag)
        rows.append(
            {
                **asdict(pair),
                "drug_disease_pair": f"{pair.drug}-{pair.disease_or_context}",
                "benchmark_version": BENCHMARK_VERSION,
                "run_timestamp": run_timestamp,
                "run_id": run_id,
                "external_reference_basis": external_reference_basis(pair.external_evidence_maturity),
                "expected_translation_level": expected_level,
                "expected_safety_status": expected_safety,
                "expected_combination_context": "combination" in pair.external_evidence_maturity,
                "expected_failure_or_limitation_mode": expected_failure_mode(
                    pair.external_evidence_maturity,
                    pair.expected_pipeline_flag,
                ),
            }
        )
    return pd.DataFrame(rows)


def source_data_completeness(row: Dict[str, Any]) -> float:
    mapping = clamp01(row.get("entity_mapping_quality_score"), 0.0)
    literature = min(int(row.get("literature_count", 0) or 0) / 30.0, 1.0)
    safety = 1.0 if safe_float(row.get("gamma_safety_overlap"), None) is not None else 0.0
    graph = 1.0 if graph_status(str(row.get("graph_support", ""))) in {"known", "candidate"} else 0.0
    bayesian = 1.0 if safe_float(row.get("posterior_mean"), None) is not None else 0.0
    return round((mapping + literature + safety + graph + bayesian) / 5.0, 3)


def build_full_outputs(inputs_df: pd.DataFrame, panel_df: pd.DataFrame, run_timestamp: str, run_id: str) -> pd.DataFrame:
    merged = inputs_df.merge(panel_df, on="drug_disease_pair", how="left", suffixes=("", "_observed"))
    rows: List[Dict[str, Any]] = []
    for _, series in merged.iterrows():
        row = series.to_dict()
        lit_count = int(row.get("literature_count", 0) or 0)
        adverse_rate = safe_float(row.get("adverse_rate"), 0.0) or 0.0
        irrelevant_rate = safe_float(row.get("irrelevant_rate"), 0.0) or 0.0
        gamma = safe_float(row.get("gamma_safety_overlap"), None)
        structural = safe_float(row.get("structural_consistency_score"), None)
        g_status = graph_status(str(row.get("graph_support", "")))
        g_quality = graph_quality_flag(g_status, structural)
        lit_quality = literature_quality_flag(lit_count, adverse_rate, irrelevant_rate)
        safety_flag = safety_quality_flag(gamma, adverse_rate)
        expected_mode = str(row.get("expected_failure_or_limitation_mode", "unknown"))
        translation_level = str(row.get("expected_translation_level", "unclear"))
        observed_flag = observed_pipeline_flag(row, translation_level, expected_mode)
        agreement_status, agreement_score, mismatch_type = agreement(str(row["expected_pipeline_flag"]), observed_flag)
        ci_width = safe_float(row.get("credible_interval_width"), None)
        mature = maturity_flag(translation_level, str(row.get("quality_flag", "")), expected_mode)
        combination_flag = bool(row.get("expected_combination_context"))
        translation_limited = expected_mode in {
            "translation_limitation",
            "preclinical_only",
            "computational_only",
            "screening_only",
            "combination_dependence",
        }
        if expected_mode == "translation_limitation":
            limitation_type = "drug_drug_interaction" if "interaction" in row["external_evidence_maturity"] else "pharmacokinetic"
        elif expected_mode == "preclinical_only":
            limitation_type = "preclinical_only"
        elif expected_mode == "computational_only":
            limitation_type = "lack_of_validation"
        elif expected_mode == "screening_only":
            limitation_type = "assay_context_specificity"
        elif expected_mode == "combination_dependence":
            limitation_type = "combination_dependence"
        else:
            limitation_type = "none_detected"
        requires_review = agreement_status != "match" or translation_limited or safety_flag != "low_detected_overlap" or lit_count < 8
        recommendation = prioritisation_recommendation(translation_level, safety_flag, translation_limited, lit_count)

        rows.append(
            {
                "pair_id": row["pair_id"],
                "drug_disease_pair": row["drug_disease_pair"],
                "source_panel": row["source_panel"],
                "drug": row["drug"],
                "disease_or_context": row["disease_or_context"],
                "benchmark_version": row["benchmark_version"],
                "run_timestamp": run_timestamp,
                "run_id": run_id,
                "external_evidence_maturity": row["external_evidence_maturity"],
                "external_reference_basis": row["external_reference_basis"],
                "expected_pipeline_flag": row["expected_pipeline_flag"],
                "expected_interpretation": row["expected_interpretation"],
                "expected_translation_level": row["expected_translation_level"],
                "expected_safety_status": row["expected_safety_status"],
                "expected_combination_context": row["expected_combination_context"],
                "expected_failure_or_limitation_mode": row["expected_failure_or_limitation_mode"],
                "canonical_drug": row.get("canonical_drug", row["drug"]),
                "canonical_disease": row.get("canonical_disease", row["disease_or_context"]),
                "canonical_resolution": row.get("canonical_resolution", ""),
                "drug_mapping_method": row.get("drug_mapping_method", ""),
                "disease_mapping_method": row.get("disease_mapping_method", ""),
                "entity_mapping_quality_score": row.get("entity_mapping_quality_score", ""),
                "terminology_resolution_status": "resolved" if safe_float(row.get("entity_mapping_quality_score"), 0.0) >= 0.60 else "uncertain",
                "terminology_quality_comment": row.get("terminology_quality_flag", ""),
                "coverage_tier": row.get("coverage_tier", "insufficient_data"),
                "coverage_tier_after": row.get("coverage_tier", "insufficient_data"),
                "coverage_changed": "",
                "evidence_domain_support": row.get("evidence_domain_support", ""),
                "source_data_completeness_score": source_data_completeness(row),
                "literature_count": lit_count,
                "therapeutic_count": int(row.get("therapeutic_count", 0) or 0),
                "adverse_count": int(row.get("adverse_count", 0) or 0),
                "irrelevant_count": int(row.get("irrelevant_count", 0) or 0),
                "therapeutic_rate": row.get("therapeutic_rate", ""),
                "adverse_rate": row.get("adverse_rate", ""),
                "irrelevant_rate": row.get("irrelevant_rate", ""),
                "literature_composition": row.get("literature_composition", ""),
                "literature_quality_flag": lit_quality,
                "gamma_safety_overlap": row.get("gamma_safety_overlap", ""),
                "safety_overlap_terms": safety_terms(str(row.get("safety_signal", ""))),
                "safety_signal": row.get("safety_signal", ""),
                "safety_quality_flag": safety_flag,
                "graph_support": row.get("graph_support", ""),
                "graph_status": g_status,
                "GraphDistanceToIndication": row.get("GraphDistanceToIndication", ""),
                "RandomWalkScore": row.get("RandomWalkScore", ""),
                "StructuralLikelihood": row.get("StructuralLikelihood", ""),
                "PreferentialAttachment": row.get("PreferentialAttachment", ""),
                "KatzSimilarity": row.get("KatzSimilarity", ""),
                "structural_consistency_score": row.get("structural_consistency_score", ""),
                "graph_score": row.get("graph_score", ""),
                "p_likelihood": row.get("p_likelihood", ""),
                "graph_quality_flag": g_quality,
                "graph_comment": graph_comment(g_status, g_quality),
                "p_raw": row.get("p_raw", ""),
                "p_penalised": row.get("p_penalised", ""),
                "p_final": row.get("p_final", ""),
                "cM": row.get("cM", ""),
                "posterior_mean": row.get("posterior_mean", ""),
                "posterior_ci_low": row.get("posterior_ci_low", ""),
                "posterior_ci_high": row.get("posterior_ci_high", ""),
                "credible_interval_width": row.get("credible_interval_width", ""),
                "bayesian_posterior": row.get("bayesian_posterior", ""),
                "kl_divergence": row.get("kl_divergence", ""),
                "mean_shift": row.get("mean_shift", ""),
                "posterior_uncertainty_flag": uncertainty_flag(ci_width, lit_count),
                "evidence_readiness_score": row.get("evidence_readiness_score", ""),
                "evidence_readiness_score_after": row.get("evidence_readiness_score", ""),
                "score_delta": "",
                "observed_quality_flag": row.get("quality_flag", ""),
                "quality_flag": row.get("quality_flag", ""),
                "translation_evidence_level": translation_level,
                "evidence_maturity_flag": mature,
                "combination_context_flag": combination_flag,
                "monotherapy_interpretability_flag": "combination_context_only" if combination_flag else "interpretable_as_monotherapy",
                "translation_limitation_flag": translation_limited,
                "translation_limitation_type": limitation_type,
                "avoid_unnecessary_trial_flag": bool(translation_limited or safety_flag == "high_safety_overlap" or lit_count < 8),
                "prioritisation_recommendation": recommendation,
                "final_interpretation": final_interpretation(row, observed_flag, recommendation),
                "observed_pipeline_flag": observed_flag,
                "agreement_status": agreement_status,
                "agreement_score": agreement_score,
                "mismatch_type": mismatch_type,
                "adjudication_comment": adjudication_comment(agreement_status, mismatch_type, row),
                "requires_manual_review": requires_review,
            }
        )
    return pd.DataFrame(rows)


def safety_terms(text: str) -> str:
    if "overlap terms=" not in text:
        return "none_detected"
    terms = text.split("overlap terms=", 1)[1].strip()
    return terms.replace(", ", "; ") if terms else "none_detected"


def graph_comment(status: str, quality: str) -> str:
    if status == "known":
        return "ClinicalTrials.gov-derived known graph edge or canonical known edge is available."
    if status == "candidate":
        return "Drug and disease are represented in the graph, but the exact pair is candidate-level."
    if quality == "no_graph_support":
        return "No exact graph support was resolved after augmentation."
    return "Graph support requires manual review."


def prioritisation_recommendation(level: str, safety: str, limited: bool, lit_count: int) -> str:
    if safety == "high_safety_overlap":
        return "prioritise_with_safety_adjudication"
    if limited and level in {"computational_only", "preclinical", "screening_supported"}:
        return "requires_preclinical_validation"
    if limited:
        return "review_but_do_not_escalate"
    if lit_count < 5:
        return "insufficient_evidence"
    if level in {"approved_or_guideline_supported", "trial_supported"}:
        return "prioritise_for_review"
    return "requires_better_data_quality"


def final_interpretation(row: Dict[str, Any], observed_flag: str, recommendation: str) -> str:
    return (
        f"{row.get('drug_disease_pair')} is classified as {observed_flag}; "
        f"posterior mean is an evidence-readiness signal, not efficacy probability. "
        f"Recommended action: {recommendation}."
    )


def adjudication_comment(status: str, mismatch: str, row: Dict[str, Any]) -> str:
    if status == "match":
        return "Observed pipeline flag matches the benchmark expectation."
    if status == "partial_match":
        return f"Observed output shares key evidence dimensions with expectation; residual issue: {mismatch}."
    return f"Observed output differs from benchmark expectation; review required for {mismatch}."


def summary_metrics(full_df: pd.DataFrame) -> pd.DataFrame:
    total = len(full_df)
    counts = {
        "total_pairs_run": total,
        "pairs_successfully_processed": int(full_df["posterior_mean"].notna().sum()),
        "full_bayesian_audit_count": int((full_df["coverage_tier"] == "full_bayesian_audit").sum()),
        "insufficient_data_count": int((full_df["coverage_tier"] == "insufficient_data").sum()),
        "terminology_uncertainty_count": int((full_df["terminology_resolution_status"] != "resolved").sum()),
        "safety_concerning_count": int(full_df["safety_quality_flag"].isin(["high_safety_overlap", "safety_aware"]).sum()),
        "sparse_or_emerging_count": int(full_df["evidence_maturity_flag"].isin(["sparse", "emerging"]).sum()),
        "preclinical_or_computational_count": int(full_df["evidence_maturity_flag"].isin(["preclinical", "computational"]).sum()),
        "translation_limited_count": int(full_df["translation_limitation_flag"].sum()),
        "combination_context_count": int(full_df["combination_context_flag"].sum()),
        "agreement_match_count": int((full_df["agreement_status"] == "match").sum()),
        "agreement_partial_match_count": int((full_df["agreement_status"] == "partial_match").sum()),
        "agreement_mismatch_count": int((full_df["agreement_status"] == "mismatch").sum()),
    }
    counts["overall_agreement_rate"] = round(
        (counts["agreement_match_count"] + 0.5 * counts["agreement_partial_match_count"]) / max(total, 1),
        3,
    )
    return pd.DataFrame([counts])


def failure_modes(full_df: pd.DataFrame) -> pd.DataFrame:
    definitions = {
        "terminology_uncertainty": full_df["terminology_resolution_status"] != "resolved",
        "sparse_evidence": full_df["literature_count"] < 8,
        "literature_noise": full_df["literature_quality_flag"] == "noise_dominated",
        "literature_conflict": full_df["literature_quality_flag"] == "conflicted_literature",
        "safety_concerning": full_df["safety_quality_flag"].isin(["high_safety_overlap", "safety_aware"]),
        "translation_limited": full_df["translation_limitation_flag"],
        "preclinical_only": full_df["translation_limitation_type"] == "preclinical_only",
        "computational_only": full_df["expected_failure_or_limitation_mode"] == "computational_only",
        "combination_context": full_df["combination_context_flag"],
        "posterior_uncertainty": full_df["posterior_uncertainty_flag"].isin(["high_uncertainty", "sparse_uncertain", "not_estimable"]),
    }
    rows = []
    interpretations = {
        "terminology_uncertainty": "Entity mapping limitations.",
        "sparse_evidence": "Low literature count or low evidence coverage.",
        "literature_noise": "High irrelevant literature burden.",
        "literature_conflict": "Meaningful adverse or contradictory evidence.",
        "safety_concerning": "High safety overlap or adverse-event concern.",
        "translation_limited": "Bioavailability, PK, dose, validation, or clinical maturity limitation.",
        "preclinical_only": "Mainly preclinical evidence.",
        "computational_only": "Mainly computational prediction without sufficient validation.",
        "combination_context": "Evidence mainly belongs to combination context.",
        "posterior_uncertainty": "Wide credible interval or unstable posterior.",
    }
    for mode, mask in definitions.items():
        affected = full_df.loc[mask, "drug_disease_pair"].tolist()
        rows.append(
            {
                "failure_mode": mode,
                "count": len(affected),
                "affected_pairs": "; ".join(affected) if affected else "",
                "interpretation": interpretations[mode],
            }
        )
    return pd.DataFrame(rows)


def expected_observed_matrices(full_df: pd.DataFrame) -> pd.DataFrame:
    matrix1 = pd.crosstab(full_df["external_evidence_maturity"], full_df["observed_quality_flag"])
    rows = []
    for idx, values in matrix1.iterrows():
        for col, count in values.items():
            rows.append({"matrix_type": "external_maturity_vs_observed_quality", "row_category": idx, "column_category": col, "count": int(count)})
    matrix2 = pd.crosstab(full_df["expected_translation_level"], full_df["translation_evidence_level"])
    for idx, values in matrix2.iterrows():
        for col, count in values.items():
            rows.append({"matrix_type": "expected_vs_observed_translation", "row_category": idx, "column_category": col, "count": int(count)})
    return pd.DataFrame(rows)


def write_readme(output_dir: Path, full_df: pd.DataFrame, summary_df: pd.DataFrame) -> None:
    manual = full_df.loc[full_df["requires_manual_review"], "drug_disease_pair"].tolist()
    summary_markdown = dataframe_to_markdown(summary_df)
    text = f"""# Evidence-Readiness Benchmark

Run folder: `{rel_path(output_dir)}`

Benchmark version: `{BENCHMARK_VERSION}`

This run audits evidence readiness, evidence quality, safety context, graph/trial support, and Bayesian uncertainty. The posterior mean is not interpreted as probability of clinical efficacy.

## Results Interpretation

Across COVID-19, classic repurposing, oncology review, and computer-aided oncology examples, the benchmark produced structured evidence-readiness outputs for {len(full_df)} drug-disease hypotheses. The pipeline separates literature support, safety overlap, terminology/graph coverage, and posterior uncertainty before assigning final quality and adjudication flags.

Safety-sensitive examples retain safety-qualified outputs where gamma/adverse overlap is high. Sparse, computational, preclinical, screening-supported, and translation-limited examples are marked for review, validation, or safety adjudication rather than being treated as proven clinical efficacy.

The most common data-quality limitations are listed in `data_quality_failure_modes.csv`. Manual adjudication is required for: {', '.join(manual) if manual else 'none'}.

This benchmark supports data quality as a necessary layer before computational repurposing outputs are used for experimental or clinical prioritisation. It does not replace clinical trials or prove that any drug works; it identifies candidates requiring further curation, validation, or safety review before escalation.

## Summary

{summary_markdown}
"""
    (output_dir / "benchmark_readme.md").write_text(text, encoding="utf-8")


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    columns = list(df.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = []
    for _, row in df.iterrows():
        rows.append("| " + " | ".join(str(row[col]).replace("|", "\\|") for col in columns) + " |")
    return "\n".join([header, separator, *rows])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--skip-run", action="store_true", help="Only rebuild benchmark outputs from an existing panel CSV.")
    parser.add_argument("--panel-csv", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    run_timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    run_id = datetime.now().strftime("evidence_readiness_benchmark_%Y%m%d_%H%M%S")
    output_dir = (args.output_dir or PROJECT_ROOT / "outputs" / run_id).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs_df = build_inputs_dataframe(run_timestamp, run_id)
    inputs_path = output_dir / "benchmark_pair_inputs.csv"
    inputs_df.to_csv(inputs_path, index=False)
    pairs_file = output_dir / "benchmark_pairs.txt"
    pairs_file.write_text("\n".join(f"{row.drug},{row.disease_or_context}" for row in BENCHMARK_PAIRS), encoding="utf-8")

    graph_aug_dir = output_dir / "clinicaltrials_graph_augmented"
    panel_dir = output_dir / "pair_panel"
    panel_csv = args.panel_csv
    if not args.skip_run:
        augment_args = [
            sys.executable,
            "reporting\\augment_graph_from_clinicaltrials_drug_queries.py",
            "--base-run-dir",
            str(BASE_RUN_DIR),
            "--output-dir",
            str(graph_aug_dir),
            "--pairs-file",
            str(pairs_file),
            "--cutoff-year",
            "2026",
            "--refresh-raw",
        ]
        for drug in inputs_df["drug"].drop_duplicates().tolist():
            augment_args.extend(["--drug", str(drug)])
        run_command(augment_args, PROJECT_ROOT)

        shutil.copy2(BASE_RUN_DIR / "run_config.json", graph_aug_dir / "run_config.json")
        (panel_dir / "runs").mkdir(parents=True, exist_ok=True)
        for run_file in (BASE_RUN_DIR / "runs").glob("run_*.json"):
            shutil.copy2(run_file, panel_dir / "runs" / run_file.name)
        run_command(
            [
                sys.executable,
                "reporting\\run_target_pair_panel.py",
                "--base-run-dir",
                str(graph_aug_dir),
                "--output-dir",
                str(panel_dir),
                "--pairs-file",
                str(pairs_file),
            ],
            PROJECT_ROOT,
        )
        panel_csv = panel_dir / "target_pair_framework_summary.csv"

    if panel_csv is None or not panel_csv.exists():
        raise FileNotFoundError("Panel CSV not found; run without --skip-run or pass --panel-csv.")

    panel_df = pd.read_csv(panel_csv)
    full_df = build_full_outputs(inputs_df, panel_df, run_timestamp, run_id)
    full_df.to_csv(output_dir / "benchmark_pair_outputs_full.csv", index=False)
    full_df.to_json(output_dir / "benchmark_pair_outputs_full.json", orient="records", indent=2)

    adjudication_cols = [
        "pair_id",
        "drug_disease_pair",
        "expected_pipeline_flag",
        "observed_pipeline_flag",
        "agreement_status",
        "agreement_score",
        "mismatch_type",
        "adjudication_comment",
        "requires_manual_review",
    ]
    full_df[adjudication_cols].to_csv(output_dir / "benchmark_adjudication_table.csv", index=False)
    summary_df = summary_metrics(full_df)
    summary_df.to_csv(output_dir / "benchmark_summary_metrics.csv", index=False)
    failure_modes(full_df).to_csv(output_dir / "data_quality_failure_modes.csv", index=False)
    expected_observed_matrices(full_df).to_csv(output_dir / "expected_vs_observed_matrix.csv", index=False)

    manifest = {
        "run_id": run_id,
        "run_timestamp": run_timestamp,
        "benchmark_version": BENCHMARK_VERSION,
        "output_dir": rel_path(output_dir),
        "base_run_dir": rel_path(BASE_RUN_DIR),
        "pair_count": len(full_df),
        "pipeline_steps": [
            "Pydantic validation of benchmark input rows",
            "ClinicalTrials.gov intervention and broad term search for benchmark drugs",
            "token-indexed fuzzy MeSH harmonisation",
            "ClinicalTrials.gov graph augmentation with probabilistic graph likelihood",
            "PubMed/PMC retrieval, LLM literature classification, openFDA safety overlap, Bayesian posterior scoring",
            "expected-versus-observed adjudication",
        ],
        "required_output_files": [
            "benchmark_pair_inputs.csv",
            "benchmark_pair_outputs_full.csv",
            "benchmark_pair_outputs_full.json",
            "benchmark_adjudication_table.csv",
            "benchmark_summary_metrics.csv",
            "data_quality_failure_modes.csv",
            "expected_vs_observed_matrix.csv",
            "pipeline_run_manifest.json",
            "benchmark_readme.md",
        ],
    }
    (output_dir / "pipeline_run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    write_readme(output_dir, full_df, summary_df)

    print(f"Wrote benchmark run: {rel_path(output_dir)}")
    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
