"""Augment the trial graph with ClinicalTrials.gov records fetched by drug name.

The archived manuscript run is condition/therapeutic-area driven. This helper
adds a complementary intervention-driven pull for target drugs, maps the new
trial conditions/interventions with the existing MeSH mapper, appends the new
known edges to the archived processed pairs, and rebuilds graph features in a
separate output directory.
"""

from __future__ import annotations

import argparse
import json
import random
import re
import shutil
import sys
import time
from collections import Counter
from datetime import datetime
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import networkx as nx

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = PROJECT_ROOT / "code"
for path in (PROJECT_ROOT, CODE_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from condition_drug_pairs import ConditionDrugPairBuilder  # noqa: E402
from data_extraction import ClinicalTrialFetcher  # noqa: E402


DEFAULT_DRUGS: Tuple[str, ...] = (
    "Remdesivir",
    "Molnupiravir",
    "Baricitinib",
    "BIO101",
    "Valproic acid",
    "Disulfiram",
    "Carfilzomib",
    "Amodiaquine",
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


def rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def safe_filename(text: str) -> str:
    return re.sub(r"[^\w\-\.]+", "_", text.strip().lower()).strip("_")


def norm_entity(text: Any) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip().lower())


def flat_key(text: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm_entity(text))


def read_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def unique_by_nct_id(trials: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    out: List[Dict[str, Any]] = []
    for trial in trials:
        nct_id = str(trial.get("nctId", "")).strip()
        if not nct_id or nct_id in seen:
            continue
        seen.add(nct_id)
        out.append(trial)
    return out


def fetch_drug_query_trials(
    fetcher: ClinicalTrialFetcher,
    drug: str,
    *,
    query_parameter: str,
    trials_per_drug: Optional[int],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Fetch eligible trials using a ClinicalTrials.gov query parameter."""
    trials: List[Dict[str, Any]] = []
    next_page_token: Optional[str] = None
    audit: Dict[str, Any] = Counter()
    audit.update(
        {
            "query_drug": drug,
            "query_parameter": query_parameter,
            "raw_returned": 0,
            "eligible_status_phase_drug": 0,
            "excluded_after_cutoff": 0,
            "excluded_missing_nct_id": 0,
            "excluded_duplicate_nct_id": 0,
            "eligible_saved": 0,
            "stopped_at_trials_per_drug_cap": 0,
        }
    )

    while True:
        params = {
            query_parameter: drug,
            "pageSize": fetcher.page_size,
        }
        if next_page_token:
            params["pageToken"] = next_page_token

        result = fetcher._get_with_retries(fetcher.base_url, params=params)
        if result is None:
            audit["request_failed"] = 1
            break

        studies = result.get("studies", []) or []
        audit["raw_returned"] += len(studies)
        next_page_token = result.get("nextPageToken")

        for study in studies:
            if not fetcher.is_eligible_trial(study):
                continue
            audit["eligible_status_phase_drug"] += 1

            if not fetcher.is_within_cutoff(study):
                audit["excluded_after_cutoff"] += 1
                continue

            nct_id = (
                study.get("protocolSection", {})
                .get("identificationModule", {})
                .get("nctId", "")
            )
            if not nct_id:
                audit["excluded_missing_nct_id"] += 1
                continue
            if nct_id in fetcher.seen_nct_ids:
                audit["excluded_duplicate_nct_id"] += 1
                continue

            extracted = fetcher.extract_trial_data(study)
            extracted["sourceTherapeuticArea"] = f"drug_query:{drug}"
            extracted["sourceQueryDrug"] = drug
            extracted["sourceQueryParameter"] = query_parameter
            trials.append(extracted)
            fetcher.seen_nct_ids.add(nct_id)
            audit["eligible_saved"] += 1

            if trials_per_drug is not None and len(trials) >= int(trials_per_drug):
                audit["stopped_at_trials_per_drug_cap"] = 1
                return trials, dict(audit)

        time.sleep(random.uniform(0.2, 0.6))
        if not next_page_token:
            break

    return trials, dict(audit)


def fetch_all_drug_queries(
    drugs: List[str],
    raw_dir: Path,
    args: argparse.Namespace,
) -> pd.DataFrame:
    raw_dir.mkdir(parents=True, exist_ok=True)
    fetcher = ClinicalTrialFetcher(
        output_dir=str(raw_dir),
        page_size=args.page_size,
        cutoff_year=args.cutoff_year,
        include_missing_dates=args.include_missing_dates,
    )

    audit_rows: List[Dict[str, Any]] = []
    for drug in drugs:
        output_path = raw_dir / f"{safe_filename(drug)}_drug_query_upto_{args.cutoff_year}.json"
        if output_path.exists() and not args.refresh_raw:
            trials = read_json(output_path, [])
            trials = trials if isinstance(trials, list) else []
            audit_rows.append(
                {
                    "query_drug": drug,
                    "query_parameter": "reused_existing_raw_file",
                    "raw_returned": "",
                    "eligible_status_phase_drug": "",
                    "excluded_after_cutoff": "",
                    "excluded_missing_nct_id": "",
                    "excluded_duplicate_nct_id": "",
                    "eligible_saved": len(trials),
                    "stopped_at_trials_per_drug_cap": "",
                }
            )
            print(f"[Reuse] {len(trials)} existing trials -> {rel_path(output_path)}")
            continue

        print(f"[ClinicalTrials.gov] intervention query: {drug}")
        all_trials, audit = fetch_drug_query_trials(
            fetcher,
            drug,
            query_parameter=args.query_parameter,
            trials_per_drug=args.trials_per_drug,
        )
        audit_rows.append(audit)
        if args.fallback_query_term and args.query_parameter != "query.term":
            print(f"[ClinicalTrials.gov] broad query.term recall pass: {drug}")
            more_trials, fallback_audit = fetch_drug_query_trials(
                fetcher,
                drug,
                query_parameter="query.term",
                trials_per_drug=args.trials_per_drug,
            )
            fallback_audit["fallback_from"] = args.query_parameter
            audit_rows.append(fallback_audit)
            all_trials.extend(more_trials)

        trials = unique_by_nct_id(all_trials)
        write_json(output_path, trials)
        print(f"[Saved] {len(trials)} eligible trials -> {rel_path(output_path)}")

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(raw_dir.parent / "clinicaltrials_drug_query_audit.csv", index=False)
    return audit_df


def tokenise(text: str) -> List[str]:
    return [token for token in text.split() if token]


def jaccard(left: List[str], right: List[str]) -> float:
    left_set = set(left)
    right_set = set(right)
    if not left_set or not right_set:
        return 0.0
    return len(left_set & right_set) / len(left_set | right_set)


def fast_match_term(
    builder: ConditionDrugPairBuilder,
    raw_term: str,
    mesh_map: Dict[str, str],
    token_index: Dict[str, List[str]],
) -> Tuple[Optional[str], Dict[str, Any]]:
    """Match exact first, then token-indexed candidates only.

    This avoids difflib scanning every MeSH key for each fetched trial term.
    """
    normalised = builder.normalize_for_match(raw_term)
    if not normalised:
        return None, {"method": "empty", "normalized": normalised}
    if normalised in mesh_map:
        return mesh_map[normalised], {
            "method": "exact",
            "normalized": normalised,
            "matched_key": normalised,
            "score": 1.0,
            "sequence_ratio": 1.0,
            "token_jaccard_score": 1.0,
        }

    tokens = tokenise(normalised)
    if not tokens:
        return None, {"method": "no_tokens", "normalized": normalised}

    candidates: List[str] = []
    seen = set()
    for token in tokens:
        for candidate in token_index.get(token, []):
            if candidate not in seen:
                candidates.append(candidate)
                seen.add(candidate)

    if not candidates:
        return None, {"method": "no_candidates", "normalized": normalised}

    best_key = None
    best_score = 0.0
    best_ratio = 0.0
    best_jaccard = 0.0
    for candidate in candidates:
        candidate_tokens = tokenise(candidate)
        token_jaccard = jaccard(tokens, candidate_tokens)
        if token_jaccard < builder.token_jaccard_min:
            continue
        ratio = SequenceMatcher(None, normalised, candidate).ratio()
        score = 0.65 * token_jaccard + 0.35 * ratio
        if score > best_score:
            best_key = candidate
            best_score = score
            best_ratio = ratio
            best_jaccard = token_jaccard

    if best_key is None:
        return None, {"method": "token_score_failed", "normalized": normalised}

    return mesh_map[best_key], {
        "method": "token_score",
        "normalized": normalised,
        "matched_key": best_key,
        "score": round(best_score, 4),
        "sequence_ratio": round(best_ratio, 4),
        "token_jaccard_score": round(best_jaccard, 4),
    }


def split_trial_items(builder: ConditionDrugPairBuilder, values: Any) -> Tuple[List[str], int]:
    if not isinstance(values, list):
        return [], 0
    out: List[str] = []
    placebo_count = 0
    for value in values:
        for part in re.split(r"\s*(?:,|;|\+|/|\band\b)\s*", str(value or "")):
            part = part.strip()
            if not part:
                continue
            if not builder.include_placebo and part.lower() == "placebo":
                placebo_count += 1
                continue
            out.append(part)
    return out, placebo_count


def map_drug_query_trials_fast(
    raw_dir: Path,
    processed_dir: Path,
    mesh_path: Path,
    args: argparse.Namespace,
) -> Tuple[pd.DataFrame, Path, Path]:
    processed_dir.mkdir(parents=True, exist_ok=True)
    matched_path = processed_dir / "condition_drug_pairs.json"
    unmatched_path = processed_dir / "unmatched_pairs.json"
    builder = ConditionDrugPairBuilder(
        input_dir=rel_path(raw_dir),
        output_dir=rel_path(processed_dir),
        output_filename=matched_path.name,
        unmatched_filename=unmatched_path.name,
        mesh_path=rel_path(mesh_path),
        fuzzy_cutoff=args.fuzzy_cutoff,
        token_jaccard_min=args.token_jaccard_min,
        include_placebo=False,
        keep_unmatched_debug_fields=True,
    )

    rows: List[Dict[str, Any]] = []
    unmatched: List[Dict[str, Any]] = []
    audit_rows: List[Dict[str, Any]] = []
    for path in sorted(raw_dir.glob("*.json")):
        trials = read_json(path, [])
        if not isinstance(trials, list):
            continue
        for trial in trials:
            conditions, _ = split_trial_items(builder, trial.get("conditions", []))
            interventions, placebo_count = split_trial_items(builder, trial.get("interventions", []))
            phases = trial.get("phases", trial.get("phase", []))
            if isinstance(phases, str):
                phases = [phases]
            if not isinstance(phases, list):
                phases = []

            base = {
                "nct_id": trial.get("nctId", ""),
                "source_file": path.name,
                "source_therapeutic_area": trial.get("sourceTherapeuticArea", ""),
                "source_query_drug": trial.get("sourceQueryDrug", ""),
                "phases": "|".join(map(str, phases)),
                "status": trial.get("status", ""),
                "placebo_terms_excluded": placebo_count,
            }
            if not conditions:
                audit_rows.append({**base, "mapping_status": "missing_condition"})
                continue
            if not interventions:
                audit_rows.append({**base, "mapping_status": "missing_or_placebo_only_intervention"})
                continue

            for raw_condition in conditions:
                condition_match, condition_debug = fast_match_term(
                    builder,
                    raw_condition,
                    builder.condition_term_map,
                    builder._cond_token_index,
                )
                condition_debug = condition_debug or {}
                mapped_condition = builder.clean_output_name(condition_match) if condition_match else ""
                if not condition_match:
                    unmatched.append(
                        {
                            "raw_condition": raw_condition,
                            "raw_intervention": None,
                            "reason": "unmatched condition",
                            "debug": condition_debug,
                            "nct_id": trial.get("nctId", ""),
                            "source_file": path.name,
                        }
                    )
                    audit_rows.append(
                        {
                            **base,
                            "raw_disease": raw_condition,
                            "mapped_disease": "",
                            "raw_drug": "",
                            "mapped_drug": "",
                            "mapping_status": "unmapped_condition",
                            "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                        }
                    )
                    continue

                for raw_drug in interventions:
                    drug_match, drug_debug = fast_match_term(
                        builder,
                        raw_drug,
                        builder.drug_term_map,
                        builder._drug_token_index,
                    )
                    drug_debug = drug_debug or {}
                    mapped_drug = builder.clean_output_name(drug_match) if drug_match else ""
                    if not drug_match:
                        unmatched.append(
                            {
                                "raw_condition": raw_condition,
                                "mapped_condition": mapped_condition,
                                "raw_intervention": raw_drug,
                                "reason": "unmatched intervention",
                                "debug": drug_debug,
                                "nct_id": trial.get("nctId", ""),
                                "source_file": path.name,
                            }
                        )
                        audit_rows.append(
                            {
                                **base,
                                "raw_disease": raw_condition,
                                "mapped_disease": mapped_condition,
                                "raw_drug": raw_drug,
                                "mapped_drug": "",
                                "mapping_status": "unmapped_intervention",
                                "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                                "drug_debug": json.dumps(drug_debug, ensure_ascii=False),
                            }
                        )
                        continue

                    matched_row = {
                        "condition": mapped_condition,
                        "intervention": mapped_drug,
                        "phases": phases,
                        "status": str(trial.get("status", "")).strip(),
                        "raw_condition": raw_condition,
                        "raw_intervention": raw_drug,
                        "condition_cleaned_term": builder.normalize_for_match(raw_condition),
                        "condition_mapping_method": str(condition_debug.get("method", "")),
                        "condition_mapping_confidence": float(condition_debug.get("score", 0.0) or 0.0),
                        "condition_mapping_score": float(condition_debug.get("score", 0.0) or 0.0),
                        "condition_token_jaccard_score": float(
                            condition_debug.get("token_jaccard_score", 0.0) or 0.0
                        ),
                        "condition_mesh_id": "",
                        "intervention_cleaned_term": builder.normalize_for_match(raw_drug),
                        "intervention_mapping_method": str(drug_debug.get("method", "")),
                        "intervention_mapping_confidence": float(drug_debug.get("score", 0.0) or 0.0),
                        "intervention_mapping_score": float(drug_debug.get("score", 0.0) or 0.0),
                        "intervention_token_jaccard_score": float(
                            drug_debug.get("token_jaccard_score", 0.0) or 0.0
                        ),
                        "intervention_mesh_id": "",
                        "nct_id": trial.get("nctId", ""),
                        "source_file": path.name,
                        "source_therapeutic_area": trial.get("sourceTherapeuticArea", ""),
                        "source_query_drug": trial.get("sourceQueryDrug", ""),
                    }
                    rows.append(matched_row)
                    audit_rows.append(
                        {
                            **base,
                            "raw_disease": raw_condition,
                            "mapped_disease": mapped_condition,
                            "raw_drug": raw_drug,
                            "mapped_drug": mapped_drug,
                            "mapping_status": "matched",
                            "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                            "drug_debug": json.dumps(drug_debug, ensure_ascii=False),
                        }
                    )

    write_json(matched_path, rows)
    write_json(unmatched_path, unmatched)
    return pd.DataFrame(audit_rows), matched_path, unmatched_path


def add_raw_query_drug_fallback_edges(
    raw_dir: Path,
    new_processed_dir: Path,
    mesh_path: Path,
    query_drugs: List[str],
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Add traceable fallback edges when a queried investigational drug is not in MeSH."""
    builder = ConditionDrugPairBuilder(
        input_dir=rel_path(raw_dir),
        output_dir=rel_path(new_processed_dir),
        mesh_path=rel_path(mesh_path),
        fuzzy_cutoff=args.fuzzy_cutoff,
        token_jaccard_min=args.token_jaccard_min,
        include_placebo=False,
        keep_unmatched_debug_fields=True,
    )
    query_flats = {flat_key(drug): drug for drug in query_drugs}
    rows: List[Dict[str, Any]] = []

    for path in sorted(raw_dir.glob("*.json")):
        trials = read_json(path, [])
        if not isinstance(trials, list):
            continue
        for trial in trials:
            source_query = str(trial.get("sourceQueryDrug", "")).strip()
            source_query_flat = flat_key(source_query)
            conditions = builder.normalize_items(trial.get("conditions", []))
            interventions = builder.normalize_items(trial.get("interventions", []))
            phases = trial.get("phases", trial.get("phase", []))
            if isinstance(phases, str):
                phases = [phases]
            if not isinstance(phases, list):
                phases = []

            for raw_condition in conditions:
                condition_match, condition_debug = fast_match_term(
                    builder,
                    raw_condition,
                    builder.condition_term_map,
                    builder._cond_token_index,
                )
                if not condition_match:
                    continue
                mapped_condition = builder.clean_output_name(condition_match)
                condition_debug = condition_debug or {}

                for raw_drug in interventions:
                    raw_flat = flat_key(raw_drug)
                    matched_query = None
                    for query_flat, query_label in query_flats.items():
                        if query_flat and (query_flat in raw_flat or raw_flat in query_flat):
                            matched_query = query_label
                            break
                    if matched_query is None and source_query_flat and source_query_flat in raw_flat:
                        matched_query = source_query
                    if matched_query is None:
                        continue

                    drug_match, _ = fast_match_term(
                        builder,
                        raw_drug,
                        builder.drug_term_map,
                        builder._drug_token_index,
                    )
                    if drug_match:
                        continue

                    fallback_name = matched_query.strip() or raw_drug.strip()
                    rows.append(
                        {
                            "condition": mapped_condition,
                            "intervention": fallback_name,
                            "phases": phases,
                            "status": str(trial.get("status", "")).strip(),
                            "raw_condition": raw_condition,
                            "raw_intervention": raw_drug,
                            "condition_cleaned_term": builder.normalize_for_match(raw_condition),
                            "condition_mapping_method": str(condition_debug.get("method", "")),
                            "condition_mapping_confidence": float(condition_debug.get("score", 0.0) or 0.0),
                            "condition_mapping_score": float(condition_debug.get("score", 0.0) or 0.0),
                            "condition_token_jaccard_score": float(
                                condition_debug.get("token_jaccard_score", 0.0) or 0.0
                            ),
                            "condition_mesh_id": "",
                            "intervention_cleaned_term": builder.normalize_for_match(raw_drug),
                            "intervention_mapping_method": "raw_query_drug_fallback",
                            "intervention_mapping_confidence": 0.50,
                            "intervention_mapping_score": 0.50,
                            "intervention_token_jaccard_score": 0.50,
                            "intervention_mesh_id": "",
                            "nct_id": trial.get("nctId", ""),
                            "source_file": path.name,
                            "source_therapeutic_area": trial.get("sourceTherapeuticArea", ""),
                            "source_query_drug": source_query,
                            "fallback_reason": "queried_drug_not_mapped_to_mesh",
                        }
                    )

    fallback_df = pd.DataFrame(rows)
    fallback_path = new_processed_dir / "raw_query_drug_fallback_edges.csv"
    fallback_df.to_csv(fallback_path, index=False)
    return fallback_df


def dedupe_rows(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    out: List[Dict[str, Any]] = []
    for row in rows:
        key = (
            norm_entity(row.get("nct_id")),
            norm_entity(row.get("intervention")),
            norm_entity(row.get("condition")),
            norm_entity(row.get("raw_intervention")),
            norm_entity(row.get("raw_condition")),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def combine_processed_pairs(
    base_run_dir: Path,
    new_processed_dir: Path,
    output_processed_dir: Path,
    fallback_df: pd.DataFrame,
) -> Tuple[Path, pd.DataFrame]:
    base_rows = read_json(base_run_dir / "processed_data" / "condition_drug_pairs.json", [])
    new_rows = read_json(new_processed_dir / "condition_drug_pairs.json", [])
    if not isinstance(base_rows, list):
        base_rows = []
    if not isinstance(new_rows, list):
        new_rows = []

    fallback_rows = fallback_df.to_dict(orient="records") if not fallback_df.empty else []
    combined = dedupe_rows([*base_rows, *new_rows, *fallback_rows])
    output_path = output_processed_dir / "condition_drug_pairs.json"
    write_json(output_path, combined)

    summary = pd.DataFrame(
        [
            {"metric": "base_processed_rows", "value": len(base_rows)},
            {"metric": "new_drug_query_mapped_rows", "value": len(new_rows)},
            {"metric": "raw_query_drug_fallback_rows", "value": len(fallback_rows)},
            {"metric": "combined_deduped_rows", "value": len(combined)},
            {
                "metric": "combined_unique_pairs",
                "value": len({(norm_entity(r.get("intervention")), norm_entity(r.get("condition"))) for r in combined}),
            },
        ]
    )
    summary.to_csv(output_processed_dir.parent / "processed_pair_augmentation_summary.csv", index=False)
    return output_path, summary


def disease_variants(text: Any) -> List[str]:
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
    return list(dict.fromkeys(v for v in variants if v))


def resolve_entity_pair(
    requested_drug: str,
    requested_disease: str,
    combined_rows: List[Dict[str, Any]],
) -> Tuple[str, str, str]:
    drug_flat = flat_key(requested_drug)
    disease_flats = disease_variants(requested_disease)

    exact_pairs = []
    for row in combined_rows:
        intervention = str(row.get("intervention") or "").strip()
        condition = str(row.get("condition") or "").strip()
        raw_intervention = str(row.get("raw_intervention") or "")
        raw_condition = str(row.get("raw_condition") or "")
        intervention_text = f"{intervention} | {raw_intervention}"
        condition_text = f"{condition} | {raw_condition}"
        if drug_flat in flat_key(intervention_text) and any(v in flat_key(condition_text) for v in disease_flats):
            if intervention and condition:
                exact_pairs.append((intervention, condition))
    if exact_pairs:
        pair, _ = Counter(exact_pairs).most_common(1)[0]
        return pair[0], pair[1], "pair_level"

    drug_candidates = []
    disease_candidates = []
    for row in combined_rows:
        intervention = str(row.get("intervention") or "").strip()
        condition = str(row.get("condition") or "").strip()
        raw_intervention = str(row.get("raw_intervention") or "")
        raw_condition = str(row.get("raw_condition") or "")
        if intervention and drug_flat in flat_key(f"{intervention} | {raw_intervention}"):
            drug_candidates.append(intervention)
        if condition and any(v in flat_key(f"{condition} | {raw_condition}") for v in disease_flats):
            disease_candidates.append(condition)

    if drug_candidates and disease_candidates:
        drug, _ = Counter(drug_candidates).most_common(1)[0]
        disease, _ = Counter(disease_candidates).most_common(1)[0]
        return drug, disease, "entity_level"

    return requested_drug, requested_disease, "literal"


def build_graph_from_rows(rows: Iterable[Dict[str, Any]]) -> Tuple[nx.Graph, set[str], set[str], Dict[str, set[str]]]:
    graph = nx.Graph()
    drugs: set[str] = set()
    diseases: set[str] = set()
    known_indications: Dict[str, set[str]] = {}
    for row in rows:
        drug = norm_entity(row.get("intervention"))
        disease = norm_entity(row.get("condition"))
        if not drug or not disease or "placebo" in drug:
            continue
        graph.add_node(drug, bipartite="drug")
        graph.add_node(disease, bipartite="disease")
        graph.add_edge(drug, disease)
        drugs.add(drug)
        diseases.add(disease)
        known_indications.setdefault(drug, set()).add(disease)
    return graph, drugs, diseases, known_indications


def targeted_katz(graph: nx.Graph, drug: str, disease: str, alpha: float = 0.005) -> float:
    if not graph.has_node(drug) or not graph.has_node(disease):
        return 0.0
    score = alpha if graph.has_edge(drug, disease) else 0.0
    # Add a tiny length-3 contribution for shared intermediate structure.
    drug_neighbors = set(graph.neighbors(drug))
    disease_neighbors = set(graph.neighbors(disease))
    length3_paths = 0
    for node in drug_neighbors:
        length3_paths += len(set(graph.neighbors(node)) & disease_neighbors)
    return round(score + (alpha**3) * length3_paths, 6)


def compute_pair_features(
    graph: nx.Graph,
    known_indications: Dict[str, set[str]],
    degree: Dict[str, float],
    eigen: Dict[str, float],
    between: Dict[str, float],
    pageranks: Dict[str, Dict[str, float]],
    drug: str,
    disease: str,
) -> Dict[str, Any]:
    try:
        min_dist = min(
            nx.shortest_path_length(graph, source=disease, target=indication)
            for indication in known_indications.get(drug, [])
            if graph.has_node(disease) and nx.has_path(graph, disease, indication)
        )
        graph_distance_score = 1 / (1 + min_dist)
    except ValueError:
        graph_distance_score = 0.0

    centrality_drug = 0.33 * degree.get(drug, 0.0) + 0.33 * eigen.get(drug, 0.0) + 0.34 * between.get(drug, 0.0)
    centrality_disease = (
        0.33 * degree.get(disease, 0.0) + 0.33 * eigen.get(disease, 0.0) + 0.34 * between.get(disease, 0.0)
    )
    return {
        "Drug": drug,
        "Disease": disease,
        "GraphDistanceToIndication": round(graph_distance_score, 4),
        "RandomWalkScore": round(pageranks.get(drug, {}).get(disease, 0.0), 6),
        "StructuralLikelihood": round((1 + centrality_drug) * (1 + centrality_disease), 4),
        "PreferentialAttachment": round(degree.get(drug, 0.0) * degree.get(disease, 0.0), 6),
        "KatzSimilarity": targeted_katz(graph, drug, disease),
    }


def build_augmented_graph(
    matched_path: Path,
    graph_dir: Path,
    base_run_dir: Path,
    candidate_pairs: List[Tuple[str, str]],
    args: argparse.Namespace,
) -> pd.DataFrame:
    graph_dir.mkdir(parents=True, exist_ok=True)
    known_path = graph_dir / "graph_features_known.csv"
    unknown_path = graph_dir / "graph_features_unknown.csv"
    combined_rows = read_json(matched_path, [])
    combined_rows = combined_rows if isinstance(combined_rows, list) else []

    graph, drugs, diseases, known_indications = build_graph_from_rows(combined_rows)
    nx.write_graphml(graph, graph_dir / "bipartite.graphml")

    base_known_path = base_run_dir / "graph" / "graph_features_known.csv"
    base_unknown_path = base_run_dir / "graph" / "graph_features_unknown.csv"
    base_known_df = pd.read_csv(base_known_path) if base_known_path.exists() else pd.DataFrame()
    base_known_pairs = {
        (norm_entity(row["Drug"]), norm_entity(row["Disease"]))
        for _, row in base_known_df.iterrows()
        if "Drug" in row and "Disease" in row
    }

    combined_known_pairs = {
        (norm_entity(row.get("intervention")), norm_entity(row.get("condition")))
        for row in combined_rows
        if row.get("intervention") and row.get("condition")
    }
    requested_candidate_pairs = []
    for drug, disease in candidate_pairs:
        canonical_drug, canonical_disease, resolution = resolve_entity_pair(drug, disease, combined_rows)
        drug_key = norm_entity(canonical_drug)
        disease_key = norm_entity(canonical_disease)
        if drug_key in drugs and disease_key in diseases and (drug_key, disease_key) not in combined_known_pairs:
            requested_candidate_pairs.append((drug_key, disease_key, resolution))

    new_known_pairs = sorted(combined_known_pairs - base_known_pairs)
    feature_drugs = sorted({drug for drug, _ in new_known_pairs} | {drug for drug, _, _ in requested_candidate_pairs})

    degree = nx.degree_centrality(graph)
    try:
        eigen = nx.eigenvector_centrality(graph, max_iter=1000)
    except Exception:
        eigen = {node: 0.0 for node in graph.nodes()}
    between_sample = min(int(args.betweenness_sample), len(graph.nodes())) if args.betweenness_sample else None
    if between_sample and between_sample < len(graph.nodes()):
        between = nx.betweenness_centrality(graph, k=between_sample, seed=42)
    else:
        between = nx.betweenness_centrality(graph)
    pageranks = {
        drug: nx.pagerank(graph, alpha=0.85, personalization={drug: 1.0})
        for drug in feature_drugs
        if graph.has_node(drug)
    }

    new_known_rows = []
    for drug, disease in new_known_pairs:
        row = compute_pair_features(graph, known_indications, degree, eigen, between, pageranks, drug, disease)
        row["Label"] = 1
        new_known_rows.append(row)

    candidate_rows = []
    candidate_seen = set()
    for drug, disease, resolution in requested_candidate_pairs:
        if (drug, disease) in candidate_seen:
            continue
        candidate_seen.add((drug, disease))
        row = compute_pair_features(graph, known_indications, degree, eigen, between, pageranks, drug, disease)
        row["Label"] = 0
        candidate_rows.append(row)

    known_df = pd.concat([base_known_df, pd.DataFrame(new_known_rows)], ignore_index=True)
    known_df = known_df.drop_duplicates(["Drug", "Disease"], keep="first")
    known_df.to_csv(known_path, index=False)

    if base_unknown_path.exists():
        shutil.copy2(base_unknown_path, unknown_path)
        if candidate_rows:
            pd.DataFrame(candidate_rows).to_csv(unknown_path, mode="a", header=False, index=False)
    else:
        pd.DataFrame(candidate_rows).to_csv(unknown_path, index=False)

    unknown_row_count = sum(1 for _ in unknown_path.open("r", encoding="utf-8")) - 1 if unknown_path.exists() else 0
    audit = pd.DataFrame(
        [
            {"metric": "known_pair_rows", "value": len(known_df)},
            {"metric": "new_known_pair_rows", "value": len(new_known_rows)},
            {"metric": "target_candidate_rows_added", "value": len(candidate_rows)},
            {"metric": "unknown_pair_rows", "value": unknown_row_count},
            {"metric": "feature_rows", "value": len(known_df) + unknown_row_count},
            {"metric": "unique_drugs_in_augmented_graph", "value": len(drugs)},
            {"metric": "unique_diseases_in_augmented_graph", "value": len(diseases)},
        ]
    )
    audit.to_csv(graph_dir.parent / "augmented_graph_audit.csv", index=False)
    return audit


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
    return None


def parse_candidate_pairs(args: argparse.Namespace) -> List[Tuple[str, str]]:
    if args.pairs_file and args.pairs_file.exists():
        pairs: List[Tuple[str, str]] = []
        for line in args.pairs_file.read_text(encoding="utf-8").splitlines():
            parsed = parse_pair_line(line)
            if parsed:
                pairs.append(parsed)
        return pairs
    return list(TARGET_PAIRS)


def parse_drugs(args: argparse.Namespace) -> List[str]:
    if args.drugs_file:
        drugs = [
            line.strip()
            for line in args.drugs_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
    elif args.drug:
        drugs = args.drug
    elif args.pairs_file and args.pairs_file.exists():
        drugs = [drug for drug, _ in parse_candidate_pairs(args)]
    else:
        drugs = list(DEFAULT_DRUGS)
    return list(dict.fromkeys(drugs))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run-dir", type=Path, default=PROJECT_ROOT / "outputs" / "20260610_bayesian")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / "target_pair_graph_augmented_20260618")
    parser.add_argument("--drug", action="append", help="Drug name to query. Repeat for multiple drugs.")
    parser.add_argument("--drugs-file", type=Path)
    parser.add_argument("--pairs-file", type=Path)
    parser.add_argument("--query-parameter", default="query.intr")
    parser.add_argument("--fallback-query-term", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--refresh-raw", action="store_true")
    parser.add_argument("--page-size", type=int, default=100)
    parser.add_argument("--trials-per-drug", type=int)
    parser.add_argument("--cutoff-year", type=int, default=2026)
    parser.add_argument("--include-missing-dates", action="store_true")
    parser.add_argument("--mesh-path", type=Path, default=PROJECT_ROOT / "mesh_data" / "desc2026.xml")
    parser.add_argument("--fuzzy-cutoff", type=float, default=0.80)
    parser.add_argument("--token-jaccard-min", type=float, default=0.60)
    parser.add_argument("--betweenness-sample", type=int, default=256)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    raw_dir = output_dir / "raw_trials_drug_queries"
    new_processed_dir = output_dir / "processed_data_drug_queries"
    combined_processed_dir = output_dir / "processed_data"
    graph_dir = output_dir / "graph"
    for path in [raw_dir, new_processed_dir, combined_processed_dir, graph_dir]:
        path.mkdir(parents=True, exist_ok=True)

    drugs = parse_drugs(args)
    print(f"Drug queries: {', '.join(drugs)}")
    audit_df = fetch_all_drug_queries(drugs, raw_dir, args)

    mapping_audit_df, new_matched_path, new_unmatched_path = map_drug_query_trials_fast(
        raw_dir=raw_dir,
        processed_dir=new_processed_dir,
        mesh_path=args.mesh_path.resolve(),
        args=args,
    )
    mapping_audit_df.to_csv(output_dir / "drug_query_terminology_mapping_audit.csv", index=False)

    fallback_df = add_raw_query_drug_fallback_edges(
        raw_dir=raw_dir,
        new_processed_dir=new_processed_dir,
        mesh_path=args.mesh_path.resolve(),
        query_drugs=drugs,
        args=args,
    )
    combined_path, combine_summary = combine_processed_pairs(
        base_run_dir=args.base_run_dir.resolve(),
        new_processed_dir=new_processed_dir,
        output_processed_dir=combined_processed_dir,
        fallback_df=fallback_df,
    )

    graph_audit_df = build_augmented_graph(
        combined_path,
        graph_dir,
        args.base_run_dir.resolve(),
        parse_candidate_pairs(args),
        args,
    )
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "base_run_dir": rel_path(args.base_run_dir.resolve()),
        "output_dir": rel_path(output_dir),
        "drug_queries": drugs,
        "cutoff_year": args.cutoff_year,
        "query_parameter": args.query_parameter,
        "raw_trials_dir": rel_path(raw_dir),
        "new_matched_path": rel_path(new_matched_path),
        "new_unmatched_path": rel_path(new_unmatched_path),
        "combined_matched_path": rel_path(combined_path),
        "graph_dir": rel_path(graph_dir),
        "audit_files": [
            rel_path(output_dir / "clinicaltrials_drug_query_audit.csv"),
            rel_path(output_dir / "drug_query_terminology_mapping_audit.csv"),
            rel_path(output_dir / "processed_pair_augmentation_summary.csv"),
            rel_path(output_dir / "augmented_graph_audit.csv"),
        ],
        "notes": [
            "The archived manuscript run is not modified.",
            "Drug-query records use ClinicalTrials.gov intervention search by default.",
            "raw_query_drug_fallback edges are added only for queried drugs that fail MeSH drug mapping but have mapped conditions.",
        ],
    }
    write_json(output_dir / "clinicaltrials_drug_query_graph_manifest.json", manifest)

    print("\nClinicalTrials.gov drug-query audit")
    print(audit_df.to_string(index=False))
    print("\nProcessed-pair augmentation summary")
    print(combine_summary.to_string(index=False))
    print("\nAugmented graph audit")
    print(graph_audit_df.to_string(index=False))
    print(f"\nWrote augmented graph to: {rel_path(graph_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
