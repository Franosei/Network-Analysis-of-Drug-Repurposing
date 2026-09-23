"""
Run the ClinicalTrials.gov -> MeSH mapping -> graph -> quality-report pipeline.

This is an additive orchestration script. It leaves the original modules intact,
but writes a manuscript-style rerun folder with numbered audit outputs:

  01_clinical_trial_extraction_audit.csv
  02_terminology_mapping_audit.csv
  03_graph_construction_audit.csv
  04_literature_retrieval_audit.csv
  05_semantic_classification_audit.csv
  06_safety_overlap_audit.csv
  07_bayesian_uncertainty_audit.csv
  08_full_evidence_quality_ledger.csv
  09_quality_category_counts.csv
  10_summary_dashboard.csv
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import re
import shutil
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import requests

from data_quality.metric_definitions import LEDGER_COLUMNS, METRIC_DEFINITIONS
from data_quality.quality_flags import coverage_tier as rule_coverage_tier


PROJECT_ROOT = Path(__file__).resolve().parent
CODE_DIR = PROJECT_ROOT / "code"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from condition_drug_pairs import ConditionDrugPairBuilder  # noqa: E402
from data_extraction import ClinicalTrialFetcher  # noqa: E402
from evidence_quality_report import (  # noqa: E402
    beta_metrics,
    clamp01,
    generate_reports,
    norm_entity,
    read_json,
    round_or_none,
    timestamp_from_path,
)
from scalable_graph_builder import build_targeted_full_graph  # noqa: E402


DEFAULT_AREAS = [
    "cancer",
    "diabetes",
    "asthma",
    "Alzheimer's disease",
    "depression",
    "hypertension",
    "HIV/AIDS",
    "rheumatoid arthritis",
    "epilepsy",
    "obesity",
]

MESH_2026_URL = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.xml"
DEFAULT_PANEL_PATH = "config/case_study_panel_publication.csv"
# Legacy fallback — used when publication panel is absent
_LEGACY_PANEL_PATH = "config/case_study_panel.csv"
DEFAULT_WEIGHTS = {
    "GraphDistanceToIndication": 0.1148,
    "RandomWalkScore": 0.247,
    "StructuralLikelihood": -0.1154,
    "PreferentialAttachment": -0.041,
    "KatzSimilarity": 1.6515,
}

# ── Case study definitions ────────────────────────────────────────────────────
# Three classes: successful, failed/harmful, emerging.
# Narrative explains the scientific rationale for each selection.
CASE_STUDIES: List[Dict[str, Any]] = [
    {
        "drug": "dexamethasone",
        "disease": "covid-19",
        "case_type": "successful_repurposing",
        "expected_quality": "High evidence quality",
        "narrative": (
            "Dexamethasone-COVID-19 (RECOVERY trial). Corticosteroid reduced 28-day mortality "
            "in patients requiring oxygen/ventilation. Strong RCT evidence; adverse signal "
            "reflects steroid toxicity in mild disease where it is contra-indicated. "
            "Model of evidence-quality-supported repurposing: therapeutic and adverse "
            "signals can be correctly stratified."
        ),
    },
    {
        "drug": "hydroxychloroquine",
        "disease": "covid-19",
        "case_type": "failed_repurposing",
        "expected_quality": "Safety-conflicted evidence",
        "narrative": (
            "Hydroxychloroquine-COVID-19. Early in-vitro enthusiasm; multiple large RCTs "
            "(SOLIDARITY, RECOVERY, WHO Solidarity) showed no clinical benefit and "
            "raised cardiac-arrhythmia concerns. Demonstrates how high early literature "
            "volume with safety conflict produces a misleading repurposing signal when "
            "evidence quality is not audited."
        ),
    },
    {
        "drug": "ivermectin",
        "disease": "covid-19",
        "case_type": "weak_evidence_high_attention",
        "expected_quality": "Literature noise dominated",
        "narrative": (
            "Ivermectin-COVID-19. Massive public attention; major RCTs (TOGETHER, ACTIV-6, "
            "PRINCIPLE) showed no clinically meaningful benefit. Illustrates the key "
            "thesis: high literature volume does not equal high evidence quality when "
            "evidence is noisy, conflicting, or methodologically weak."
        ),
    },
    {
        "drug": "colchicine",
        "disease": "coronary artery disease",
        "case_type": "emerging_repurposing",
        "expected_quality": "Moderate evidence quality",
        "narrative": (
            "Colchicine-CAD. Anti-inflammatory mechanism plausible via NLRP3 inflammasome "
            "inhibition. COLCOT and LoDoCo2 trials showed event reduction post-MI. "
            "FDA-approved (2023) for CV risk reduction. Model of emerging repurposing: "
            "growing evidence maturity but monitoring for GI and pulmonary adverse effects."
        ),
    },
    {
        "drug": "metformin",
        "disease": "uterine cervical dysplasia",
        "case_type": "novel_hypothesis",
        "expected_quality": "Insufficient evidence",
        "narrative": (
            "Metformin-cervical dysplasia. mTOR/AMPK mechanistic hypothesis; limited "
            "clinical trial evidence. Model of how the framework identifies 'not yet "
            "evidence-ready' candidates that need further investigation before prioritisation."
        ),
    },
    {
        "drug": "metformin",
        "disease": "mitochondrial diseases",
        "case_type": "novel_hypothesis",
        "expected_quality": "Insufficient evidence",
        "narrative": (
            "Metformin-mitochondrial diseases. AMPK activation and mitophagy induction "
            "proposed; small case series only. Shows how structural network proximity "
            "can surface mechanistically plausible but evidence-sparse candidates."
        ),
    },
]


def str_to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"true", "1", "yes", "y"}:
        return True
    if value in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def parse_areas(value: str) -> List[str]:
    if not value:
        return list(DEFAULT_AREAS)
    return [item.strip() for item in value.split(",") if item.strip()]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def _package_version(distribution: str) -> str:
    """Return an installed distribution version without making it a hard dependency."""
    try:
        from importlib.metadata import version

        return version(distribution)
    except Exception:
        return "unavailable"


def ensure_mesh(mesh_path: Path, download_mesh: bool, log_lines: List[str]) -> None:
    if mesh_path.exists():
        log_lines.append(f"MeSH file found: {rel(mesh_path)}")
        return

    if not download_mesh:
        raise FileNotFoundError(
            f"MeSH descriptor file not found at {mesh_path}. "
            "Run with --download_mesh true or place desc2026.xml there."
        )

    mesh_path.parent.mkdir(parents=True, exist_ok=True)
    log_lines.append(f"Downloading MeSH descriptors from {MESH_2026_URL}")
    with requests.get(MESH_2026_URL, stream=True, timeout=60) as response:
        response.raise_for_status()
        with mesh_path.open("wb") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    log_lines.append(f"Downloaded MeSH file to {rel(mesh_path)}")


def load_mesh_descriptor_ids(mesh_path: Path) -> Dict[str, str]:
    """Return descriptor-name lookup to MeSH descriptor UI."""
    out: Dict[str, str] = {}
    if not mesh_path.exists():
        return out
    root = ET.parse(mesh_path).getroot()
    for descriptor in root.findall("DescriptorRecord"):
        canonical = (descriptor.findtext("DescriptorName/String") or "").strip()
        mesh_id = (descriptor.findtext("DescriptorUI") or "").strip()
        if not canonical or not mesh_id:
            continue
        keys = {
            canonical,
            canonical.lower(),
            re.sub(r"\s+", " ", canonical.lower()).strip(),
            re.sub(r"[^a-z0-9]+", "", canonical.lower()),
        }
        for key in keys:
            out[key] = mesh_id
    return out


def mesh_id_for(term: Any, mesh_ids: Dict[str, str]) -> str:
    text = str(term or "").strip()
    if not text:
        return ""
    return (
        mesh_ids.get(text)
        or mesh_ids.get(text.lower())
        or mesh_ids.get(re.sub(r"\s+", " ", text.lower()).strip())
        or mesh_ids.get(re.sub(r"[^a-z0-9]+", "", text.lower()))
        or ""
    )


def load_case_study_panel(path: Path) -> pd.DataFrame:
    if path.exists():
        panel = pd.read_csv(path)
    else:
        # Try legacy panel before falling back to hard-coded list
        legacy = path.parent / _LEGACY_PANEL_PATH.split("/")[-1]
        if legacy.exists():
            panel = pd.read_csv(legacy)
        else:
            panel = pd.DataFrame(CASE_STUDIES)
    required = {"drug", "disease"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"Case-study panel is missing required columns: {sorted(missing)}")
    panel["drug"] = panel["drug"].astype(str).str.strip()
    panel["disease"] = panel["disease"].astype(str).str.strip()
    panel = panel[(panel["drug"] != "") & (panel["disease"] != "")].copy()
    panel["case_type"] = panel.get("case_type", "panel_pair")
    panel["expected_quality"] = panel.get("expected_quality", "")
    panel["priority"] = panel.get("priority", 99)
    panel["narrative"] = panel.get("narrative", "")
    panel = panel.drop_duplicates(subset=["drug", "disease"], keep="first")
    return panel.reset_index(drop=True)


def trial_filter_flags(fetcher: ClinicalTrialFetcher, study: Dict[str, Any]) -> Dict[str, Any]:
    protocol = study.get("protocolSection", {}) or {}
    status = str(protocol.get("statusModule", {}).get("overallStatus", "")).upper()
    phases = {str(p).upper() for p in (protocol.get("designModule", {}).get("phases", []) or [])}
    interventions = protocol.get("armsInterventionsModule", {}).get("interventions", []) or []
    conditions = protocol.get("conditionsModule", {}).get("conditions", []) or []
    has_drug = any(str(i.get("type", "")).upper() == fetcher.valid_intervention_type for i in interventions)
    year, year_source = fetcher._extract_reference_year(study)

    return {
        "status": status,
        "phases": sorted(phases),
        "has_valid_status": status in fetcher.valid_statuses,
        "has_valid_phase": bool(phases.intersection(fetcher.valid_phases)),
        "has_drug_intervention": has_drug,
        "within_cutoff": fetcher.is_within_cutoff(study),
        "reference_year": year,
        "reference_year_source": year_source,
        "missing_status": not bool(status),
        "missing_phase": not bool(phases),
        "missing_intervention": not bool(interventions),
        "missing_condition": not bool(conditions),
    }


def fetch_area_with_audit(
    fetcher: ClinicalTrialFetcher,
    area: str,
    trials_per_area: Optional[int],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    trials: List[Dict[str, Any]] = []
    next_page_token: Optional[str] = None
    audit = Counter()

    while True:
        params = {"query.term": area, "pageSize": fetcher.page_size}
        if next_page_token:
            params["pageToken"] = next_page_token

        result = fetcher._get_with_retries(fetcher.base_url, params=params)
        if result is None:
            audit["api_failures"] += 1
            break

        audit["api_pages"] += 1
        studies = result.get("studies", []) or []
        audit["api_records_seen"] += len(studies)
        next_page_token = result.get("nextPageToken")

        for study in studies:
            flags = trial_filter_flags(fetcher, study)
            audit["missing_status"] += int(flags["missing_status"])
            audit["missing_phase"] += int(flags["missing_phase"])
            audit["missing_intervention"] += int(flags["missing_intervention"])
            audit["missing_condition"] += int(flags["missing_condition"])

            if not flags["has_valid_status"]:
                audit["excluded_status"] += 1
                continue
            if not flags["has_valid_phase"]:
                audit["excluded_phase"] += 1
                continue
            if not flags["has_drug_intervention"]:
                audit["excluded_no_drug_intervention"] += 1
                continue
            if not flags["within_cutoff"]:
                audit["excluded_after_cutoff_or_missing_date"] += 1
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
            extracted["sourceTherapeuticArea"] = area
            trials.append(extracted)
            fetcher.seen_nct_ids.add(nct_id)
            audit["eligible_saved"] += 1

            if trials_per_area is not None and len(trials) >= int(trials_per_area):
                audit["stopped_at_trials_per_area_cap"] = 1
                audit_row = dict(audit)
                audit_row["therapeutic_area"] = area
                return trials, audit_row

        time.sleep(random.uniform(0.2, 0.6))
        if not next_page_token:
            break

    audit_row = dict(audit)
    audit_row["therapeutic_area"] = area
    return trials, audit_row


def fetch_trials_stage(args: argparse.Namespace, raw_dir: Path, log_lines: List[str]) -> pd.DataFrame:
    raw_dir.mkdir(parents=True, exist_ok=True)
    areas = parse_areas(args.therapeutic_areas)
    snapshot_time = datetime.now(timezone.utc).isoformat()

    fetcher = ClinicalTrialFetcher(
        output_dir=rel(raw_dir),
        page_size=args.page_size,
        cutoff_year=args.cutoff_year,
        include_missing_dates=args.include_missing_dates,
    )

    audit_rows: List[Dict[str, Any]] = []
    for area in areas:
        log_lines.append(f"Fetching ClinicalTrials.gov records for {area}")
        trials, audit = fetch_area_with_audit(fetcher, area, args.trials_per_area)
        audit_rows.append(audit)

        filename = fetcher.sanitize_filename(area) + f"_upto_{args.cutoff_year}.json"
        output_path = raw_dir / filename
        with output_path.open("w", encoding="utf-8") as f:
            json.dump(trials, f, indent=2)
        log_lines.append(f"Saved {len(trials)} eligible trials to {rel(output_path)}")

    audit_df = pd.DataFrame(audit_rows)
    audit_df["registry_snapshot_timestamp_utc"] = snapshot_time
    audit_df["historical_snapshot_status"] = "current_record_filtered_by_date_not_historical"
    audit_df["cutoff_semantics"] = (
        "study fields were retrieved at snapshot time and filtered by available dates; "
        "later amendments may be present"
    )
    metadata_cols = {
        "therapeutic_area",
        "registry_snapshot_timestamp_utc",
        "historical_snapshot_status",
        "cutoff_semantics",
    }
    numeric_cols = [col for col in audit_df.columns if col not in metadata_cols]
    total_row = {
        "therapeutic_area": "ALL_AREAS",
        "registry_snapshot_timestamp_utc": snapshot_time,
        "historical_snapshot_status": "current_record_filtered_by_date_not_historical",
        "cutoff_semantics": (
            "study fields were retrieved at snapshot time and filtered by available dates; "
            "later amendments may be present"
        ),
    }
    for col in numeric_cols:
        total_row[col] = pd.to_numeric(audit_df[col], errors="coerce").fillna(0).sum()
    audit_df = pd.concat([audit_df, pd.DataFrame([total_row])], ignore_index=True)
    return audit_df


def fetch_full_registry_stage(
    args: argparse.Namespace, raw_dir: Path, log_lines: List[str]
) -> pd.DataFrame:
    """Download/resume the whole registry and export drug trials in JSON chunks."""
    from reporting.run_exact_pair_full_registry import (  # noqa: PLC0415
        download_registry,
        initialise_database,
        metadata_get,
        unpack_metadata,
    )

    raw_dir.mkdir(parents=True, exist_ok=True)
    snapshot_time = datetime.now(timezone.utc).isoformat()
    configured = str(args.full_registry_db or "").strip()
    database_path = (
        (PROJECT_ROOT / configured).resolve()
        if configured
        else raw_dir / "clinicaltrials_full_registry.sqlite"
    )
    # A completed registry database is a snapshot, not a timeless source.  By
    # default start a fresh acquisition so an old database cannot silently
    # masquerade as the current registry.  Resumption is explicit.
    connection = initialise_database(
        database_path, refresh=not bool(getattr(args, "resume_full_registry", False))
    )
    try:
        session = requests.Session()
        session.headers.update(
            {"User-Agent": "full-registry-bayesian-pipeline/2.0 (research audit)"}
        )
        registry_audit = download_registry(
            connection,
            session,
            page_size=max(1, int(args.page_size)),
            max_pages=None,
            pause=0.05,
        )
        complete = unpack_metadata(metadata_get(connection, "registry_complete"), False)
        if not complete:
            raise RuntimeError("ClinicalTrials.gov full-registry snapshot is incomplete")

        for old_chunk in raw_dir.glob("full_registry_drug_trials_*.json"):
            old_chunk.unlink()
        chunk_size = max(1, int(args.registry_export_chunk_size))
        chunk: List[Dict[str, Any]] = []
        chunk_number = 0
        exported = 0
        current_nct = None
        current: Optional[Dict[str, Any]] = None

        def flush_study() -> None:
            nonlocal current, chunk, chunk_number, exported
            if current is None:
                return
            current["conditions"] = sorted(current.pop("_conditions"))
            current["interventions"] = sorted(current.pop("_interventions"))
            chunk.append(current)
            exported += 1
            current = None
            if len(chunk) >= chunk_size:
                chunk_number += 1
                path = raw_dir / f"full_registry_drug_trials_{chunk_number:04d}.json"
                path.write_text(json.dumps(chunk, ensure_ascii=False), encoding="utf-8")
                chunk = []
                print(f"Exported {exported:,} full-registry drug trials", flush=True)

        query = """
            SELECT p.nct_id, s.overall_status, s.phases_json, s.start_date,
                   s.first_posted, p.raw_disease, p.raw_drug
            FROM study_pairs p
            JOIN studies s ON s.nct_id = p.nct_id
            ORDER BY p.nct_id
        """
        for nct_id, status, phases_json, start_date, first_posted, disease, drug in connection.execute(query):
            if nct_id != current_nct:
                flush_study()
                current_nct = nct_id
                current = {
                    "nctId": nct_id,
                    "status": status,
                    "phases": json.loads(phases_json or "[]"),
                    "startDate": start_date,
                    "studyFirstPostDate": first_posted,
                    "sourceTherapeuticArea": "FULL_REGISTRY",
                    "_conditions": set(),
                    "_interventions": set(),
                }
            current["_conditions"].add(str(disease))
            current["_interventions"].add(str(drug))
        flush_study()
        if chunk:
            chunk_number += 1
            path = raw_dir / f"full_registry_drug_trials_{chunk_number:04d}.json"
            path.write_text(json.dumps(chunk, ensure_ascii=False), encoding="utf-8")

        log_lines.append(
            f"Full ClinicalTrials.gov snapshot: {database_path}; exported {exported} drug trials in {chunk_number} chunks"
        )
        return pd.DataFrame(
            [
                {
                    "therapeutic_area": "FULL_REGISTRY",
                    "registry_database": rel(database_path),
                    "registry_complete": True,
                    "registry_studies": registry_audit.get("studies_returned", 0),
                    "interventional_with_drug": registry_audit.get("interventional_with_drug", 0),
                    "study_pair_rows": registry_audit.get("study_pair_rows", 0),
                    "exported_drug_trials": exported,
                    "export_chunks": chunk_number,
                    "registry_snapshot_timestamp_utc": snapshot_time,
                    "historical_snapshot_status": "complete_current_registry_snapshot_not_historical_archive",
                }
            ]
        )
    finally:
        connection.close()


def classify_mapping(method: Optional[str], score: Optional[float] = None) -> Tuple[str, float]:
    method_norm = norm_entity(method)
    if method_norm == "exact":
        return "exact_match", 1.0
    if method_norm == "fuzzy":
        confidence = float(score or 0.0)
        if confidence >= 0.80:
            return "fuzzy_high_confidence", confidence
        return "fuzzy_low_confidence", confidence
    if method_norm == "token_score":
        confidence = float(score or 0.0)
        if confidence >= 0.80:
            return "fuzzy_high_confidence", confidence
        return "fuzzy_low_confidence", confidence
    if method_norm in {"empty", "no_tokens", "no_candidates", "token_score_failed"}:
        return "unmapped", 0.0
    return "unmapped", 0.0


def debug_score(debug: Dict[str, Any], default: float = 0.0) -> float:
    try:
        return float(debug.get("score", default))
    except Exception:
        return default


def debug_token_jaccard(debug: Dict[str, Any], default: float = 0.0) -> float:
    try:
        return float(debug.get("token_jaccard_score", default))
    except Exception:
        return default


def split_trial_items(builder: ConditionDrugPairBuilder, items: Any) -> Tuple[List[str], int]:
    if not isinstance(items, list):
        return [], 0

    result: List[str] = []
    placebo_count = 0
    for item in items:
        if not item:
            continue
        # Preserve conjunctions, punctuation and slash-containing biomedical
        # names. ClinicalTrials.gov already supplies list elements.
        parts = [str(item)]
        for part in parts:
            part = part.strip()
            if not part:
                continue
            if not builder.include_placebo and part.lower() == "placebo":
                placebo_count += 1
                continue
            result.append(part)
    return result, placebo_count


def mapping_stage(
    raw_dir: Path,
    processed_dir: Path,
    mesh_path: Path,
    args: argparse.Namespace,
    log_lines: List[str],
) -> Tuple[pd.DataFrame, Path, Path]:
    processed_dir.mkdir(parents=True, exist_ok=True)
    matched_path = processed_dir / "condition_drug_pairs.json"
    unmatched_path = processed_dir / "unmatched_pairs.json"

    builder = ConditionDrugPairBuilder(
        input_dir=rel(raw_dir),
        output_dir=rel(processed_dir),
        output_filename=matched_path.name,
        unmatched_filename=unmatched_path.name,
        mesh_path=rel(mesh_path),
        fuzzy_cutoff=args.fuzzy_cutoff,
        token_jaccard_min=args.token_jaccard_min,
        include_placebo=False,
        keep_unmatched_debug_fields=True,
        build_token_indexes=not args.mesh_exact_only,
    )
    mesh_ids = load_mesh_descriptor_ids(mesh_path)

    condition_cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]] = {}
    drug_cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]] = {}

    def cached_match(
        raw_term: str,
        *,
        mesh_map: Dict[str, str],
        token_index: Dict[str, List[str]],
        cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]],
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        cache_key = str(raw_term)
        if cache_key in cache:
            return cache[cache_key]
        if args.mesh_exact_only:
            normalized = builder.normalize_for_match(raw_term)
            matched = mesh_map.get(normalized)
            debug = (
                {
                    "raw_term": str(raw_term),
                    "method": "exact",
                    "normalized": normalized,
                    "matched_key": normalized,
                    "score": 1.0,
                    "sequence_ratio": 1.0,
                    "token_jaccard_score": 1.0,
                    "removed_modifiers": [],
                    "ambiguity_flag": False,
                }
                if matched
                else {
                    "raw_term": str(raw_term),
                    "method": "no_exact_mesh_term",
                    "normalized": normalized,
                    "score": 0.0,
                    "removed_modifiers": [],
                    "ambiguity_flag": False,
                }
            )
        else:
            matched, debug_value = builder.match_term(
                raw_term=raw_term,
                mesh_map=mesh_map,
                token_index=token_index,
            )
            debug = debug_value or {}
        cache[cache_key] = (matched, debug)
        return cache[cache_key]

    rows: List[Dict[str, Any]] = []
    unmatched: List[Dict[str, Any]] = []
    audit_rows: List[Dict[str, Any]] = []

    for path in sorted(raw_dir.glob("*.json")):
        trials = read_json(path, [])
        if not isinstance(trials, list):
            continue

        for trial in trials:
            if args.full_registry:
                conditions = [
                    str(value).strip()
                    for value in trial.get("conditions", [])
                    if str(value).strip()
                ]
                raw_interventions = [
                    str(value).strip()
                    for value in trial.get("interventions", [])
                    if str(value).strip()
                ]
                interventions = [
                    value for value in raw_interventions if "placebo" not in value.casefold()
                ]
                placebo_count = len(raw_interventions) - len(interventions)
            else:
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
                cleaned_condition = builder.normalize_for_match(raw_condition)
                condition_match, condition_debug = cached_match(
                    raw_term=raw_condition,
                    mesh_map=builder.condition_term_map,
                    token_index=builder._cond_token_index,
                    cache=condition_cache,
                )
                condition_debug = condition_debug or {}
                condition_method, condition_confidence = classify_mapping(
                    condition_debug.get("method"),
                    condition_debug.get("score"),
                )
                mapped_condition = builder.clean_output_name(condition_match) if condition_match else ""
                condition_mesh_id = mesh_id_for(mapped_condition, mesh_ids)
                condition_mapping_score = debug_score(condition_debug, condition_confidence)
                condition_token_jaccard = debug_token_jaccard(condition_debug)

                if not condition_match:
                    unmatched_entry = {
                        "raw_condition": raw_condition,
                        "raw_intervention": None,
                        "reason": "unmatched condition",
                        "debug": condition_debug,
                        "nct_id": trial.get("nctId", ""),
                        "source_file": path.name,
                    }
                    unmatched.append(unmatched_entry)
                    audit_rows.append(
                        {
                            **base,
                            "original_term": raw_condition,
                            "cleaned_term": cleaned_condition,
                            "mapped_term": "",
                            "mesh_id": "",
                            "entity_type": "disease",
                            "mapping_method": condition_method,
                            "mapping_score": condition_mapping_score,
                            "token_jaccard_score": condition_token_jaccard,
                            "mapping_status": "unmapped",
                            "failure_reason": "unmatched condition",
                            "original_disease_term": raw_condition,
                            "disease_cleaned_term": cleaned_condition,
                            "raw_disease": raw_condition,
                            "mapped_disease": "",
                            "disease_mesh_id": "",
                            "disease_mapping_method": condition_method,
                            "disease_mapping_confidence": condition_confidence,
                            "disease_mapping_score": condition_mapping_score,
                            "disease_token_jaccard_score": condition_token_jaccard,
                            "disease_mapping_status": "unmapped",
                            "raw_drug": "",
                            "original_drug_term": "",
                            "drug_cleaned_term": "",
                            "mapped_drug": "",
                            "drug_mesh_id": "",
                            "drug_mapping_method": "not_attempted",
                            "drug_mapping_confidence": 0.0,
                            "drug_mapping_score": 0.0,
                            "drug_token_jaccard_score": 0.0,
                            "drug_mapping_status": "not_attempted",
                            "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                            "drug_debug": "",
                        }
                    )
                    continue

                for raw_drug in interventions:
                    cleaned_drug = builder.normalize_for_match(raw_drug)
                    drug_match, drug_debug = cached_match(
                        raw_term=raw_drug,
                        mesh_map=builder.drug_term_map,
                        token_index=builder._drug_token_index,
                        cache=drug_cache,
                    )
                    drug_debug = drug_debug or {}
                    drug_method, drug_confidence = classify_mapping(
                        drug_debug.get("method"),
                        drug_debug.get("score"),
                    )
                    mapped_drug = builder.clean_output_name(drug_match) if drug_match else ""
                    drug_mesh_id = mesh_id_for(mapped_drug, mesh_ids)
                    drug_mapping_score = debug_score(drug_debug, drug_confidence)
                    drug_token_jaccard = debug_token_jaccard(drug_debug)

                    if not drug_match:
                        unmatched_entry = {
                            "raw_condition": raw_condition,
                            "mapped_condition": builder.clean_output_name(condition_match),
                            "raw_intervention": raw_drug,
                            "reason": "unmatched intervention",
                            "debug": drug_debug,
                            "nct_id": trial.get("nctId", ""),
                            "source_file": path.name,
                        }
                        unmatched.append(unmatched_entry)
                        audit_rows.append(
                            {
                                **base,
                                "original_term": raw_drug,
                                "cleaned_term": cleaned_drug,
                                "mapped_term": "",
                                "mesh_id": "",
                                "entity_type": "drug",
                                "mapping_method": drug_method,
                                "mapping_score": drug_mapping_score,
                                "token_jaccard_score": drug_token_jaccard,
                                "mapping_status": "unmapped",
                                "failure_reason": "unmatched intervention",
                                "original_disease_term": raw_condition,
                                "disease_cleaned_term": cleaned_condition,
                                "raw_disease": raw_condition,
                                "mapped_disease": mapped_condition,
                                "disease_mesh_id": condition_mesh_id,
                                "disease_mapping_method": condition_method,
                                "disease_mapping_confidence": condition_confidence,
                                "disease_mapping_score": condition_mapping_score,
                                "disease_token_jaccard_score": condition_token_jaccard,
                                "disease_mapping_status": "mapped",
                                "raw_drug": raw_drug,
                                "original_drug_term": raw_drug,
                                "drug_cleaned_term": cleaned_drug,
                                "mapped_drug": "",
                                "drug_mesh_id": "",
                                "drug_mapping_method": drug_method,
                                "drug_mapping_confidence": drug_confidence,
                                "drug_mapping_score": drug_mapping_score,
                                "drug_token_jaccard_score": drug_token_jaccard,
                                "drug_mapping_status": "unmapped",
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
                        "condition_cleaned_term": cleaned_condition,
                        "condition_mapping_method": condition_method,
                        "condition_mapping_confidence": condition_confidence,
                        "condition_mapping_score": condition_mapping_score,
                        "condition_token_jaccard_score": condition_token_jaccard,
                        "condition_mesh_id": condition_mesh_id,
                        "intervention_cleaned_term": cleaned_drug,
                        "intervention_mapping_method": drug_method,
                        "intervention_mapping_confidence": drug_confidence,
                        "intervention_mapping_score": drug_mapping_score,
                        "intervention_token_jaccard_score": drug_token_jaccard,
                        "intervention_mesh_id": drug_mesh_id,
                        "nct_id": trial.get("nctId", ""),
                        "source_file": path.name,
                    }
                    rows.append(matched_row)
                    audit_rows.append(
                        {
                            **base,
                            "original_term": f"{raw_drug} | {raw_condition}",
                            "cleaned_term": f"{cleaned_drug} | {cleaned_condition}",
                            "mapped_term": f"{mapped_drug} | {mapped_condition}",
                            "mesh_id": f"{drug_mesh_id} | {condition_mesh_id}",
                            "entity_type": "drug_disease_pair",
                            "mapping_method": "matched",
                            "mapping_score": min(drug_mapping_score, condition_mapping_score),
                            "token_jaccard_score": min(drug_token_jaccard, condition_token_jaccard),
                            "original_disease_term": raw_condition,
                            "disease_cleaned_term": cleaned_condition,
                            "raw_disease": raw_condition,
                            "mapped_disease": mapped_condition,
                            "disease_mesh_id": condition_mesh_id,
                            "disease_mapping_method": condition_method,
                            "disease_mapping_confidence": condition_confidence,
                            "disease_mapping_score": condition_mapping_score,
                            "disease_token_jaccard_score": condition_token_jaccard,
                            "disease_mapping_status": "mapped",
                            "raw_drug": raw_drug,
                            "original_drug_term": raw_drug,
                            "drug_cleaned_term": cleaned_drug,
                            "mapped_drug": mapped_drug,
                            "drug_mesh_id": drug_mesh_id,
                            "drug_mapping_method": drug_method,
                            "drug_mapping_confidence": drug_confidence,
                            "drug_mapping_score": drug_mapping_score,
                            "drug_token_jaccard_score": drug_token_jaccard,
                            "drug_mapping_status": "mapped",
                            "mapping_status": "matched",
                            "failure_reason": "",
                            "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                            "drug_debug": json.dumps(drug_debug, ensure_ascii=False),
                        }
                    )

    with matched_path.open("w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2, ensure_ascii=False)
    with unmatched_path.open("w", encoding="utf-8") as f:
        json.dump(unmatched, f, indent=2, ensure_ascii=False)

    log_lines.append(f"Saved {len(rows)} matched drug-disease rows to {rel(matched_path)}")
    log_lines.append(f"Saved {len(unmatched)} unmatched mapping rows to {rel(unmatched_path)}")
    log_lines.append(
        f"MeSH unique terms evaluated: {len(condition_cache)} diseases, {len(drug_cache)} drugs; "
        f"exact_only={args.mesh_exact_only}"
    )
    return pd.DataFrame(audit_rows), matched_path, unmatched_path


def mapping_stage_full_registry(
    raw_dir: Path,
    processed_dir: Path,
    mesh_path: Path,
    args: argparse.Namespace,
    log_lines: List[str],
) -> Tuple[pd.DataFrame, Path, Path]:
    """Memory-bounded high-precision MeSH mapping for full registry chunks."""
    processed_dir.mkdir(parents=True, exist_ok=True)
    matched_path = processed_dir / "condition_drug_pairs.json"
    unmatched_path = processed_dir / "unmatched_pairs.json"
    builder = ConditionDrugPairBuilder(
        input_dir=rel(raw_dir),
        output_dir=rel(processed_dir),
        output_filename=matched_path.name,
        unmatched_filename=unmatched_path.name,
        mesh_path=rel(mesh_path),
        fuzzy_cutoff=args.fuzzy_cutoff,
        token_jaccard_min=args.token_jaccard_min,
        include_placebo=False,
        keep_unmatched_debug_fields=True,
        build_token_indexes=not args.mesh_exact_only,
    )
    mesh_ids = load_mesh_descriptor_ids(mesh_path)
    condition_cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]] = {}
    drug_cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]] = {}

    def map_one(
        raw: str,
        mesh_map: Dict[str, str],
        token_index: Dict[str, List[str]],
        cache: Dict[str, Tuple[Optional[str], Dict[str, Any]]],
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        if raw in cache:
            return cache[raw]
        normalized = builder.normalize_for_match(raw)
        if args.mesh_exact_only:
            matched = mesh_map.get(normalized)
            debug = {
                "method": "exact" if matched else "no_exact_mesh_term",
                "normalized": normalized,
                "matched_key": normalized if matched else "",
                "score": 1.0 if matched else 0.0,
                "token_jaccard_score": 1.0 if matched else 0.0,
            }
        else:
            matched, detail = builder.match_term(raw, mesh_map, token_index)
            debug = detail or {}
        cache[raw] = (matched, debug)
        return cache[raw]

    matched_rows: List[Dict[str, Any]] = []
    entity_audit: Dict[Tuple[str, str], Dict[str, Any]] = {}
    unmatched_entities: Dict[Tuple[str, str], Dict[str, Any]] = {}
    pair_attempts = 0
    placebo_excluded = 0
    trials_processed = 0

    chunk_paths = sorted(raw_dir.glob("full_registry_drug_trials_*.json"))
    if not chunk_paths:
        raise FileNotFoundError(f"No full-registry JSON chunks found under {raw_dir}")
    for chunk_index, path in enumerate(chunk_paths, start=1):
        trials = read_json(path, [])
        for trial in trials if isinstance(trials, list) else []:
            trials_processed += 1
            raw_conditions = sorted(
                {str(value).strip() for value in trial.get("conditions", []) if str(value).strip()}
            )
            all_drugs = sorted(
                {str(value).strip() for value in trial.get("interventions", []) if str(value).strip()}
            )
            raw_drugs = [value for value in all_drugs if "placebo" not in value.casefold()]
            placebo_excluded += len(all_drugs) - len(raw_drugs)
            mapped_conditions = []
            for raw in raw_conditions:
                mapped, debug = map_one(
                    raw,
                    builder.condition_term_map,
                    builder._cond_token_index,
                    condition_cache,
                )
                method, confidence = classify_mapping(debug.get("method"), debug.get("score"))
                entity_audit.setdefault(
                    ("disease", raw),
                    {
                        "entity_type": "disease",
                        "original_term": raw,
                        "cleaned_term": debug.get("normalized", ""),
                        "mapped_term": mapped or "",
                        "mesh_id": mesh_id_for(mapped, mesh_ids) if mapped else "",
                        "mapping_method": method,
                        "mapping_score": confidence,
                        "mapping_status": "mapped" if mapped else "unmapped",
                    },
                )
                if mapped:
                    mapped_conditions.append(
                        (
                            raw,
                            builder.clean_output_name(mapped),
                            mesh_id_for(mapped, mesh_ids),
                            debug,
                            method,
                            confidence,
                        )
                    )
                else:
                    unmatched_entities[("disease", raw)] = {
                        "entity_type": "disease",
                        "raw_term": raw,
                        "reason": "unmatched condition",
                        "debug": debug,
                    }

            mapped_drugs = []
            for raw in raw_drugs:
                mapped, debug = map_one(
                    raw,
                    builder.drug_term_map,
                    builder._drug_token_index,
                    drug_cache,
                )
                method, confidence = classify_mapping(debug.get("method"), debug.get("score"))
                entity_audit.setdefault(
                    ("drug", raw),
                    {
                        "entity_type": "drug",
                        "original_term": raw,
                        "cleaned_term": debug.get("normalized", ""),
                        "mapped_term": mapped or "",
                        "mesh_id": mesh_id_for(mapped, mesh_ids) if mapped else "",
                        "mapping_method": method,
                        "mapping_score": confidence,
                        "mapping_status": "mapped" if mapped else "unmapped",
                    },
                )
                if mapped:
                    mapped_drugs.append(
                        (
                            raw,
                            builder.clean_output_name(mapped),
                            mesh_id_for(mapped, mesh_ids),
                            debug,
                            method,
                            confidence,
                        )
                    )
                else:
                    unmatched_entities[("drug", raw)] = {
                        "entity_type": "drug",
                        "raw_term": raw,
                        "reason": "unmatched intervention",
                        "debug": debug,
                    }

            pair_attempts += len(raw_conditions) * len(raw_drugs)
            phases = trial.get("phases", [])
            if not isinstance(phases, list):
                phases = [phases] if phases else []
            for (
                raw_condition,
                condition,
                condition_mesh_id,
                condition_debug,
                condition_method,
                condition_confidence,
            ) in mapped_conditions:
                for (
                    raw_drug,
                    drug,
                    drug_mesh_id,
                    drug_debug,
                    drug_method,
                    drug_confidence,
                ) in mapped_drugs:
                    matched_rows.append(
                        {
                            "condition": condition,
                            "intervention": drug,
                            "raw_condition": raw_condition,
                            "raw_intervention": raw_drug,
                            "condition_cleaned_term": condition_debug.get("normalized", ""),
                            "intervention_cleaned_term": drug_debug.get("normalized", ""),
                            "condition_mesh_id": condition_mesh_id,
                            "intervention_mesh_id": drug_mesh_id,
                            "condition_mapping_method": condition_method,
                            "condition_mapping_confidence": condition_confidence,
                            "condition_mapping_score": debug_score(condition_debug, condition_confidence),
                            "condition_token_jaccard_score": debug_token_jaccard(condition_debug),
                            "intervention_mapping_method": drug_method,
                            "intervention_mapping_confidence": drug_confidence,
                            "intervention_mapping_score": debug_score(drug_debug, drug_confidence),
                            "intervention_token_jaccard_score": debug_token_jaccard(drug_debug),
                            "condition_debug": json.dumps(condition_debug, ensure_ascii=False),
                            "drug_debug": json.dumps(drug_debug, ensure_ascii=False),
                            "nct_id": trial.get("nctId", ""),
                            "phases": phases,
                            "status": trial.get("status", ""),
                            "source_file": path.name,
                        }
                    )
        print(
            f"MeSH mapping chunk {chunk_index}/{len(chunk_paths)}: "
            f"{trials_processed:,} trials, {len(matched_rows):,} mapped rows",
            flush=True,
        )

    matched_path.write_text(
        json.dumps(matched_rows, ensure_ascii=False), encoding="utf-8"
    )
    unmatched_path.write_text(
        json.dumps(list(unmatched_entities.values()), ensure_ascii=False),
        encoding="utf-8",
    )
    audit_df = pd.DataFrame(entity_audit.values())
    summary_rows = pd.DataFrame(
        [
            {
                "entity_type": "summary",
                "original_term": "FULL_REGISTRY",
                "cleaned_term": "",
                "mapped_term": "",
                "mesh_id": "",
                "mapping_method": "exact_entry_term" if args.mesh_exact_only else "legacy_multistage",
                "mapping_score": "",
                "mapping_status": (
                    f"trials={trials_processed};pair_attempts={pair_attempts};"
                    f"mapped_rows={len(matched_rows)};placebo_excluded={placebo_excluded}"
                ),
            }
        ]
    )
    audit_df = pd.concat([audit_df, summary_rows], ignore_index=True)
    log_lines.append(
        f"Full-registry MeSH mapping: {trials_processed} trials, {pair_attempts} raw pair attempts, "
        f"{len(matched_rows)} mapped rows, {len(unmatched_entities)} unique unmapped entities"
    )
    return audit_df, matched_path, unmatched_path


def graph_stage(
    matched_path: Path,
    graph_dir: Path,
    target_pairs: List[Tuple[str, str]],
    args: argparse.Namespace,
    log_lines: List[str],
) -> Tuple[pd.DataFrame, Path, Path, Path]:
    """Build the canonical leakage-controlled graph for every pipeline mode."""
    graph_dir.mkdir(parents=True, exist_ok=True)
    audit, known_path, unknown_path, weights_path = build_targeted_full_graph(
        matched_path=matched_path,
        graph_dir=graph_dir,
        target_pairs=target_pairs,
        random_seed=args.random_seed,
        max_training_positives=args.max_graph_training_positives,
    )
    log_lines.append(
        f"Saved leakage-controlled graph features to {rel(graph_dir)}; "
        f"target edges withheld and weights written to {rel(weights_path)}"
    )
    return audit, known_path, unknown_path, weights_path


def literature_retrieval_audit(literature_dir: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for path in sorted(literature_dir.glob("classified_*.json")) if literature_dir.exists() else []:
        records = read_json(path, [])
        if not isinstance(records, list):
            records = []
        rows.append(
            {
                "literature_file": rel(path),
                "timestamp": timestamp_from_path(path),
                "records_retrieved_or_saved": len(records),
                "records_usable": len(records),
                "records_with_abstracts": sum(1 for r in records if str(r.get("abstract", "")).strip()),
                "records_with_intro": sum(1 for r in records if str(r.get("introduction", "")).strip()),
                "records_with_conclusion": sum(1 for r in records if str(r.get("conclusion", "")).strip()),
                "records_excluded_missing_text": 0,
            }
        )
    return pd.DataFrame(rows)


def semantic_classification_audit(literature_dir: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for path in sorted(literature_dir.glob("classified_*.json")) if literature_dir.exists() else []:
        records = read_json(path, [])
        if not isinstance(records, list):
            records = []
        counts = Counter(norm_entity(row.get("category")) for row in records)
        total = sum(counts.values())
        benefit = counts.get("benefit", 0)
        null_count = counts.get("null", 0)
        harm = counts.get("harm", 0)
        conflicting = counts.get("conflicting", 0)
        therapeutic = benefit + counts.get("therapeutic", 0)
        adverse = null_count + harm + conflicting + counts.get("adverse", 0)
        rows.append(
            {
                "literature_file": rel(path),
                "classified_total": total,
                "benefit": benefit,
                "null": null_count,
                "harm": harm,
                "conflicting": conflicting,
                "therapeutic": therapeutic,
                "adverse": adverse,
                "irrelevant": counts.get("irrelevant", 0),
                "irrelevant_noise_rate": counts.get("irrelevant", 0) / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def load_run_records(runs_dir: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for path in sorted(runs_dir.glob("run_*.json")) if runs_dir.exists() else []:
        payload = read_json(path, {})
        if not isinstance(payload, dict):
            continue
        components = payload.get("components") or {}
        if not isinstance(components, dict):
            components = {}
        rows.append(
            {
                "run_log": rel(path),
                "timestamp": payload.get("timestamp", timestamp_from_path(path)),
                "drug": payload.get("drug", ""),
                "disease": payload.get("disease", ""),
                "components": components,
            }
        )
    return rows


def safety_overlap_audit(runs_dir: Path) -> pd.DataFrame:
    rows = []
    for record in load_run_records(runs_dir):
        components = record["components"]
        matching_effects = components.get("matching_effects")
        if not isinstance(matching_effects, list):
            matching_effects = []
        rows.append(
            {
                "run_log": record["run_log"],
                "drug": record["drug"],
                "disease": record["disease"],
                "gamma_safety_overlap": components.get("gamma"),
                "safety_data_status": components.get("safety_data_status", "UNKNOWN"),
                "literature_data_status": components.get("literature_data_status", "UNKNOWN"),
                "safety_overlap_term_count": len(matching_effects) if matching_effects else "",
                "top_safety_overlap_terms": "; ".join(map(str, matching_effects[:10])),
                "p_penalised": components.get("p_penalised"),
                "p_final": components.get("p_final"),
            }
        )
    return pd.DataFrame(rows)


def bayesian_uncertainty_audit(runs_dir: Path) -> pd.DataFrame:
    rows = []
    for record in load_run_records(runs_dir):
        components = record["components"]
        uncertainty = beta_metrics(components)
        rows.append(
            {
                "run_log": record["run_log"],
                "drug": record["drug"],
                "disease": record["disease"],
                "p_raw": components.get("p_raw"),
                "p_penalised": components.get("p_penalised"),
                "p_final": components.get("p_final"),
                "M": components.get("M"),
                "cM": components.get("cM"),
                "p_likelihood": components.get("p_likelihood"),
                "posterior_mean": components.get("post_mean"),
                "posterior_variance": round_or_none(uncertainty.get("posterior_variance"), 8),
                "posterior_ci_low": round_or_none(uncertainty.get("posterior_ci_low")),
                "posterior_ci_high": round_or_none(uncertainty.get("posterior_ci_high")),
                "credible_interval_width": round_or_none(uncertainty.get("credible_interval_width")),
                "kl_divergence": round_or_none(uncertainty.get("kl_divergence")),
                "mean_shift": round_or_none(uncertainty.get("mean_shift")),
            }
        )
    return pd.DataFrame(rows)


def safe_run_slug(text: str) -> str:
    return re.sub(r"[^a-z0-9_,.-]+", "_", text.strip().lower()).strip("_")


def combined_graph_df(known_graph_path: Path, unknown_graph_path: Path) -> pd.DataFrame:
    frames = []
    for path in (known_graph_path, unknown_graph_path):
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["Drug"] = df["Drug"].astype(str).str.strip().str.lower()
    df["Disease"] = df["Disease"].astype(str).str.strip().str.lower()
    return df


def load_existing_classified_prior(
    literature_dir: Path,
    runs_dir: Path,
    drug: str,
    disease: str,
) -> Dict[str, Any]:
    """Rebuild prior inputs from saved classifications without calling an LLM."""
    safe_drug = safe_run_slug(drug)
    safe_disease = safe_run_slug(disease)
    candidates = sorted(
        literature_dir.glob(f"classified_pubmed_{safe_drug}_{safe_disease}_*.json")
    )
    if not candidates:
        raise FileNotFoundError(
            f"No saved classified literature found for {drug} / {disease} in {literature_dir}"
        )
    records = read_json(candidates[-1], [])
    if not isinstance(records, list) or not records:
        raise ValueError(f"Saved classified literature is empty: {candidates[-1]}")

    aliases = {
        "therapeutic": "benefit",
        "adverse": "conflicting",
        "no_benefit": "null",
        "safety": "harm",
    }
    counts = Counter()
    evidence_keys = set()
    for record in records:
        category = aliases.get(str(record.get("category", "")).strip().lower(), str(record.get("category", "")).strip().lower())
        if category not in {"benefit", "null", "harm", "conflicting", "irrelevant"}:
            category = "irrelevant"
        counts[category] += 1
        key = str(record.get("doi") or record.get("pmid") or record.get("Title") or "").strip().casefold()
        if key:
            evidence_keys.add(key)

    total = sum(counts.values())
    benefit = counts.get("benefit", 0)
    adverse = counts.get("null", 0) + counts.get("harm", 0) + counts.get("conflicting", 0)
    counts["therapeutic"] = benefit
    counts["adverse"] = adverse

    safety = {}
    run_candidates = sorted(runs_dir.glob(f"run_{safe_drug}_{safe_disease}_*.json"))
    if run_candidates:
        payload = read_json(run_candidates[-1], {})
        safety = payload.get("components", {}) if isinstance(payload, dict) else {}
    gamma = safety.get("gamma")
    gamma = clamp01(gamma) if gamma is not None else None
    relation = safety.get("safety_relation")
    penalty_scale = float(safety.get("safety_penalty_scale") or 0.5)
    penalised = max((benefit - 2 * adverse) / total, 0.0) if total else 0.5
    p_final = penalised * (1.0 - penalty_scale * gamma) if relation and gamma is not None else penalised

    return {
        "prior": benefit / total if total else 0.5,
        "penalised_prior": penalised,
        "enhanced_prior": clamp01(p_final),
        "gamma": gamma,
        "safety_data_status": safety.get("safety_data_status", "LEGACY_IMPORTED"),
        "literature_data_status": "LEGACY_IMPORTED",
        "literature_query": safety.get("literature_query", ""),
        "classifier_schema_version": "legacy_imported",
        "raw_counts": counts,
        "total_articles": total,
        "evidence_units": len(evidence_keys) if evidence_keys else total,
        "records_retrieved": safety.get("records_retrieved", total),
        "records_excluded_by_exact_verification": safety.get("records_excluded_by_exact_verification", 0),
        "matching_effects": safety.get("matching_effects", []),
        "safety_relation": relation,
        "safety_penalty_scale": penalty_scale,
    }


def run_bayesian_panel_stage(
    panel_df: pd.DataFrame,
    known_graph_path: Path,
    unknown_graph_path: Path,
    output_dir: Path,
    args: argparse.Namespace,
    log_lines: List[str],
) -> Tuple[pd.DataFrame, Path, Path]:
    """Refresh PubMed/PMC literature, safety overlap, and Bayesian scoring for panel pairs."""
    from bayesian_predictor import (  # noqa: PLC0415
        LLMConfig,
        PubMedSearchConfig,
        concentration_c,
        make_classifier,
        PredictorConfig,
        beta_params_from_prob,
        clamp01,
        sigmoid,
    )

    runs_dir = output_dir / "runs"
    literature_dir = output_dir / "literatures"
    source_literature_dir = literature_dir
    source_runs_dir = runs_dir
    if args.reuse_classified_literature and str(args.input_run_dir or "").strip():
        input_root = (PROJECT_ROOT / args.input_run_dir).resolve()
        source_literature_dir = input_root / "literatures"
        source_runs_dir = input_root / "runs"
    plots_dir = output_dir / "plots"
    for path in (runs_dir, literature_dir, plots_dir):
        path.mkdir(parents=True, exist_ok=True)

    cfg = PredictorConfig(
        pubmed_max_articles=args.pubmed_max_articles,
        pubmed_filter_level=args.pubmed_filter_level,
        pubmed_years_back=args.pubmed_years_back,
        pubmed_use_cache=not args.force_pubmed_refresh,
        llm_model=args.llm_model,
        llm_batch_size=args.llm_batch_size,
        llm_delay_s=args.llm_delay_s,
        llm_max_retries=args.llm_max_retries,
        cmax=args.cmax,
        tau=args.tau,
        likelihood_strength=args.likelihood_strength,
        likelihood_intercept=args.likelihood_intercept,
        plots_dir=str(plots_dir),
        logs_dir=str(runs_dir),
        save_run_logs=True,
    )
    search_cfg = PubMedSearchConfig(
        max_results=args.pubmed_max_articles,
        years_back=args.pubmed_years_back,
        enable_pmc_augmentation=args.pubmed_enable_pmc,
    )
    llm_cfg = LLMConfig(
        model=args.llm_model,
        delay_s=args.llm_delay_s,
        batch_size=args.llm_batch_size,
        max_retries=args.llm_max_retries,
    )
    # Use make_classifier to preserve current predictor behavior, then override the
    # save_dir per call so publication artifacts stay inside the run folder.
    classifier = None
    if not args.reuse_classified_literature:
        classifier = make_classifier(cfg)
        classifier.search_cfg = search_cfg
        classifier.pubmed.cfg = search_cfg
        classifier.llm_cfg = llm_cfg

    feature_df = combined_graph_df(known_graph_path, unknown_graph_path)
    active_weights = dict(DEFAULT_WEIGHTS)
    active_intercept = float(args.likelihood_intercept)
    weights_source = "legacy_default_weights"
    if str(args.graph_weights_path or "").strip():
        graph_weights_path = Path(args.graph_weights_path)
        if not graph_weights_path.is_absolute():
            graph_weights_path = (PROJECT_ROOT / graph_weights_path).resolve()
        payload = read_json(graph_weights_path, {})
        candidate_weights = payload.get("raw_feature_weights", {}) if isinstance(payload, dict) else {}
        if candidate_weights:
            active_weights = {
                feature: float(candidate_weights[feature])
                for feature in DEFAULT_WEIGHTS
                if feature in candidate_weights
            }
            active_intercept = float(
                payload.get("raw_feature_intercept", active_intercept)
            )
            weights_source = rel(graph_weights_path)
    log_lines.append(
        f"Bayesian graph weights source: {weights_source}; intercept={active_intercept}"
    )
    rows: List[Dict[str, Any]] = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    for _, panel_row in panel_df.iterrows():
        drug = str(panel_row["drug"]).strip()
        disease = str(panel_row["disease"]).strip()
        drug_l = drug.lower()
        disease_l = disease.lower()
        log_lines.append(f"Refreshing Bayesian evidence for panel pair: {drug} -> {disease}")

        if args.reuse_classified_literature:
            prior_result = load_existing_classified_prior(
                source_literature_dir, source_runs_dir, drug, disease
            )
            # Carry the exact classified record file into the new run so the
            # literature audit stages remain populated and auditable.  This is
            # a byte-for-byte provenance copy, not a new retrieval or label.
            imported_files = sorted(
                source_literature_dir.glob(
                    f"classified_pubmed_{safe_run_slug(drug)}_{safe_run_slug(disease)}_*.json"
                )
            )
            if imported_files:
                destination = literature_dir / imported_files[-1].name
                if imported_files[-1].resolve() != destination.resolve():
                    shutil.copyfile(imported_files[-1], destination)
                log_lines.append(
                    f"Imported saved classified literature file for audit: {rel(destination)}"
                )
            log_lines.append("Reused saved literature classifications; no LLM classification was run.")
        else:
            prior_result = classifier.build_semantic_prior(
                drug=drug,
                disease=disease,
                max_articles=args.pubmed_max_articles,
                filter_level=args.pubmed_filter_level,
                save_dir=str(literature_dir),
                use_cache=not args.force_pubmed_refresh,
            )

        counts = dict(prior_result.get("raw_counts", {}))
        articles_retrieved = int(prior_result.get("total_articles", 0))
        m_articles = int(prior_result.get("evidence_units", articles_retrieved))
        p_raw = clamp01(prior_result.get("prior", 0.5))
        p_pen = clamp01(prior_result.get("penalised_prior", 0.5))
        p_final = clamp01(prior_result.get("enhanced_prior", 0.5))
        gamma = prior_result.get("gamma", None)
        gamma = clamp01(gamma) if gamma is not None else None
        c_m = concentration_c(m_articles, cmax=args.cmax, tau=args.tau)
        prior_a, prior_b = beta_params_from_prob(p_final, c_m)

        graph_score = active_intercept
        feature_values: Dict[str, float] = {}
        row_match = pd.DataFrame()
        if not feature_df.empty:
            row_match = feature_df[(feature_df["Drug"] == drug_l) & (feature_df["Disease"] == disease_l)]
        if not row_match.empty:
            graph_row = row_match.iloc[0]
            for feature, weight in active_weights.items():
                value = float(graph_row.get(feature, 0.0))
                feature_values[feature] = value
                graph_score += float(weight) * value

        p_likelihood = clamp01(sigmoid(graph_score))
        like_a, like_b = beta_params_from_prob(p_likelihood, args.likelihood_strength)
        post_a = prior_a + (like_a - 1.0)
        post_b = prior_b + (like_b - 1.0)
        post_mean = post_a / (post_a + post_b)

        components = {
            "p_raw": p_raw,
            "p_penalised": p_pen,
            "p_final": p_final,
            "gamma": gamma,
            "safety_data_status": prior_result.get("safety_data_status", "UNKNOWN"),
            "literature_data_status": prior_result.get("literature_data_status", "UNKNOWN"),
            "literature_query": prior_result.get("literature_query", ""),
            "classifier_schema_version": prior_result.get("classifier_schema_version", "unknown"),
            "counts": counts,
            "M": m_articles,
            "articles_retrieved": articles_retrieved,
            "evidence_units": m_articles,
            "records_retrieved": prior_result.get("records_retrieved", m_articles),
            "records_excluded_by_exact_verification": prior_result.get(
                "records_excluded_by_exact_verification", 0
            ),
            "cM": c_m,
            "prior_a": prior_a,
            "prior_b": prior_b,
            "graph_score": graph_score,
            "p_likelihood": p_likelihood,
            "likelihood_strength": float(args.likelihood_strength),
            "like_a": like_a,
            "like_b": like_b,
            "feature_values": feature_values,
            "post_a": post_a,
            "post_b": post_b,
            "post_mean": float(post_mean),
            "matching_effects": prior_result.get("matching_effects", []),
            "safety_relation": prior_result.get("safety_relation"),
            "safety_penalty_scale": prior_result.get("safety_penalty_scale"),
        }

        payload = {
            "timestamp": timestamp,
            "drug": drug,
            "disease": disease,
            "case_type": panel_row.get("case_type", ""),
            "config": {
                "pubmed_max_articles": args.pubmed_max_articles,
                "pubmed_filter_level": args.pubmed_filter_level,
                "pubmed_years_back": args.pubmed_years_back,
                "pubmed_use_cache": not args.force_pubmed_refresh,
                "llm_model": args.llm_model,
                "llm_batch_size": args.llm_batch_size,
                "llm_delay_s": args.llm_delay_s,
                "llm_max_retries": args.llm_max_retries,
                "cmax": args.cmax,
                "tau": args.tau,
                "likelihood_strength": args.likelihood_strength,
                "likelihood_intercept": active_intercept,
                "weights": active_weights,
                "weights_source": weights_source,
            },
            "components": components,
        }
        run_path = runs_dir / f"run_{safe_run_slug(drug)}_{safe_run_slug(disease)}_{timestamp}.json"
        run_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        rows.append(
            {
                "drug": drug,
                "disease": disease,
                "run_log": rel(run_path),
                "articles_used": m_articles,
                "gamma": gamma,
                "safety_data_status": prior_result.get("safety_data_status", "UNKNOWN"),
                "literature_data_status": prior_result.get("literature_data_status", "UNKNOWN"),
                "posterior_mean": post_mean,
                "status": "completed",
            }
        )

    log_lines.append(f"Refreshed Bayesian panel for {len(rows)} drug-disease pairs.")
    return pd.DataFrame(rows), runs_dir, literature_dir


def pair_norm(drug: Any, disease: Any) -> Tuple[str, str]:
    return str(drug or "").strip().lower(), str(disease or "").strip().lower()


def _to_float(value: Any) -> Optional[float]:
    """Convert a value to float, returning None if conversion fails."""
    if value is None:
        return None
    try:
        f = float(value)
        return None if np.isnan(f) else f
    except (TypeError, ValueError):
        return None


def build_trial_pair_summary(
    terminology_audit_path: Path,
    matched_path: Optional[Path] = None,
) -> pd.DataFrame:
    df = pd.read_csv(terminology_audit_path) if terminology_audit_path.exists() else pd.DataFrame()
    matched = (
        df[df["mapping_status"] == "matched"].copy()
        if not df.empty and "mapping_status" in df.columns
        else pd.DataFrame()
    )
    # Full-registry mapping stores one provenance row per unique entity
    # (mapping_status=mapped) and the complete pair provenance in the matched
    # JSON.  When the audit CSV is a refresh=false stub, recover that pair
    # table instead of silently reporting zero trials.
    if matched.empty and matched_path is not None and matched_path.exists():
        raw_rows = read_json(matched_path, [])
        if isinstance(raw_rows, list) and raw_rows:
            matched = pd.DataFrame(raw_rows)
            if "mapped_drug" not in matched.columns:
                matched["mapped_drug"] = matched.get("intervention", "")
            if "mapped_disease" not in matched.columns:
                matched["mapped_disease"] = matched.get("condition", "")
            if "drug_mapping_method" not in matched.columns:
                matched["drug_mapping_method"] = matched.get(
                    "intervention_mapping_method", ""
                )
            if "disease_mapping_method" not in matched.columns:
                matched["disease_mapping_method"] = matched.get(
                    "condition_mapping_method", ""
                )
            matched["mapping_status"] = "matched"
    if matched.empty:
        return pd.DataFrame()
    matched["drug_key"] = matched["mapped_drug"].astype(str).str.lower().str.strip()
    matched["disease_key"] = matched["mapped_disease"].astype(str).str.lower().str.strip()

    def first_nonempty(series: Any) -> str:
        if series is None:
            return ""
        for value in series:
            text = str(value or "").strip()
            if text and text.lower() != "nan":
                return text
        return ""

    def dist(series: pd.Series) -> str:
        counts = Counter()
        for value in series.dropna():
            for item in str(value).split("|"):
                item = item.strip()
                if item:
                    counts[item] += 1
        return json.dumps(dict(sorted(counts.items())), sort_keys=True)

    grouped = matched.groupby(["drug_key", "disease_key"], dropna=False)
    rows = []
    for (drug, disease), sub in grouped:
        rows.append(
            {
                "drug_key": drug,
                "disease_key": disease,
                "trial_count": sub["nct_id"].nunique() if "nct_id" in sub else len(sub),
                "phase_distribution": dist(sub["phases"]) if "phases" in sub else "{}",
                "status_distribution": dist(sub["status"]) if "status" in sub else "{}",
                "drug_mapping_method": sub["drug_mapping_method"].mode().iloc[0] if "drug_mapping_method" in sub and not sub["drug_mapping_method"].mode().empty else "",
                "disease_mapping_method": sub["disease_mapping_method"].mode().iloc[0] if "disease_mapping_method" in sub and not sub["disease_mapping_method"].mode().empty else "",
                "drug_mapping_score": pd.to_numeric(sub.get("drug_mapping_score"), errors="coerce").mean() if "drug_mapping_score" in sub else np.nan,
                "disease_mapping_score": pd.to_numeric(sub.get("disease_mapping_score"), errors="coerce").mean() if "disease_mapping_score" in sub else np.nan,
                "drug_mesh_id": first_nonempty(sub.get("intervention_mesh_id"))
                if "intervention_mesh_id" in sub
                else first_nonempty(sub.get("drug_mesh_id")),
                "disease_mesh_id": first_nonempty(sub.get("condition_mesh_id"))
                if "condition_mesh_id" in sub
                else first_nonempty(sub.get("disease_mesh_id")),
            }
        )
    return pd.DataFrame(rows)


def build_publication_ledger(
    pair_level_path: Path,
    terminology_audit_path: Path,
    validation_path: Optional[Path],
    panel_df: pd.DataFrame,
    matched_path: Optional[Path] = None,
) -> pd.DataFrame:
    pair_df = pd.read_csv(pair_level_path) if pair_level_path.exists() else pd.DataFrame()
    # A reuse run deliberately keeps the mapping audit as provenance, but the
    # audit can be a small refresh=false note.  Pass the actual matched pair
    # table so trial counts and mapping identifiers are recovered rather than
    # silently reported as zero/legacy.
    processed_path = matched_path or (terminology_audit_path.parent.parent / "processed_data" / "condition_drug_pairs.json")
    trial_df = build_trial_pair_summary(terminology_audit_path, processed_path)
    validation_df = pd.read_csv(validation_path) if validation_path and validation_path.exists() else pd.DataFrame()

    trial_lookup = {
        (row["drug_key"], row["disease_key"]): row.to_dict()
        for _, row in trial_df.iterrows()
    } if not trial_df.empty else {}
    validation_lookup = {
        pair_norm(row.get("drug"), row.get("disease")): row.to_dict()
        for _, row in validation_df.iterrows()
    } if not validation_df.empty else {}

    rows: List[Dict[str, Any]] = []
    for _, row in pair_df.iterrows():
        key = pair_norm(row.get("drug"), row.get("disease"))
        trial = trial_lookup.get(key, {})
        validation = validation_lookup.get(key, {})
        has_mapping = bool(trial)
        has_graph = pd.notna(row.get("GraphDistanceToIndication")) or pd.notna(row.get("RandomWalkScore"))
        has_lit = float(row.get("records_retrieved", 0) or 0) > 0
        safety_status = str(row.get("safety_data_status", "UNKNOWN") or "UNKNOWN").upper()
        has_safety = pd.notna(row.get("gamma_safety_overlap")) and safety_status == "COMPLETE"
        has_bayes = pd.notna(row.get("posterior_mean"))
        tier = validation.get("coverage_tier") or rule_coverage_tier(
            has_mapping=has_mapping,
            has_graph=has_graph,
            has_literature=has_lit,
            has_safety=has_safety,
            has_bayesian=has_bayes,
        )
        # ── derived / computed columns ─────────────────────────────────────
        t_count = _to_float(row.get("therapeutic_count", 0))
        a_count = _to_float(row.get("adverse_count", 0))
        irr_count = _to_float(row.get("irrelevant_count", 0))
        total_classified = t_count + a_count + irr_count
        therapeutic_ratio = round(t_count / total_classified, 4) if total_classified > 0 else None
        adverse_burden = round(a_count / total_classified, 4) if total_classified > 0 else None
        articles_ret = _to_float(row.get("records_retrieved", 0))
        usable = _to_float(row.get("records_usable", 0))
        lit_completeness = None
        if articles_ret > 0:
            # Exhaustive mode has no arbitrary article-count ceiling.  This
            # component measures retrieval usability, not an article-count
            # threshold.
            lit_completeness = round(usable / articles_ret, 4)

        post_mean_val = _to_float(row.get("posterior_mean"))
        prior_mean_val = _to_float(row.get("p_final"))
        mean_shift = round(post_mean_val - prior_mean_val, 4) if post_mean_val is not None and prior_mean_val is not None else None

        ci_w = _to_float(row.get("credible_interval_width"))
        if ci_w is None:
            uncertainty_level = None
        elif ci_w < 0.15:
            uncertainty_level = "low"
        elif ci_w <= 0.35:
            uncertainty_level = "moderate"
        else:
            uncertainty_level = "high"

        # Structural consistency: average of available graph scores normalised to [0,1]
        graph_scores = [_to_float(row.get(c)) for c in ("RandomWalkScore", "KatzSimilarity", "StructuralLikelihood") if _to_float(row.get(c)) is not None]
        structural_consistency = round(float(np.mean(graph_scores)), 4) if graph_scores else None

        # Mapping method taxonomy — elevate legacy entries if provenance is missing
        def _clean_method(m: Any) -> str:
            m = str(m or "").strip()
            return m if m else "mapped_legacy_no_provenance"

        drug_method = _clean_method(trial.get("drug_mapping_method", row.get("drug_mapping_method", "")))
        disease_method = _clean_method(trial.get("disease_mapping_method", row.get("disease_mapping_method", "")))

        drug_score = round_or_none(trial.get("drug_mapping_score", row.get("entity_mapping_quality_score")))
        disease_score = round_or_none(trial.get("disease_mapping_score", row.get("entity_mapping_quality_score")))

        # Derive mapping status from method
        def _method_to_status(method: str) -> str:
            if method in ("exact_match", "fuzzy_high_confidence", "fuzzy_low_confidence"):
                return "mapped"
            if method == "unmapped":
                return "unmapped"
            return "legacy"

        # Safety penalty approximation: 0.5 * gamma
        gamma_val = _to_float(row.get("gamma_safety_overlap"))
        safety_penalty = round(0.5 * gamma_val, 4) if gamma_val is not None else None

        # Final interpretation label
        qflag = str(row.get("quality_flag", "") or "")
        if "High evidence" in qflag:
            interpretation = "Evidence-ready for further clinical review."
        elif "Safety-conflicted" in qflag or "Safety-concerning" in qflag:
            interpretation = "Safety signal requires resolution before prioritisation."
        elif "Insufficient" in qflag:
            interpretation = "Evidence base too sparse for reliable assessment."
        elif "noise" in qflag.lower():
            interpretation = "Literature retrieval dominated by irrelevant results."
        elif "conflicted" in qflag.lower():
            interpretation = "Conflicting evidence signals require expert adjudication."
        elif "Terminology" in qflag:
            interpretation = "Terminology mapping uncertainty reduces confidence."
        elif "Moderate" in qflag:
            interpretation = "Moderate evidence; further validation recommended."
        else:
            interpretation = "Insufficient data for interpretation."

        rows.append(
            {
                # Provenance columns
                "drug_original": row.get("drug"),
                "disease_original": row.get("disease"),
                "drug_cleaned": row.get("drug_cleaned", row.get("drug")),
                "disease_cleaned": row.get("disease_cleaned", row.get("disease")),
                "drug_mapped": trial.get("drug_key", row.get("drug")),
                "disease_mapped": trial.get("disease_key", row.get("disease")),
                "drug_mesh_id": trial.get("drug_mesh_id", row.get("drug_mesh_id", "")),
                "disease_mesh_id": trial.get("disease_mesh_id", row.get("disease_mesh_id", "")),
                # Short-form aliases
                "drug": row.get("drug"),
                "disease": row.get("disease"),
                # Mapping audit
                "drug_mapping_method": drug_method,
                "disease_mapping_method": disease_method,
                "drug_mapping_score": drug_score,
                "disease_mapping_score": disease_score,
                "drug_token_jaccard_score": round_or_none(trial.get("drug_token_jaccard_score", row.get("drug_token_jaccard_score"))),
                "disease_token_jaccard_score": round_or_none(trial.get("disease_token_jaccard_score", row.get("disease_token_jaccard_score"))),
                "drug_mapping_status": _method_to_status(drug_method),
                "disease_mapping_status": _method_to_status(disease_method),
                # Trial evidence
                "trial_count": int(trial.get("trial_count", 0) or 0),
                "phase_distribution": trial.get("phase_distribution", "{}"),
                "status_distribution": trial.get("status_distribution", "{}"),
                # Literature
                "articles_retrieved": row.get("records_retrieved"),
                "usable_articles": row.get("records_usable"),
                "therapeutic_count": row.get("therapeutic_count"),
                "adverse_count": row.get("adverse_count"),
                "irrelevant_count": row.get("irrelevant_count"),
                "therapeutic_ratio": therapeutic_ratio,
                "adverse_burden": adverse_burden,
                "irrelevant_noise_rate": row.get("irrelevant_retrieval_noise_rate"),
                "literature_completeness_score": lit_completeness,
                "literature_data_status": row.get("literature_data_status", "UNKNOWN"),
                # Safety
                "safety_overlap_gamma": row.get("gamma_safety_overlap"),
                "safety_data_status": safety_status,
                "safety_penalty": safety_penalty,
                # Graph
                "graph_distance": row.get("GraphDistanceToIndication"),
                "random_walk_score": row.get("RandomWalkScore"),
                "katz_similarity": row.get("KatzSimilarity"),
                "preferential_attachment": row.get("PreferentialAttachment"),
                "structural_likelihood": row.get("StructuralLikelihood"),
                "structural_consistency_score": structural_consistency,
                # Bayesian
                "prior_mean": row.get("p_final"),
                "posterior_mean": row.get("posterior_mean"),
                "credible_interval_lower": row.get("posterior_ci_low"),
                "credible_interval_upper": row.get("posterior_ci_high"),
                "credible_interval_width": row.get("credible_interval_width"),
                "kl_divergence": row.get("kl_divergence"),
                "mean_shift": mean_shift,
                # Composite
                "evidence_readiness_score": row.get("evidence_readiness_score"),
                "uncertainty_level": uncertainty_level,
                "quality_flag": row.get("quality_flag"),
                "coverage_tier": tier,
                "final_interpretation": interpretation,
            }
        )

    ledger = pd.DataFrame(rows)
    if ledger.empty:
        ledger = pd.DataFrame(columns=LEDGER_COLUMNS)
    for col in LEDGER_COLUMNS:
        if col not in ledger.columns:
            ledger[col] = ""
    ledger = ledger[LEDGER_COLUMNS]

    if not panel_df.empty:
        panel_keys = {pair_norm(row["drug"], row["disease"]) for _, row in panel_df.iterrows()}
        ledger["_panel_order"] = ledger.apply(
            lambda row: 0 if pair_norm(row["drug"], row["disease"]) in panel_keys else 1,
            axis=1,
        )
        ledger = ledger.sort_values(["_panel_order", "drug", "disease"]).drop(columns=["_panel_order"])
    return ledger


def validate_publication_outputs(
    *,
    ledger: pd.DataFrame,
    terminology_audit_path: Path,
    panel_df: pd.DataFrame,
    output_dir: Path,
    strict: bool,
    log_lines: List[str],
) -> pd.DataFrame:
    checks: List[Dict[str, Any]] = []

    def add_check(name: str, status: str, detail: str, severity: str = "warning") -> None:
        checks.append({"check": name, "status": status, "severity": severity, "detail": detail})

    terminology = pd.read_csv(terminology_audit_path) if terminology_audit_path.exists() else pd.DataFrame()
    if terminology.empty:
        add_check("mapping_provenance_present", "fail", "Terminology audit table is empty.", "error")
    else:
        required_mapping_cols = {
            "original_term",
            "cleaned_term",
            "mapped_term",
            "mesh_id",
            "mapping_method",
            "mapping_score",
            "token_jaccard_score",
            "mapping_status",
            "failure_reason",
        }
        missing = required_mapping_cols - set(terminology.columns)
        add_check(
            "mapping_provenance_present",
            "fail" if missing else "pass",
            f"Missing columns: {sorted(missing)}" if missing else "All provenance columns present.",
            "error" if missing else "info",
        )

    duplicate_count = 0
    if not ledger.empty:
        duplicate_count = int(ledger.duplicated(subset=["drug", "disease"], keep=False).sum())
    add_check(
        "duplicate_drug_disease_pairs",
        "fail" if duplicate_count else "pass",
        f"Duplicate ledger rows: {duplicate_count}",
        "error" if duplicate_count else "info",
    )

    panel_missing = []
    panel_no_lit = []
    panel_no_safety = []
    panel_no_uncertainty = []
    ledger_index = {pair_norm(row["drug"], row["disease"]): row for _, row in ledger.iterrows()}
    for _, panel_row in panel_df.iterrows():
        key = pair_norm(panel_row["drug"], panel_row["disease"])
        row = ledger_index.get(key)
        if row is None:
            panel_missing.append(f"{panel_row['drug']} -> {panel_row['disease']}")
            continue
        if pd.isna(row.get("posterior_mean")) or row.get("posterior_mean") == "":
            panel_missing.append(f"{panel_row['drug']} -> {panel_row['disease']}")
        articles_value = pd.to_numeric(pd.Series([row.get("articles_retrieved")]), errors="coerce").fillna(0).iloc[0]
        if float(articles_value) == 0:
            panel_no_lit.append(f"{panel_row['drug']} -> {panel_row['disease']}")
        if pd.isna(row.get("safety_overlap_gamma")) or row.get("safety_overlap_gamma") == "":
            panel_no_safety.append(f"{panel_row['drug']} -> {panel_row['disease']}")
        if pd.isna(row.get("credible_interval_width")) or row.get("credible_interval_width") == "":
            panel_no_uncertainty.append(f"{panel_row['drug']} -> {panel_row['disease']}")

    add_check("bayesian_scores_for_panel", "fail" if panel_missing else "pass", "; ".join(panel_missing[:20]) or "All panel pairs have Bayesian rows.", "error" if panel_missing else "info")
    add_check("literature_nonzero_for_panel", "fail" if panel_no_lit else "pass", "; ".join(panel_no_lit[:20]) or "All panel pairs have nonzero literature counts.", "error" if panel_no_lit else "info")
    add_check("safety_gamma_for_panel", "fail" if panel_no_safety else "pass", "; ".join(panel_no_safety[:20]) or "All panel pairs have safety gamma.", "warning" if panel_no_safety else "info")
    add_check("posterior_uncertainty_for_panel", "fail" if panel_no_uncertainty else "pass", "; ".join(panel_no_uncertainty[:20]) or "All panel pairs have posterior uncertainty.", "error" if panel_no_uncertainty else "info")

    expected_rows = len(ledger)
    quality_counts_path = output_dir / "09_quality_category_counts.csv"
    if quality_counts_path.exists():
        counts_df = pd.read_csv(quality_counts_path)
        count_sum = int(pd.to_numeric(counts_df.get("count", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
        add_check(
            "output_table_row_count_consistency",
            "pass" if count_sum == expected_rows else "fail",
            f"quality_category_counts sum={count_sum}; ledger rows={expected_rows}",
            "error" if count_sum != expected_rows else "info",
        )

    add_check(
        "ledger_required_columns",
        "pass" if set(LEDGER_COLUMNS).issubset(ledger.columns) else "fail",
        "Ledger has all required publication columns." if set(LEDGER_COLUMNS).issubset(ledger.columns) else "Ledger is missing required columns.",
        "error",
    )
    add_check("ledger_row_count", "pass" if expected_rows > 0 else "fail", f"Ledger rows: {expected_rows}", "error" if expected_rows == 0 else "info")

    validation = pd.DataFrame(checks)
    validation.to_csv(output_dir / "validation_checks.csv", index=False)
    failures = validation[(validation["status"] == "fail") & (validation["severity"] == "error")]
    if not failures.empty:
        log_lines.append(f"Validation found {len(failures)} error-level issue(s). See validation_checks.csv.")
        if strict:
            raise RuntimeError("Strict validation failed. See validation_checks.csv.")
    else:
        log_lines.append("Validation checks passed without error-level failures.")
    return validation


def _run_full_validation(
    output_dir: Path,
    panel_df: pd.DataFrame,
    log_lines: List[str],
) -> None:
    """Call the standalone validation module and write validation_report.json."""
    try:
        import sys as _sys
        _val_dir = str(PROJECT_ROOT / "validation")
        if _val_dir not in _sys.path:
            _sys.path.insert(0, _val_dir)
        from validate_publication_run import validate_run  # noqa: PLC0415

        panel_path: Optional[Path] = None
        for candidate in [
            output_dir / "case_study_panel.csv",
            PROJECT_ROOT / "config" / "case_study_panel_publication.csv",
            PROJECT_ROOT / "config" / "case_study_panel.csv",
        ]:
            if candidate.exists():
                panel_path = candidate
                break

        validate_run(output_dir=output_dir, panel_csv=panel_path, strict=False)
        log_lines.append("Full validation completed; validation_report.json written.")
    except Exception as exc:
        log_lines.append(f"Validation module error (non-fatal): {exc}")


def generate_manuscript_figures(
    output_dir: Path,
    ledger_path: Path,
    audit_dir: Path,
    panel_csv: Optional[Path],
    log_lines: List[str],
) -> None:
    """Delegate all figure generation to the visualisation module."""
    try:
        import sys as _sys
        _vis_dir = str(PROJECT_ROOT / "visualisation")
        if _vis_dir not in _sys.path:
            _sys.path.insert(0, _vis_dir)
        from make_all_publication_figures import generate_all_figures  # noqa: PLC0415

        figures_dir = output_dir / "manuscript_figures"
        supp_figs_dir = output_dir / "supplementary_figures"
        generate_all_figures(
            ledger_path=ledger_path,
            audit_dir=audit_dir,
            runs_dir=PROJECT_ROOT / "runs",
            output_dir=figures_dir,
            supp_dir=supp_figs_dir,
            panel_csv=panel_csv,
        )
        n_main = len(list(figures_dir.glob("Figure*.png")))
        n_supp = len(list(supp_figs_dir.glob("SuppFig*.png")))
        log_lines.append(f"Wrote {n_main} main figures and {n_supp} supplementary figures.")
    except Exception as exc:
        log_lines.append(f"Figure generation error: {exc}")
        import traceback
        log_lines.append(traceback.format_exc())


def _load_graph_lookup(graph_dir: Path) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """Return {(drug_lower, disease_lower): feature_row} from known+unknown CSVs."""
    lookup: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for name in ("graph_features_known.csv", "graph_features_unknown.csv"):
        path = graph_dir / name
        if not path.exists():
            path = PROJECT_ROOT / "graph" / name
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
            for _, row in df.iterrows():
                key = (str(row.get("Drug", "")).strip().lower(), str(row.get("Disease", "")).strip().lower())
                if key not in lookup:
                    lookup[key] = row.to_dict()
        except Exception:
            pass
    return lookup


def _load_lit_lookup(literature_dir: Path) -> Dict[str, Dict[str, Any]]:
    """Return {flat_drug_disease_key: classification_counts} from literature files."""
    lookup: Dict[str, Dict[str, Any]] = {}
    if not literature_dir.exists():
        return lookup

    for path in literature_dir.glob("classified_*.json"):
        records = read_json(path, [])
        if not isinstance(records, list):
            continue
        counts: Counter = Counter()
        for r in records:
            cat = norm_entity(r.get("category", ""))
            if cat in {"therapeutic", "adverse", "irrelevant"}:
                counts[cat] += 1
        counts["total"] = sum(counts.values())
        # Derive drug+disease from filename
        stem = path.stem
        for prefix in ("classified_pubmed_", "classified_articles_"):
            stem = stem.replace(prefix, "")
        stem = re.sub(r"_\d{8}_\d{6}$", "", stem)
        flat = re.sub(r"[^a-z0-9]", "", stem.lower())
        # Store by flat key; keep the file with most records if duplicate
        if flat not in lookup or counts["total"] > lookup[flat]["total"]:
            lookup[flat] = {
                "file": path.name,
                "therapeutic": counts.get("therapeutic", 0),
                "adverse": counts.get("adverse", 0),
                "irrelevant": counts.get("irrelevant", 0),
                "total": counts["total"],
            }
    return lookup


def _find_lit_record(drug: str, disease: str, lit_lookup: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    drug_flat = re.sub(r"[^a-z0-9]", "", drug.lower())
    disease_flat = re.sub(r"[^a-z0-9]", "", disease.lower())
    for flat_key_str, record in lit_lookup.items():
        if drug_flat and disease_flat and drug_flat in flat_key_str and disease_flat in flat_key_str:
            return record
    return None


def validate_pairs_stage(
    matched_path: Path,
    graph_dir: Path,
    runs_dir: Path,
    literature_dir: Path,
    log_lines: List[str],
) -> pd.DataFrame:
    """
    Stage 11 — Disease and drug-disease pair validation.

    For every unique (drug, disease) pair in the processed matched pairs:
    - checks presence in graph features (known vs unknown)
    - checks presence of a Bayesian run log
    - checks presence of a literature file
    - computes a per-pair coverage tier aligned with the publication ledger taxonomy
    """
    matched_rows = read_json(matched_path, []) if matched_path.exists() else []
    if not isinstance(matched_rows, list):
        matched_rows = []

    graph_lookup = _load_graph_lookup(graph_dir)
    lit_lookup = _load_lit_lookup(literature_dir)

    # Build run-log index
    run_index: Dict[Tuple[str, str], str] = {}
    if runs_dir.exists():
        for rpath in sorted(runs_dir.glob("run_*.json")):
            payload = read_json(rpath, {})
            if isinstance(payload, dict) and payload.get("drug") and payload.get("disease"):
                k = (payload["drug"].strip().lower(), payload["disease"].strip().lower())
                if k not in run_index:
                    run_index[k] = rpath.name

    # Unique drug-disease pairs from matched data
    seen: set = set()
    validation_rows: List[Dict[str, Any]] = []
    for row in matched_rows:
        drug = str(row.get("intervention", "")).strip()
        disease = str(row.get("condition", "")).strip()
        if not drug or not disease:
            continue
        key = (drug.lower(), disease.lower())
        if key in seen:
            continue
        seen.add(key)

        in_graph = key in graph_lookup
        graph_row = graph_lookup.get(key, {})
        label = graph_row.get("Label", None)
        in_known = (label == 1) if in_graph else False

        run_file = run_index.get(key)
        lit_record = _find_lit_record(drug, disease, lit_lookup)

        has_run = run_file is not None
        has_lit = lit_record is not None

        if has_run and in_graph:
            coverage_tier = "full_bayesian_audit"
        elif has_lit and in_graph:
            coverage_tier = "literature_and_graph"
        elif in_graph and not has_lit:
            coverage_tier = "graph_only"
        elif has_lit and not in_graph:
            coverage_tier = "literature_only"
        else:
            coverage_tier = "matched_pairs_only"

        validation_rows.append({
            "drug": drug,
            "disease": disease,
            "in_graph_features": in_graph,
            "is_known_pair_label1": bool(in_known),
            "has_bayesian_run_log": has_run,
            "has_literature_file": has_lit,
            "run_log_file": run_file or "",
            "literature_file": lit_record["file"] if lit_record else "",
            "literature_therapeutic": lit_record["therapeutic"] if lit_record else "",
            "literature_adverse": lit_record["adverse"] if lit_record else "",
            "literature_total": lit_record["total"] if lit_record else "",
            "graph_rw_score": round(float(graph_row.get("RandomWalkScore", 0)), 6) if in_graph else "",
            "graph_katz": round(float(graph_row.get("KatzSimilarity", 0)), 6) if in_graph else "",
            "coverage_tier": coverage_tier,
            "condition_mapping_method": row.get("condition_mapping_method", ""),
            "drug_mapping_method": row.get("intervention_mapping_method", row.get("drug_mapping_method", "")),
        })

    df = pd.DataFrame(validation_rows)
    log_lines.append(
        f"Pair validation: {len(validation_rows)} unique matched pairs | "
        f"full_bayesian={sum(1 for r in validation_rows if r['coverage_tier']=='full_bayesian_audit')} | "
        f"lit+graph={sum(1 for r in validation_rows if r['coverage_tier']=='literature_and_graph')} | "
        f"graph_only={sum(1 for r in validation_rows if r['coverage_tier']=='graph_only')}"
    )
    return df


def case_study_validation_stage(
    graph_dir: Path,
    runs_dir: Path,
    literature_dir: Path,
    panel_df: pd.DataFrame,
    log_lines: List[str],
) -> pd.DataFrame:
    """
    Stage 12 — Case study validation.

    Extracts and audits evidence for each predefined case study pair.
    Shows which evidence tiers are populated and what the key metrics are.
    """
    graph_lookup = _load_graph_lookup(graph_dir)
    lit_lookup = _load_lit_lookup(literature_dir)

    # Latest run-log per (drug, disease)
    run_data: Dict[Tuple[str, str], Dict[str, Any]] = {}
    if runs_dir.exists():
        for rpath in sorted(runs_dir.glob("run_*.json")):
            payload = read_json(rpath, {})
            if not isinstance(payload, dict):
                continue
            drug_r = str(payload.get("drug", "")).strip()
            disease_r = str(payload.get("disease", "")).strip()
            if not drug_r or not disease_r:
                continue
            k = (drug_r.lower(), disease_r.lower())
            ts = payload.get("timestamp", "")
            if k not in run_data or ts > run_data[k].get("timestamp", ""):
                run_data[k] = payload

    rows: List[Dict[str, Any]] = []
    for _, cs in panel_df.iterrows():
        drug = str(cs["drug"])
        disease = str(cs["disease"])
        key = (drug.lower(), disease.lower())

        in_graph = key in graph_lookup
        graph_row = graph_lookup.get(key, {})
        run_payload = run_data.get(key, {})
        comp = run_payload.get("components", {}) if run_payload else {}
        lit_record = _find_lit_record(drug, disease, lit_lookup)

        # Evidence availability flags
        has_run = bool(run_payload)
        has_lit_file = lit_record is not None

        # Literature counts — prefer run log, fall back to literature file
        if has_run and isinstance(comp.get("counts"), dict):
            counts = comp["counts"]
            therapeutic = int(counts.get("therapeutic", 0))
            adverse = int(counts.get("adverse", 0))
            irrelevant = int(counts.get("irrelevant", 0))
            m_articles = int(comp.get("M", therapeutic + adverse + irrelevant))
        elif has_lit_file:
            therapeutic = lit_record["therapeutic"]
            adverse = lit_record["adverse"]
            irrelevant = lit_record["irrelevant"]
            m_articles = lit_record["total"]
        else:
            therapeutic = adverse = irrelevant = m_articles = 0

        total_lit = therapeutic + adverse + irrelevant
        therapeutic_rate = round(therapeutic / total_lit, 4) if total_lit > 0 else None
        adverse_rate = round(adverse / total_lit, 4) if total_lit > 0 else None
        noise_rate = round(irrelevant / total_lit, 4) if total_lit > 0 else None

        # Safety / Bayesian from run log
        gamma = comp.get("gamma") if has_run else None
        p_raw = comp.get("p_raw") if has_run else None
        p_penalised = comp.get("p_penalised") if has_run else None
        p_final = comp.get("p_final") if has_run else None
        post_mean = comp.get("post_mean") if has_run else None
        p_likelihood = comp.get("p_likelihood") if has_run else None

        # Credible interval width
        post_a = float(comp.get("post_a", 0)) if has_run else 0
        post_b = float(comp.get("post_b", 0)) if has_run else 0
        ci_width = None
        if post_a > 0 and post_b > 0:
            try:
                from scipy.stats import beta as _beta
                ci_width = round(float(_beta.ppf(0.975, post_a, post_b)) - float(_beta.ppf(0.025, post_a, post_b)), 4)
            except Exception:
                pass

        # Evidence coverage tier
        if has_run and in_graph:
            evidence_tier = "full_bayesian_audit"
        elif has_lit_file and in_graph:
            evidence_tier = "literature_and_graph"
        elif in_graph:
            evidence_tier = "graph_only"
        elif has_lit_file:
            evidence_tier = "literature_only"
        else:
            evidence_tier = "not_in_system"

        rows.append({
            "drug": drug,
            "disease": disease,
            "case_type": cs.get("case_type", ""),
            "evidence_tier": evidence_tier,
            "in_graph_features": in_graph,
            "graph_label": int(graph_row.get("Label", -1)) if in_graph else "",
            "graph_rw_score": round(float(graph_row.get("RandomWalkScore", 0)), 6) if in_graph else "",
            "graph_katz": round(float(graph_row.get("KatzSimilarity", 0)), 6) if in_graph else "",
            "graph_structural_likelihood": round(float(graph_row.get("StructuralLikelihood", 0)), 6) if in_graph else "",
            "has_bayesian_run_log": has_run,
            "has_literature_file": has_lit_file,
            "m_articles_used": m_articles if m_articles > 0 else None,
            "therapeutic_count": therapeutic if total_lit > 0 else None,
            "adverse_count": adverse if total_lit > 0 else None,
            "irrelevant_count": irrelevant if total_lit > 0 else None,
            "therapeutic_rate": therapeutic_rate,
            "adverse_rate": adverse_rate,
            "irrelevant_noise_rate": noise_rate,
            "gamma_safety_overlap": round(float(gamma), 4) if gamma is not None else None,
            "p_raw": round(float(p_raw), 4) if p_raw is not None else None,
            "p_penalised": round(float(p_penalised), 4) if p_penalised is not None else None,
            "p_final_safety_adjusted": round(float(p_final), 4) if p_final is not None else None,
            "p_likelihood_graph": round(float(p_likelihood), 4) if p_likelihood is not None else None,
            "posterior_mean": round(float(post_mean), 4) if post_mean is not None else None,
            "credible_interval_width": ci_width,
            "expected_quality_category": cs.get("expected_quality", ""),
            "narrative": cs.get("narrative", ""),
        })
        log_lines.append(
            f"Case study {drug}/{disease}: tier={evidence_tier} | "
            f"lit={m_articles} articles (T={therapeutic}, A={adverse}, I={irrelevant}) | "
            f"post_mean={round(float(post_mean), 4) if post_mean is not None else 'N/A'}"
        )

    return pd.DataFrame(rows)


def generate_legacy_manuscript_tables(
    output_dir: Path,
    clinical_audit_path: Path,
    terminology_audit_path: Path,
    graph_audit_path: Path,
    case_study_path: Path,
    pair_level_path: Path,
    log_lines: List[str],
) -> None:
    """
    Generate four publication-ready manuscript tables from audit CSVs.

    Table 1 — Clinical trial extraction and MeSH mapping summary
    Table 2 — Graph and literature coverage by coverage tier
    Table 3 — Case study comparison (main results table)
    Table 4 — Evidence quality category distribution
    """
    tables_dir = output_dir / "manuscript_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # ── Table 1: Clinical trial and mapping summary ───────────────────────────
    rows_t1: List[Dict[str, Any]] = []
    if terminology_audit_path.exists():
        ta = pd.read_csv(terminology_audit_path)
        if not ta.empty:
            # Aggregate to condition vs drug mapping stats
            for entity in ("condition", "drug"):
                sub = ta[ta.get("entity_type", pd.Series()) == entity] if "entity_type" in ta.columns else pd.DataFrame()
                if sub.empty:
                    continue
                total = sub["count"].sum() if "count" in sub.columns else 0
                exact = sub.loc[sub["mapping_method"] == "exact_match", "count"].sum() if "mapping_method" in sub.columns else 0
                fuzzy_hi = sub.loc[sub["mapping_method"] == "fuzzy_high_confidence", "count"].sum() if "mapping_method" in sub.columns else 0
                unmapped = sub.loc[sub["mapping_method"].isin(["unmapped", "not_attempted"]), "count"].sum() if "mapping_method" in sub.columns else 0
                rows_t1.append({
                    "Entity": entity.capitalize(),
                    "Total_terms": int(total),
                    "Exact_match_n": int(exact),
                    "Exact_match_pct": round(100 * exact / total, 1) if total else 0,
                    "Fuzzy_high_confidence_n": int(fuzzy_hi),
                    "Fuzzy_high_confidence_pct": round(100 * fuzzy_hi / total, 1) if total else 0,
                    "Unmapped_n": int(unmapped),
                    "Unmapped_pct": round(100 * unmapped / total, 1) if total else 0,
                })
    if rows_t1:
        pd.DataFrame(rows_t1).to_csv(tables_dir / "Table1_mapping_quality.csv", index=False)
        log_lines.append("Wrote manuscript_tables/Table1_mapping_quality.csv")

    # ── Table 2: Graph + literature coverage by tier ──────────────────────────
    if pair_level_path.exists():
        pl = pd.read_csv(pair_level_path)
        if not pl.empty and "quality_flag" in pl.columns:
            tier_counts = pl["quality_flag"].value_counts().reset_index()
            tier_counts.columns = ["quality_flag", "n_pairs"]
            tier_counts["pct_of_audited"] = (tier_counts["n_pairs"] / len(pl) * 100).round(1)
            tier_counts.to_csv(tables_dir / "Table2_quality_tier_distribution.csv", index=False)
            log_lines.append("Wrote manuscript_tables/Table2_quality_tier_distribution.csv")

    # ── Table 3: Case study comparison ───────────────────────────────────────
    if case_study_path.exists():
        cs = pd.read_csv(case_study_path)
        if not cs.empty:
            cols = [
                "drug", "disease", "case_type",
                "m_articles_used", "therapeutic_rate", "adverse_rate", "irrelevant_noise_rate",
                "gamma_safety_overlap", "p_raw", "p_penalised", "p_final_safety_adjusted",
                "p_likelihood_graph", "posterior_mean", "credible_interval_width",
                "evidence_tier", "expected_quality_category",
            ]
            available = [c for c in cols if c in cs.columns]
            cs[available].to_csv(tables_dir / "Table3_case_study_comparison.csv", index=False)
            log_lines.append("Wrote manuscript_tables/Table3_case_study_comparison.csv")

    # ── Table 4: Drug-disease pair evidence audit (top 20 by posterior mean) ──
    if pair_level_path.exists():
        pl = pd.read_csv(pair_level_path)
        if not pl.empty and "posterior_mean" in pl.columns:
            top_cols = [
                "drug", "disease", "therapeutic_count", "adverse_count", "irrelevant_count",
                "gamma_safety_overlap", "posterior_mean", "credible_interval_width",
                "evidence_readiness_score", "quality_flag",
            ]
            available = [c for c in top_cols if c in pl.columns]
            top20 = pl[available].dropna(subset=["posterior_mean"]).sort_values("posterior_mean", ascending=False).head(20)
            top20.to_csv(tables_dir / "Table4_top_repurposing_candidates.csv", index=False)
            log_lines.append("Wrote manuscript_tables/Table4_top_repurposing_candidates.csv")


def generate_manuscript_tables(
    output_dir: Path,
    audit_dir: Path,
    pair_level_path: Path,
    panel_csv: Optional[Path],
    log_lines: List[str],
    # legacy positional compat (ignored)
    clinical_audit_path: Optional[Path] = None,
    terminology_audit_path: Optional[Path] = None,
    graph_audit_path: Optional[Path] = None,
    case_study_path: Optional[Path] = None,
) -> None:
    """Delegate all table generation to the reporting module."""
    try:
        import sys as _sys
        _rep_dir = str(PROJECT_ROOT / "reporting")
        if _rep_dir not in _sys.path:
            _sys.path.insert(0, _rep_dir)
        from make_all_publication_tables import generate_all_tables  # noqa: PLC0415

        tables_dir = output_dir / "manuscript_tables"
        supp_tables_dir = output_dir / "supplementary_tables"
        generate_all_tables(
            ledger_path=pair_level_path,
            audit_dir=audit_dir,
            output_dir=tables_dir,
            supp_dir=supp_tables_dir,
            panel_csv=panel_csv,
        )
        n_main = len(list(tables_dir.glob("Table*.csv")))
        n_supp = len(list(supp_tables_dir.glob("SuppTable*.csv")))
        log_lines.append(f"Wrote {n_main} main tables and {n_supp} supplementary tables.")
    except Exception as exc:
        log_lines.append(f"Table generation error: {exc}")
        import traceback
        log_lines.append(traceback.format_exc())


def write_config(args: argparse.Namespace, output_dir: Path, log_lines: List[str]) -> None:
    refresh_flags = {
        "refresh_trials": bool(args.refresh_trials),
        "refresh_mesh_mapping": bool(args.refresh_mesh_mapping),
        "refresh_graph": bool(getattr(args, "refresh_graph", False)),
        "refresh_literature": bool(args.refresh_literature),
        "refresh_safety": bool(args.refresh_safety),
        "run_bayesian": bool(args.run_bayesian),
        "reuse_classified_literature": bool(args.reuse_classified_literature),
        "record_audit_trail": bool(getattr(args, "record_audit_trail", True)),
    }
    legacy_flags = [k for k, v in refresh_flags.items() if not v]

    config = {
        "run_date": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(output_dir),
        "refresh_flags": refresh_flags,
        "legacy_artifacts_reused": legacy_flags,
        "legacy_artifact_warning": (
            f"The following refresh flags were false — existing local artifacts were reused: {legacy_flags}. "
            "Any mixed-provenance outputs are flagged in the ledger."
        ) if legacy_flags else "All artifacts freshly generated.",
        "therapeutic_areas": parse_areas(args.therapeutic_areas),
        "full_registry": bool(args.full_registry),
        "full_registry_db": args.full_registry_db,
        "resume_full_registry": bool(args.resume_full_registry),
        "input_run_dir": args.input_run_dir,
        # Reuse mode rebuilds graph/Bayesian calculations from local saved
        # classifications and safety outputs; it must not claim API or LLM use.
        "apis_used": {
            "clinicaltrials_gov": bool(args.refresh_trials),
            "mesh_download": bool(args.download_mesh),
            "pubmed_pmc": bool((args.refresh_literature or args.run_bayesian) and not args.reuse_classified_literature),
            "openfda_faers": bool((args.refresh_safety or args.run_bayesian) and not args.reuse_classified_literature),
            "openai_llm": bool((args.refresh_literature or args.refresh_safety or args.run_bayesian) and not args.reuse_classified_literature),
        },
        "mesh_version": "desc2026.xml",
        "case_study_panel": args.case_study_panel,
        "random_seed": args.random_seed,
        "thresholds": {
            "fuzzy_cutoff": args.fuzzy_cutoff,
            "token_jaccard_min": args.token_jaccard_min,
            "cutoff_year": args.cutoff_year,
            "mesh_exact_only": bool(args.mesh_exact_only),
        },
        "search_parameters": {
            "pubmed_max_articles": args.pubmed_max_articles,
            "pubmed_filter_level": args.pubmed_filter_level,
            "pubmed_years_back": args.pubmed_years_back,
            "pubmed_enable_pmc": bool(args.pubmed_enable_pmc),
            "force_pubmed_refresh": args.force_pubmed_refresh,
            "therapeutic_areas": parse_areas(args.therapeutic_areas),
            "full_registry": bool(args.full_registry),
        },
        "llm": {
            "model": args.llm_model,
            "classifier_schema_version": "five-way-evidence-v1",
            "batch_size": args.llm_batch_size,
            "delay_s": args.llm_delay_s,
            "max_retries": args.llm_max_retries,
        },
        "bayesian": {
            "run_bayesian": args.run_bayesian,
            "cmax": args.cmax,
            "tau": args.tau,
            "likelihood_strength": args.likelihood_strength,
            "likelihood_intercept": args.likelihood_intercept,
            "weights": (
                read_json((PROJECT_ROOT / args.graph_weights_path).resolve(), {}).get(
                    "raw_feature_weights", DEFAULT_WEIGHTS
                )
                if str(args.graph_weights_path or "").strip()
                else DEFAULT_WEIGHTS
            ),
            "weights_source": args.graph_weights_path or "legacy defaults",
        },
        "safety_penalty_settings": {
            "penalty_scale": 0.5,
            "fda_report_limit": 50,
            "gamma_source": "openFDA FAERS + LLM semantic overlap",
        },
        "input_file_paths": {
            "mesh_descriptor": args.mesh_path,
            "case_study_panel": args.case_study_panel,
            "processed_pairs": "processed_data/condition_drug_pairs.json",
            "graph_features_known": "graph/graph_features_known.csv",
            "graph_features_unknown": "graph/graph_features_unknown.csv",
            "runs_dir": "runs/",
            "literatures_dir": "literatures/",
        },
        "metric_definitions": METRIC_DEFINITIONS,
        "software_versions": {
            "python": sys.version,
            "platform": platform.platform(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "requests": requests.__version__,
            "networkx": _package_version("networkx"),
            "scipy": _package_version("scipy"),
            "scikit_learn": _package_version("scikit-learn"),
            "xgboost": _package_version("xgboost"),
        },
        "notes": [
            "Evidence readiness scores measure evidence quality and audit readiness, not clinical efficacy.",
            "When run_bayesian=true without reuse_classified_literature, the selected case-study panel is refreshed through PubMed/PMC retrieval, LLM classification, safety overlap, graph likelihood, and posterior scoring.",
            "When reuse_classified_literature=true, saved literature classifications and safety outputs are imported locally; only graph and Bayesian calculations are rebuilt.",
            "When run_bayesian=false, literature, safety, and Bayesian audits are summarized from existing local artifacts.",
            "If legacy_artifacts_reused is non-empty, the run mixed fresh and existing data — see legacy_artifact_warning.",
        ],
    }

    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    (output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    (logs_dir / "run_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    # Also write run_log at output root for backward compatibility
    (output_dir / "run_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    # Requirements snapshot
    try:
        import importlib.metadata as _meta  # noqa: PLC0415
        snap_lines = []
        for pkg in [
            "pandas", "numpy", "requests", "networkx", "scipy", "scikit-learn",
            "xgboost", "matplotlib", "openai",
        ]:
            try:
                snap_lines.append(f"{pkg}=={_meta.version(pkg)}")
            except Exception:
                pass
        (logs_dir / "requirements_snapshot.txt").write_text("\n".join(snap_lines) + "\n", encoding="utf-8")
    except Exception:
        pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run full data-quality pipeline rerun folder generation.")
    parser.add_argument("--refresh_trials", type=str_to_bool, default=True)
    parser.add_argument("--refresh_mesh_mapping", type=str_to_bool, default=True)
    parser.add_argument("--refresh_graph", type=str_to_bool, default=True)
    parser.add_argument("--refresh_literature", type=str_to_bool, default=False)
    parser.add_argument("--refresh_safety", type=str_to_bool, default=False)
    parser.add_argument("--run_bayesian", type=str_to_bool, default=False)
    parser.add_argument(
        "--reuse_classified_literature",
        type=str_to_bool,
        default=False,
        help="Recompute priors from saved classified JSON without rerunning LLM classification.",
    )
    parser.add_argument("--record_audit_trail", type=str_to_bool, default=True)
    parser.add_argument("--download_mesh", type=str_to_bool, default=False)
    parser.add_argument("--output_dir", default=f"outputs/publication_run_{datetime.now().strftime('%Y%m%d')}")
    parser.add_argument(
        "--input_run_dir",
        default="",
        help="Existing publication run supplying raw, mapped and classified artifacts when refresh flags are false.",
    )
    parser.add_argument(
        "--full_registry",
        type=str_to_bool,
        default=False,
        help="Use every current ClinicalTrials.gov study instead of therapeutic-area queries.",
    )
    parser.add_argument(
        "--full_registry_db",
        default="",
        help="Optional completed/resumable full-registry SQLite snapshot.",
    )
    parser.add_argument(
        "--resume_full_registry",
        type=str_to_bool,
        default=False,
        help="Resume an incomplete full-registry SQLite download; default starts a fresh snapshot.",
    )
    parser.add_argument("--registry_export_chunk_size", type=int, default=10000)
    parser.add_argument("--therapeutic_areas", default=",".join(DEFAULT_AREAS))
    parser.add_argument("--trials_per_area", type=int, default=None)
    parser.add_argument("--page_size", type=int, default=100)
    parser.add_argument("--cutoff_year", type=int, default=2020)
    parser.add_argument("--include_missing_dates", type=str_to_bool, default=False)
    parser.add_argument("--mesh_path", default="mesh_data/desc2026.xml")
    parser.add_argument("--fuzzy_cutoff", type=float, default=0.80)
    parser.add_argument("--token_jaccard_min", type=float, default=0.60)
    parser.add_argument(
        "--mesh_exact_only",
        type=str_to_bool,
        default=False,
        help="Use only exact descriptor/entry-term MeSH mappings; recommended for full-registry runs.",
    )
    parser.add_argument("--case_study_panel", default=DEFAULT_PANEL_PATH)
    parser.add_argument(
        "--pubmed_max_articles",
        type=int,
        default=0,
        help="Maximum PubMed records per pair; 0 retrieves every query result.",
    )
    parser.add_argument("--pubmed_filter_level", default="exact")
    parser.add_argument("--pubmed_years_back", type=int, default=0)
    parser.add_argument(
        "--pubmed_enable_pmc",
        type=str_to_bool,
        default=False,
        help="Augment a subset with PMC full text; disabled for uniform exhaustive classification.",
    )
    parser.add_argument("--force_pubmed_refresh", type=str_to_bool, default=True)
    parser.add_argument("--llm_model", default="gpt-4o-mini")
    parser.add_argument("--llm_batch_size", type=int, default=10)
    parser.add_argument("--llm_delay_s", type=float, default=0.2)
    parser.add_argument("--llm_max_retries", type=int, default=2)
    parser.add_argument("--cmax", type=float, default=200.0)
    parser.add_argument("--tau", type=float, default=25.0)
    parser.add_argument("--likelihood_strength", type=float, default=50.0)
    parser.add_argument("--likelihood_intercept", type=float, default=0.0)
    parser.add_argument("--graph_weights_path", default="")
    parser.add_argument("--max_graph_training_positives", type=int, default=20000)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--strict_validation", type=str_to_bool, default=False)
    parser.add_argument(
        "--validate_pairs",
        type=str_to_bool,
        default=True,
        help="Run stages 11 and 12: disease/drug pair validation and case study extraction.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = (PROJECT_ROOT / args.output_dir).resolve()
    input_run_dir = (
        (PROJECT_ROOT / args.input_run_dir).resolve()
        if str(args.input_run_dir or "").strip()
        else None
    )
    if input_run_dir is not None and not input_run_dir.exists():
        raise FileNotFoundError(f"Input run directory does not exist: {input_run_dir}")
    raw_dir = output_dir / "raw_trials"
    processed_dir = output_dir / "processed_data"
    graph_dir = output_dir / "graph"
    audit_dir = output_dir / "audit_files"
    ledgers_dir = output_dir / "ledgers"
    logs_dir = output_dir / "logs"
    manuscript_tables_dir = output_dir / "manuscript_tables"
    manuscript_figures_dir = output_dir / "manuscript_figures"
    supp_tables_dir = output_dir / "supplementary_tables"
    supp_figures_dir = output_dir / "supplementary_figures"

    for d in [audit_dir, ledgers_dir, logs_dir, manuscript_tables_dir,
              manuscript_figures_dir, supp_tables_dir, supp_figures_dir]:
        d.mkdir(parents=True, exist_ok=True)

    random.seed(args.random_seed)
    np.random.seed(args.random_seed)
    panel_df = load_case_study_panel((PROJECT_ROOT / args.case_study_panel).resolve())
    (output_dir / "case_study_panel.csv").parent.mkdir(parents=True, exist_ok=True)
    panel_df.to_csv(output_dir / "case_study_panel.csv", index=False)

    log_lines = [f"Started run at {datetime.now().isoformat(timespec='seconds')}"]
    mesh_path = (PROJECT_ROOT / args.mesh_path).resolve()

    if args.refresh_mesh_mapping:
        ensure_mesh(mesh_path, args.download_mesh, log_lines)
    else:
        log_lines.append("Skipping MeSH check (refresh_mesh_mapping=false)")

    if args.refresh_trials:
        clinical_audit = (
            fetch_full_registry_stage(args, raw_dir, log_lines)
            if args.full_registry
            else fetch_trials_stage(args, raw_dir, log_lines)
        )
    else:
        if input_run_dir is not None:
            raw_dir = input_run_dir / "raw_trials"
        elif not args.full_registry:
            raw_dir = PROJECT_ROOT / "data"
        clinical_audit = pd.DataFrame(
            [{"note": "refresh_trials=false; using existing data directory", "data_dir": rel(raw_dir)}]
        )

    clinical_audit.to_csv(audit_dir / "01_clinical_trial_extraction_audit.csv", index=False)
    clinical_audit.to_csv(output_dir / "01_clinical_trial_extraction_audit.csv", index=False)

    if args.refresh_mesh_mapping:
        terminology_audit, matched_path, unmatched_path = (
            mapping_stage_full_registry(
                raw_dir, processed_dir, mesh_path, args, log_lines
            )
            if args.full_registry
            else mapping_stage(raw_dir, processed_dir, mesh_path, args, log_lines)
        )
    else:
        processed_base = (
            (input_run_dir / "processed_data")
            if input_run_dir is not None
            else (processed_dir if args.full_registry else PROJECT_ROOT / "processed_data")
        )
        matched_path = processed_base / "condition_drug_pairs.json"
        unmatched_path = processed_base / "unmatched_pairs.json"
        # Keep a schema-valid provenance record even when mapping is reused.
        # The complete pair table remains the authoritative source for the
        # target-level mapping fields in the publication ledger.
        terminology_audit = pd.DataFrame(
            [
                {
                    "original_term": "reused existing condition-drug pair table",
                    "cleaned_term": "stored canonical pair mappings",
                    "mapped_term": "stored canonical pair mappings",
                    "mesh_id": "",
                    "mapping_method": "reused_existing_mapping",
                    "mapping_score": "",
                    "token_jaccard_score": "",
                    "mapping_status": "reused",
                    "failure_reason": "",
                    "source_path": rel(matched_path),
                    "reuse_note": "refresh_mesh_mapping=false; no mapping operation was run",
                }
            ]
        )

    terminology_audit.to_csv(audit_dir / "02_terminology_mapping_audit.csv", index=False)
    terminology_audit.to_csv(output_dir / "02_terminology_mapping_audit.csv", index=False)

    if args.refresh_graph:
        if args.full_registry:
            target_pairs = [
                (str(row["drug"]), str(row["disease"]))
                for _, row in panel_df.iterrows()
            ]
            graph_audit, known_graph_path, unknown_graph_path, weights_path = (
                build_targeted_full_graph(
                    matched_path=matched_path,
                    graph_dir=graph_dir,
                    target_pairs=target_pairs,
                    random_seed=args.random_seed,
                    max_training_positives=args.max_graph_training_positives,
                )
            )
            args.graph_weights_path = rel(weights_path)
            log_lines.append(
                f"Built scalable full-registry graph with target-edge withholding; weights={rel(weights_path)}"
            )
        else:
            target_pairs = [
                (str(row["drug"]), str(row["disease"]))
                for _, row in panel_df.iterrows()
            ]
            graph_audit, known_graph_path, unknown_graph_path, weights_path = graph_stage(
                matched_path, graph_dir, target_pairs, args, log_lines
            )
            args.graph_weights_path = rel(weights_path)
    else:
        graph_base = graph_dir if args.full_registry else PROJECT_ROOT / "graph"
        known_graph_path = graph_base / "graph_features_known.csv"
        unknown_graph_path = graph_base / "graph_features_unknown.csv"
        existing_weights = graph_base / "updated_graph_weights.json"
        if existing_weights.exists():
            args.graph_weights_path = rel(existing_weights)
        graph_audit = pd.DataFrame(
            [
                {
                    "metric": "refresh_graph",
                    "value": False,
                    "note": "Using existing graph feature artifacts.",
                }
            ]
        )
        log_lines.append("Skipping graph rebuild (refresh_graph=false); using existing graph feature artifacts.")
    # Write audit CSVs to audit_files/ and mirror roots for compat
    graph_audit.to_csv(audit_dir / "03_graph_construction_audit.csv", index=False)
    graph_audit.to_csv(output_dir / "03_graph_construction_audit.csv", index=False)

    # Prefer artifacts already produced inside this run when the expensive
    # Bayesian stage is intentionally not repeated.  Falling back to the
    # repository-level legacy folders preserves the historical workflow.
    effective_runs_dir = (
        input_run_dir / "runs"
        if input_run_dir is not None and not args.run_bayesian
        else (output_dir / "runs" if (output_dir / "runs").exists() else PROJECT_ROOT / "runs")
    )
    effective_lit_dir = (
        input_run_dir / "literatures"
        if input_run_dir is not None and not args.run_bayesian
        else (output_dir / "literatures" if (output_dir / "literatures").exists() else PROJECT_ROOT / "literatures")
    )
    if args.run_bayesian or args.refresh_literature or args.refresh_safety or args.reuse_classified_literature:
        if not args.run_bayesian:
            raise ValueError(
                "--refresh_literature true or --refresh_safety true requires --run_bayesian true "
                "so literature, safety, and posterior artifacts are refreshed together."
            )
        bayesian_panel_df, effective_runs_dir, effective_lit_dir = run_bayesian_panel_stage(
            panel_df=panel_df,
            known_graph_path=known_graph_path,
            unknown_graph_path=unknown_graph_path,
            output_dir=output_dir,
            args=args,
            log_lines=log_lines,
        )
        bayesian_panel_df.to_csv(audit_dir / "00_bayesian_panel_refresh_status.csv", index=False)
        bayesian_panel_df.to_csv(output_dir / "00_bayesian_panel_refresh_status.csv", index=False)

    for csv_name, stage_fn, stage_arg in [
        ("04_literature_retrieval_audit.csv", literature_retrieval_audit, effective_lit_dir),
        ("05_semantic_classification_audit.csv", semantic_classification_audit, effective_lit_dir),
        ("06_safety_overlap_audit.csv", safety_overlap_audit, effective_runs_dir),
        ("07_bayesian_uncertainty_audit.csv", bayesian_uncertainty_audit, effective_runs_dir),
    ]:
        df_stage = stage_fn(stage_arg)
        df_stage.to_csv(audit_dir / csv_name, index=False)
        df_stage.to_csv(output_dir / csv_name, index=False)

    report_args = SimpleNamespace(
        matched=rel(matched_path),
        unmatched=rel(unmatched_path),
        known_graph=rel(known_graph_path),
        unknown_graph=rel(unknown_graph_path),
        runs_dir=rel(effective_runs_dir),
        literature_dir=rel(effective_lit_dir),
        data_dir=rel(raw_dir),
        output_dir=rel(manuscript_tables_dir),
    )
    report_paths = generate_reports(report_args)

    # ── Stages 11 & 12: Pair validation and case study validation ────────────
    effective_graph_dir = graph_dir if (graph_dir / "graph_features_known.csv").exists() else PROJECT_ROOT / "graph"

    if args.validate_pairs:
        log_lines.append("Stage 11: Disease and drug-disease pair validation")
        validation_df = validate_pairs_stage(
            matched_path=matched_path,
            graph_dir=effective_graph_dir,
            runs_dir=effective_runs_dir,
            literature_dir=effective_lit_dir,
            log_lines=log_lines,
        )
        val_path = audit_dir / "11_disease_drug_pair_validation.csv"
        validation_df.to_csv(val_path, index=False)
        validation_df.to_csv(output_dir / "11_disease_drug_pair_validation.csv", index=False)
        log_lines.append(f"Saved {len(validation_df)} validated pairs.")

        if not validation_df.empty and "coverage_tier" in validation_df.columns:
            tier_summary = (
                validation_df["coverage_tier"]
                .value_counts().rename_axis("coverage_tier").reset_index(name="pair_count")
            )
            tier_summary["unique_drugs"] = validation_df.groupby("coverage_tier")["drug"].nunique().reindex(tier_summary["coverage_tier"]).values
            tier_summary["unique_diseases"] = validation_df.groupby("coverage_tier")["disease"].nunique().reindex(tier_summary["coverage_tier"]).values
            tier_summary.to_csv(audit_dir / "11b_coverage_tier_summary.csv", index=False)
            tier_summary.to_csv(output_dir / "11b_coverage_tier_summary.csv", index=False)

        log_lines.append("Stage 12: Case study validation")
        case_study_df = case_study_validation_stage(
            graph_dir=effective_graph_dir,
            runs_dir=effective_runs_dir,
            literature_dir=effective_lit_dir,
            panel_df=panel_df,
            log_lines=log_lines,
        )
        cs_path = audit_dir / "12_case_study_validation.csv"
        case_study_df.to_csv(cs_path, index=False)
        case_study_df.to_csv(output_dir / "12_case_study_validation.csv", index=False)
        log_lines.append(f"Saved {len(case_study_df)} case study rows.")
    else:
        val_path = audit_dir / "11_disease_drug_pair_validation.csv"
        cs_path = audit_dir / "12_case_study_validation.csv"
        log_lines.append("Stages 11/12 skipped (--validate_pairs false)")

    # ── Build and save publication ledger ─────────────────────────────────────
    publication_ledger = build_publication_ledger(
        pair_level_path=report_paths["pair_level"],
        terminology_audit_path=audit_dir / "02_terminology_mapping_audit.csv",
        validation_path=val_path if val_path.exists() else None,
        panel_df=panel_df,
        matched_path=matched_path,
    )
    # Primary location: ledgers/full_evidence_quality_ledger.csv
    ledger_path = ledgers_dir / "full_evidence_quality_ledger.csv"
    publication_ledger.to_csv(ledger_path, index=False)
    # Legacy mirror at root for backward compatibility
    publication_ledger.to_csv(output_dir / "08_full_evidence_quality_ledger.csv", index=False)
    publication_ledger.to_csv(audit_dir / "08_full_evidence_quality_ledger.csv", index=False)

    quality_counts_df = publication_ledger["quality_flag"].fillna("missing").value_counts().rename_axis("quality_flag").reset_index(name="count")
    quality_counts_df.to_csv(audit_dir / "09_quality_category_counts.csv", index=False)
    quality_counts_df.to_csv(output_dir / "09_quality_category_counts.csv", index=False)
    try:
        shutil.copyfile(report_paths["summary"], audit_dir / "10_summary_dashboard.csv")
        shutil.copyfile(report_paths["summary"], output_dir / "10_summary_dashboard.csv")
    except Exception:
        pass

    log_lines.append(f"Evidence-quality ledger: {len(publication_ledger)} rows → {rel(ledger_path)}")

    # ── Inline validation (fast checks) ─────────────────────────────────────
    validate_publication_outputs(
        ledger=publication_ledger,
        terminology_audit_path=audit_dir / "02_terminology_mapping_audit.csv",
        panel_df=panel_df,
        output_dir=output_dir,
        strict=args.strict_validation,
        log_lines=log_lines,
    )

    # ── Manuscript tables and figures ─────────────────────────────────────────
    panel_csv_for_outputs: Optional[Path] = (PROJECT_ROOT / args.case_study_panel).resolve()
    if not panel_csv_for_outputs.exists():
        panel_csv_for_outputs = None

    generate_manuscript_tables(
        output_dir=output_dir,
        audit_dir=audit_dir,
        pair_level_path=ledger_path,
        panel_csv=panel_csv_for_outputs,
        log_lines=log_lines,
    )
    generate_manuscript_figures(
        output_dir=output_dir,
        ledger_path=ledger_path,
        audit_dir=audit_dir,
        panel_csv=panel_csv_for_outputs,
        log_lines=log_lines,
    )

    log_lines.append(f"Finished run at {datetime.now().isoformat(timespec='seconds')}")
    # write_config must come before validation so run_config.json and run_log.txt exist
    write_config(args, output_dir, log_lines)

    # ── Full standalone validation module ────────────────────────────────────
    _run_full_validation(output_dir=output_dir, panel_df=panel_df, log_lines=log_lines)

    # ── Completion summary ────────────────────────────────────────────────────
    n_audited = len(publication_ledger)
    tier_col = "coverage_tier"
    n_full_bayes = int((publication_ledger[tier_col] == "full_bayesian_audit").sum()) if tier_col in publication_ledger.columns else 0
    n_graph_only = int((publication_ledger[tier_col] == "graph_only").sum()) if tier_col in publication_ledger.columns else 0
    n_lit_graph = int((publication_ledger[tier_col] == "literature_and_graph").sum()) if tier_col in publication_ledger.columns else 0
    n_cs_success = len(panel_df[panel_df.get("case_type", "").str.contains("successful", na=False)]) if "case_type" in panel_df.columns else 0
    n_cs_failed = len(panel_df[panel_df.get("case_type", "").str.contains("fail|conflict|safety", na=False, regex=True)]) if "case_type" in panel_df.columns else 0
    n_tables = len(list((output_dir / "manuscript_tables").glob("Table*.csv")))
    n_figures = len(list((output_dir / "manuscript_figures").glob("Figure*.png")))
    val_report_path = output_dir / "validation_report.json"
    val_status = "PASS"
    if val_report_path.exists():
        try:
            import json as _json
            vr = _json.loads(val_report_path.read_text())
            val_status = vr.get("overall_status", "UNKNOWN")
        except Exception:
            pass

    print(f"\n{'='*60}")
    print("  PUBLICATION PIPELINE — COMPLETION SUMMARY")
    print(f"{'='*60}")
    print(f"  Audited pairs total            : {n_audited}")
    print(f"  Full Bayesian audit pairs      : {n_full_bayes}")
    print(f"  Graph-only pairs               : {n_graph_only}")
    print(f"  Literature + graph pairs       : {n_lit_graph}")
    print(f"  Case study panel (successful)  : {n_cs_success}")
    print(f"  Case study panel (failed/conf.): {n_cs_failed}")
    print(f"  Manuscript tables generated    : {n_tables}")
    print(f"  Manuscript figures generated   : {n_figures}")
    print(f"  Validation status              : {val_status}")
    print(f"  Output directory               : {rel(output_dir)}")
    print(f"{'='*60}\n")
    print("  Key output files:")
    for label, path in [
        ("Full evidence ledger", ledger_path),
        ("run_config.json",  output_dir / "run_config.json"),
        ("validation_report.json", output_dir / "validation_report.json"),
        ("run_log.txt", logs_dir / "run_log.txt"),
        ("requirements_snapshot.txt", logs_dir / "requirements_snapshot.txt"),
    ]:
        size = f"({path.stat().st_size:,} B)" if path.exists() else "(not generated)"
        print(f"    {label:<35} {size}")


if __name__ == "__main__":
    main()
