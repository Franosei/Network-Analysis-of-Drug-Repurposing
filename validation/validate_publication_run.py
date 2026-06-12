"""
Validate a completed publication pipeline run and produce validation_report.json.

Usage (standalone):
    python validation/validate_publication_run.py \
        --output_dir outputs/publication_run_20260610 \
        --panel_csv  config/case_study_panel_publication.csv

Or call validate_run() from the master runner.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _read_json(path: Path) -> Any:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _pair_key(drug: Any, disease: Any) -> str:
    return f"{str(drug).strip().lower()}|{str(disease).strip().lower()}"


class Check:
    def __init__(self) -> None:
        self.checks: List[Dict[str, Any]] = []

    def add(self, name: str, status: str, detail: str, severity: str = "warning") -> None:
        self.checks.append({
            "check": name,
            "status": status,
            "severity": severity,
            "detail": detail,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        })

    def fail(self, name: str, detail: str) -> None:
        self.add(name, "fail", detail, "error")

    def warn(self, name: str, detail: str) -> None:
        self.add(name, "warn", detail, "warning")

    def ok(self, name: str, detail: str = "") -> None:
        self.add(name, "pass", detail, "info")

    @property
    def errors(self) -> List[str]:
        return [c["check"] for c in self.checks if c["status"] == "fail"]

    @property
    def warnings(self) -> List[str]:
        return [c["check"] for c in self.checks if c["status"] == "warn"]


# ── individual checks ─────────────────────────────────────────────────────────

def check_directory_structure(output_dir: Path, c: Check) -> None:
    required_dirs = [
        "audit_files",
        "ledgers",
        "manuscript_tables",
        "manuscript_figures",
        "supplementary_tables",
        "supplementary_figures",
        "logs",
    ]
    missing = [d for d in required_dirs if not (output_dir / d).is_dir()]
    if missing:
        c.warn("directory_structure", f"Missing expected subdirectories: {missing}")
    else:
        c.ok("directory_structure", "All required subdirectories present.")


def check_run_config(output_dir: Path, c: Check) -> None:
    cfg_path = output_dir / "run_config.json"
    if not cfg_path.exists():
        c.fail("run_config_present", "run_config.json is missing.")
        return
    cfg = _read_json(cfg_path)
    if not isinstance(cfg, dict):
        c.fail("run_config_valid", "run_config.json is not a valid JSON object.")
        return
    required_keys = [
        "run_date", "refresh_flags", "therapeutic_areas",
        "apis_used", "mesh_version", "thresholds",
        "llm", "bayesian", "software_versions", "output_dir",
    ]
    missing = [k for k in required_keys if k not in cfg]
    if missing:
        c.warn("run_config_completeness", f"run_config.json missing keys: {missing}")
    else:
        c.ok("run_config_completeness", "run_config.json has all required fields.")

    # Check for legacy/fresh mixing
    refresh = cfg.get("refresh_flags", {})
    if isinstance(refresh, dict):
        false_flags = [k for k, v in refresh.items() if v is False]
        if false_flags:
            c.warn(
                "legacy_artifacts_in_run",
                f"These refresh flags were false (existing artifacts reused): {false_flags}. "
                "Confirm these are flagged in the ledger and report.",
            )
        else:
            c.ok("legacy_artifacts_in_run", "All refresh flags were true; no legacy artifacts.")


def check_run_log(output_dir: Path, c: Check) -> None:
    log_path = output_dir / "logs" / "run_log.txt"
    if not log_path.exists():
        log_path = output_dir / "run_log.txt"
    if not log_path.exists():
        c.fail("run_log_present", "run_log.txt is missing.")
    else:
        size = log_path.stat().st_size
        c.ok("run_log_present", f"run_log.txt found ({size:,} bytes).")


def check_ledger(output_dir: Path, c: Check) -> pd.DataFrame:
    ledger_path = output_dir / "ledgers" / "full_evidence_quality_ledger.csv"
    if not ledger_path.exists():
        ledger_path = output_dir / "08_full_evidence_quality_ledger.csv"

    if not ledger_path.exists():
        c.fail("ledger_present", "full_evidence_quality_ledger.csv is missing.")
        return pd.DataFrame()

    ledger = pd.read_csv(ledger_path)
    c.ok("ledger_present", f"Ledger found at {ledger_path.name} with {len(ledger)} rows.")

    if ledger.empty:
        c.fail("ledger_nonempty", "Ledger is empty.")
        return ledger

    # Check required columns
    required_cols = {
        "drug", "disease",
        "drug_mapping_method", "disease_mapping_method",
        "drug_mapping_score", "disease_mapping_score",
        "trial_count", "articles_retrieved",
        "therapeutic_count", "adverse_count", "irrelevant_count",
        "safety_overlap_gamma", "posterior_mean",
        "credible_interval_width", "kl_divergence",
        "evidence_readiness_score", "quality_flag", "coverage_tier",
    }
    missing_cols = required_cols - set(ledger.columns)
    if missing_cols:
        c.fail("ledger_columns", f"Ledger missing required columns: {sorted(missing_cols)}")
    else:
        c.ok("ledger_columns", "All required ledger columns present.")

    # Check for duplicates
    dup_count = int(ledger.duplicated(subset=["drug", "disease"], keep=False).sum())
    if dup_count > 0:
        dups = ledger[ledger.duplicated(subset=["drug", "disease"], keep=False)][["drug", "disease"]].drop_duplicates()
        c.fail("ledger_no_duplicates", f"{dup_count} duplicate drug–disease rows. Examples: {dups.head(5).to_dict('records')}")
    else:
        c.ok("ledger_no_duplicates", "No duplicate drug–disease rows.")

    # Check mapping provenance columns
    prov_cols = {
        "drug_mapping_method", "disease_mapping_method",
        "drug_mapping_score", "disease_mapping_score",
        "drug_mapping_status", "disease_mapping_status",
    }
    missing_prov = prov_cols - set(ledger.columns)
    if missing_prov:
        c.warn("mapping_provenance", f"Mapping provenance columns missing: {sorted(missing_prov)}")
    else:
        # Check for legacy unmapped
        legacy = 0
        for col in ["drug_mapping_method", "disease_mapping_method"]:
            if col in ledger.columns:
                legacy += int((ledger[col] == "mapped_legacy_no_provenance").sum())
        if legacy > 0:
            c.warn(
                "legacy_mapping_present",
                f"{legacy} terms have mapping_method=mapped_legacy_no_provenance. "
                "Rerun with --refresh_mesh_mapping true to resolve.",
            )
        else:
            c.ok("mapping_provenance", "All mapping provenance columns present; no legacy mappings.")

    return ledger


def check_panel_coverage(ledger: pd.DataFrame, panel_csv: Optional[Path], c: Check) -> None:
    if panel_csv is None or not panel_csv.exists():
        c.warn("panel_coverage", "No case-study panel CSV provided; skipping panel checks.")
        return
    panel = pd.read_csv(panel_csv)
    if panel.empty or not {"drug", "disease"}.issubset(panel.columns):
        c.warn("panel_coverage", "Panel CSV is empty or missing drug/disease columns.")
        return

    ledger_keys = {_pair_key(r["drug"], r["disease"]) for _, r in ledger.iterrows()} if not ledger.empty else set()
    missing_pairs, no_bayes, no_lit, no_safety, no_uncertainty = [], [], [], [], []

    for _, row in panel.iterrows():
        key = _pair_key(row["drug"], row["disease"])
        if key not in ledger_keys:
            missing_pairs.append(f"{row['drug']} / {row['disease']}")
            continue
        lr = ledger[ledger.apply(lambda r: _pair_key(r["drug"], r["disease"]) == key, axis=1)].iloc[0]
        if pd.isna(lr.get("posterior_mean")):
            no_bayes.append(f"{row['drug']} / {row['disease']}")
        art = pd.to_numeric(lr.get("articles_retrieved", 0), errors="coerce")
        if pd.isna(art) or float(art) == 0:
            no_lit.append(f"{row['drug']} / {row['disease']}")
        if pd.isna(lr.get("safety_overlap_gamma")):
            no_safety.append(f"{row['drug']} / {row['disease']}")
        if pd.isna(lr.get("credible_interval_width")):
            no_uncertainty.append(f"{row['drug']} / {row['disease']}")

    status = "fail" if missing_pairs else "pass"
    detail = f"Missing from ledger: {missing_pairs[:10]}" if missing_pairs else f"All {len(panel)} panel pairs found in ledger."
    c.add("case_study_pairs_present", status, detail, "error" if missing_pairs else "info")

    if no_bayes:
        c.warn("bayesian_scores_for_panel", f"Missing posterior_mean for: {no_bayes[:10]}")
    else:
        c.ok("bayesian_scores_for_panel", "All panel pairs have posterior_mean.")
    if no_lit:
        c.warn("literature_for_panel", f"Zero articles_retrieved for: {no_lit[:10]}")
    else:
        c.ok("literature_for_panel", "All panel pairs have nonzero articles_retrieved.")
    if no_safety:
        c.warn("safety_for_panel", f"Missing safety_overlap_gamma for: {no_safety[:10]}")
    else:
        c.ok("safety_for_panel", "All panel pairs have safety_overlap_gamma.")
    if no_uncertainty:
        c.warn("uncertainty_for_panel", f"Missing credible_interval_width for: {no_uncertainty[:10]}")
    else:
        c.ok("uncertainty_for_panel", "All panel pairs have credible_interval_width.")


def check_manuscript_tables(output_dir: Path, c: Check) -> None:
    tables_dir = output_dir / "manuscript_tables"
    required_tables = [
        "Table1_data_source_and_preprocessing_readiness.csv",
        "Table2_evidence_quality_metric_definitions.csv",
        "Table3_case_study_comparison.csv",
        "Table4_top_evidence_ready_candidates_for_review.csv",
        "Table5_failed_or_conflicted_repurposing_examples.csv",
    ]
    missing = [t for t in required_tables if not (tables_dir / t).exists()]
    if missing:
        c.fail("manuscript_tables_present", f"Missing manuscript tables: {missing}")
    else:
        c.ok("manuscript_tables_present", f"All {len(required_tables)} required manuscript tables present.")


def check_manuscript_figures(output_dir: Path, c: Check) -> None:
    figs_dir = output_dir / "manuscript_figures"
    required_figs = [
        "Figure1_data_quality_pipeline_flow.png",
        "Figure2_preprocessing_readiness_flow.png",
        "Figure3_evidence_coverage_tiers.png",
        "Figure4_case_study_evidence_quality_heatmap.png",
        "Figure5_evidence_readiness_vs_uncertainty.png",
        "Figure6_literature_evidence_composition_case_studies.png",
        "Figure7_bayesian_prior_likelihood_posterior_example.png",
        "Figure8_robustness_and_sensitivity_summary.png",
    ]
    missing = [f for f in required_figs if not (figs_dir / f).exists()]
    if missing:
        c.fail("manuscript_figures_present", f"Missing manuscript figures: {missing}")
    else:
        c.ok("manuscript_figures_present", f"All {len(required_figs)} required manuscript figures present.")


def check_row_count_consistency(output_dir: Path, ledger: pd.DataFrame, c: Check) -> None:
    if ledger.empty:
        return
    expected = len(ledger)
    quality_counts_path = output_dir / "audit_files" / "09_quality_category_counts.csv"
    if not quality_counts_path.exists():
        quality_counts_path = output_dir / "09_quality_category_counts.csv"
    if quality_counts_path.exists():
        counts_df = pd.read_csv(quality_counts_path)
        count_col = "count" if "count" in counts_df.columns else counts_df.columns[-1]
        count_sum = int(pd.to_numeric(counts_df[count_col], errors="coerce").fillna(0).sum())
        if count_sum != expected:
            c.warn(
                "row_count_consistency",
                f"quality_category_counts sum ({count_sum}) ≠ ledger rows ({expected}). "
                "Possible stale audit files.",
            )
        else:
            c.ok("row_count_consistency", f"Row counts consistent: {expected} rows.")
    else:
        c.warn("row_count_consistency", "09_quality_category_counts.csv not found; skipping count check.")


def check_requirements_snapshot(output_dir: Path, c: Check) -> None:
    for loc in [output_dir / "logs" / "requirements_snapshot.txt", output_dir / "requirements_snapshot.txt"]:
        if loc.exists():
            c.ok("requirements_snapshot", f"requirements_snapshot.txt found at {loc.name}.")
            return
    c.warn("requirements_snapshot", "requirements_snapshot.txt not found.")


# ── orchestrator ─────────────────────────────────────────────────────────────

def validate_run(
    output_dir: Path,
    panel_csv: Optional[Path] = None,
    strict: bool = False,
) -> Dict[str, Any]:
    c = Check()
    c.add("validation_started", "info", f"Validating run at {output_dir}", "info")

    check_directory_structure(output_dir, c)
    check_run_config(output_dir, c)
    check_run_log(output_dir, c)
    ledger = check_ledger(output_dir, c)
    check_panel_coverage(ledger, panel_csv, c)
    check_manuscript_tables(output_dir, c)
    check_manuscript_figures(output_dir, c)
    check_row_count_consistency(output_dir, ledger, c)
    check_requirements_snapshot(output_dir, c)

    error_count = len(c.errors)
    warn_count = len(c.warnings)
    overall = "PASS" if error_count == 0 else "FAIL"

    # Coverage summary
    coverage_summary: Dict[str, int] = {}
    if not ledger.empty and "coverage_tier" in ledger.columns:
        coverage_summary = ledger["coverage_tier"].fillna("unknown").value_counts().to_dict()

    report = {
        "run_output_dir": str(output_dir),
        "validation_timestamp": datetime.now().isoformat(timespec="seconds"),
        "overall_status": overall,
        "error_count": error_count,
        "warning_count": warn_count,
        "errors": c.errors,
        "warnings": c.warnings,
        "checks": c.checks,
        "ledger_row_count": len(ledger),
        "coverage_summary": coverage_summary,
    }

    report_path = output_dir / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nValidation report saved: {report_path}")
    print(f"  Overall: {overall}  |  Errors: {error_count}  |  Warnings: {warn_count}")
    if c.errors:
        print("  ERRORS:")
        for e in c.errors:
            print(f"    [FAIL] {e}")
    if c.warnings:
        print("  WARNINGS:")
        for w in c.warnings[:10]:
            print(f"    [WARN] {w}")

    if strict and error_count > 0:
        sys.exit(1)

    return report


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Validate a publication pipeline run.")
    p.add_argument("--output_dir", required=True, help="Path to the publication run output folder.")
    p.add_argument("--panel_csv", default=None, help="Path to case-study panel CSV.")
    p.add_argument("--strict", action="store_true", help="Exit with code 1 if any errors found.")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    validate_run(
        output_dir=Path(args.output_dir),
        panel_csv=Path(args.panel_csv) if args.panel_csv else None,
        strict=args.strict,
    )
