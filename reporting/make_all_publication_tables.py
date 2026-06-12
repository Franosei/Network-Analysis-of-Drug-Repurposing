"""
Generate all publication-ready manuscript tables from the evidence-quality ledger
and audit files produced by run_full_data_quality_pipeline.py.

Usage (standalone):
    python reporting/make_all_publication_tables.py \
        --ledger_path  outputs/publication_run_20260610/ledgers/full_evidence_quality_ledger.csv \
        --audit_dir    outputs/publication_run_20260610/audit_files \
        --output_dir   outputs/publication_run_20260610/manuscript_tables \
        --supp_dir     outputs/publication_run_20260610/supplementary_tables \
        --panel_csv    config/case_study_panel_publication.csv

Or call generate_all_tables() from the master runner.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _save(df: pd.DataFrame, path: Path, note: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"  Saved: {path.name}  ({len(df)} rows){' — ' + note if note else ''}")


def _panel_keys(panel: pd.DataFrame):
    if panel.empty or not {"drug", "disease"}.issubset(panel.columns):
        return set()
    return set(
        zip(panel["drug"].str.lower().str.strip(), panel["disease"].str.lower().str.strip())
    )


# ── Table 1 — data source and preprocessing readiness ───────────────────────

def table1_preprocessing_readiness(
    clinical_audit: pd.DataFrame,
    terminology_audit: pd.DataFrame,
    graph_audit: pd.DataFrame,
    output_dir: Path,
) -> None:
    rows: List[dict] = []

    def _metric(label: str, value, note: str = "") -> dict:
        return {"Metric": label, "Value": value, "Notes": note}

    # Clinical trial extraction
    if not clinical_audit.empty and "metric" in clinical_audit.columns:
        for _, r in clinical_audit.iterrows():
            rows.append(_metric(str(r.get("metric", "")), r.get("value", ""), str(r.get("note", ""))))
    else:
        rows.append(_metric("Clinical trial data", "See run_log.txt", "refresh_trials=false — existing data used"))

    # Terminology mapping
    if not terminology_audit.empty and "mapping_method" in terminology_audit.columns:
        total = int(terminology_audit.get("count", pd.Series(dtype=float)).fillna(0).sum())
        for method, grp in terminology_audit.groupby("mapping_method"):
            n = int(grp["count"].sum()) if "count" in grp.columns else len(grp)
            pct = round(100 * n / total, 1) if total else 0.0
            rows.append(_metric(f"Mapping method: {method}", f"{n} ({pct}%)", ""))
        rows.append(_metric("Total terms processed", total, ""))
    else:
        rows.append(_metric("Terminology mapping", "See 02_terminology_mapping_audit.csv", ""))

    # Graph construction
    if not graph_audit.empty and "metric" in graph_audit.columns:
        for _, r in graph_audit.iterrows():
            rows.append(_metric(str(r.get("metric", "")), r.get("value", ""), str(r.get("note", ""))))
    else:
        rows.append(_metric("Graph construction", "See 03_graph_construction_audit.csv", ""))

    df = pd.DataFrame(rows)
    df = df[df["Metric"].astype(str).str.strip() != ""]
    _save(df, output_dir / "Table1_data_source_and_preprocessing_readiness.csv")


# ── Table 2 — metric definitions ────────────────────────────────────────────

def table2_metric_definitions(output_dir: Path) -> None:
    definitions = [
        {
            "Metric": "evidence_readiness_score",
            "Source": "data_quality/composite_quality.py",
            "Formula_Rule": "Weighted sum of 7 component scores × 100; all components clamped to [0,1]",
            "Interpretation": "0–100 composite quality score. Higher = better evidence readiness.",
            "Direction": "Higher = better",
        },
        {
            "Metric": "coverage_tier",
            "Source": "data_quality/quality_flags.py",
            "Formula_Rule": "Rule-based: full_bayesian_audit if has_mapping+graph+literature+safety+Bayesian",
            "Interpretation": "Completeness of the local evidence record for a drug–disease pair.",
            "Direction": "full_bayesian_audit is most complete",
        },
        {
            "Metric": "therapeutic_ratio",
            "Source": "code/pubmed_utils.py — LLM classification",
            "Formula_Rule": "therapeutic_count / (therapeutic + adverse + irrelevant)",
            "Interpretation": "Fraction of articles classified as therapeutically relevant.",
            "Direction": "Higher = more therapeutic signal",
        },
        {
            "Metric": "adverse_burden",
            "Source": "code/pubmed_utils.py — LLM classification",
            "Formula_Rule": "adverse_count / (therapeutic + adverse + irrelevant)",
            "Interpretation": "Fraction of articles flagging adverse or safety concern.",
            "Direction": "Lower = cleaner evidence",
        },
        {
            "Metric": "irrelevant_noise_rate",
            "Source": "code/pubmed_utils.py — LLM classification",
            "Formula_Rule": "irrelevant_count / all_classified_articles",
            "Interpretation": "Proportion of retrieved articles not relevant to the drug–disease question.",
            "Direction": "Lower = less noise",
        },
        {
            "Metric": "literature_completeness_score",
            "Source": "run_full_data_quality_pipeline.py",
            "Formula_Rule": "min(articles_retrieved/50, 1) × (usable/retrieved if retrieved>0 else 0)",
            "Interpretation": "0–1 score for article retrieval completeness.",
            "Direction": "Higher = more complete",
        },
        {
            "Metric": "safety_overlap_gamma",
            "Source": "code/side_effect_updater.py — openFDA + LLM",
            "Formula_Rule": "LLM-scored semantic overlap between adverse-event terms and disease symptoms",
            "Interpretation": "Overlap between drug's known side effects and the target disease's phenotype.",
            "Direction": "Lower = less safety concern",
        },
        {
            "Metric": "safety_penalty",
            "Source": "code/side_effect_updater.py",
            "Formula_Rule": "p_penalised = p_raw × (1 − penalty_scale × gamma); penalty_scale=0.5",
            "Interpretation": "Multiplicative penalty applied to prior due to adverse-event overlap.",
            "Direction": "Smaller penalty = less safety concern",
        },
        {
            "Metric": "structural_consistency_score",
            "Source": "code/network_builder.py",
            "Formula_Rule": "Normalised composite of graph_distance, random_walk_score, katz_similarity, structural_likelihood",
            "Interpretation": "How strongly the network topology supports the drug–disease connection.",
            "Direction": "Higher = stronger structural support",
        },
        {
            "Metric": "prior_mean",
            "Source": "code/bayesian_predictor.py",
            "Formula_Rule": "p_final from LLM semantic prior after safety penalty",
            "Interpretation": "Prior belief in repurposing plausibility before graph evidence update.",
            "Direction": "Higher = stronger prior belief",
        },
        {
            "Metric": "posterior_mean",
            "Source": "code/bayesian_predictor.py — Beta posterior",
            "Formula_Rule": "α / (α + β) from Beta(α,β) fusion of prior and graph likelihood",
            "Interpretation": "Final posterior probability of repurposing plausibility.",
            "Direction": "Higher = stronger posterior belief",
        },
        {
            "Metric": "credible_interval_width",
            "Source": "code/bayesian_predictor.py",
            "Formula_Rule": "Beta.ppf(0.975, α, β) − Beta.ppf(0.025, α, β)",
            "Interpretation": "Width of the 95% posterior credible interval. Wide = uncertain.",
            "Direction": "Lower = more certain",
        },
        {
            "Metric": "kl_divergence",
            "Source": "code/bayesian_predictor.py",
            "Formula_Rule": "KL(posterior || prior) between Beta distributions",
            "Interpretation": "Information gain from prior to posterior. Low = evidence barely updated beliefs.",
            "Direction": "Higher = more informative evidence",
        },
        {
            "Metric": "mean_shift",
            "Source": "run_full_data_quality_pipeline.py",
            "Formula_Rule": "posterior_mean − prior_mean",
            "Interpretation": "Direction and magnitude of Bayesian update.",
            "Direction": "Positive = evidence supports repurposing; negative = against",
        },
        {
            "Metric": "uncertainty_level",
            "Source": "run_full_data_quality_pipeline.py",
            "Formula_Rule": "low if CI_width<0.15; moderate if 0.15–0.35; high if >0.35",
            "Interpretation": "Categorical uncertainty label for the posterior estimate.",
            "Direction": "low = well-constrained",
        },
        {
            "Metric": "drug_mapping_method",
            "Source": "code/condition_drug_pairs.py — MeSH fuzzy matching",
            "Formula_Rule": "exact_match | fuzzy_high_confidence (≥0.80) | fuzzy_low_confidence (<0.80) | unmapped | mapped_legacy_no_provenance",
            "Interpretation": "How the drug term was normalised to a MeSH descriptor.",
            "Direction": "exact_match is most reliable",
        },
        {
            "Metric": "disease_mapping_method",
            "Source": "code/condition_drug_pairs.py — MeSH fuzzy matching",
            "Formula_Rule": "exact_match | fuzzy_high_confidence (≥0.80) | fuzzy_low_confidence (<0.80) | unmapped | mapped_legacy_no_provenance",
            "Interpretation": "How the disease term was normalised to a MeSH descriptor.",
            "Direction": "exact_match is most reliable",
        },
        {
            "Metric": "quality_flag",
            "Source": "data_quality/quality_flags.py",
            "Formula_Rule": "Rule-based classification using coverage, noise rate, safety gamma, CI width, mapping status",
            "Interpretation": "Human-readable evidence-quality category for the drug–disease pair.",
            "Direction": "High evidence quality = strongest evidence basis",
        },
    ]
    df = pd.DataFrame(definitions)
    _save(df, output_dir / "Table2_evidence_quality_metric_definitions.csv")


# ── Table 3 — case study comparison ─────────────────────────────────────────

def table3_case_study_comparison(
    ledger: pd.DataFrame,
    panel: pd.DataFrame,
    output_dir: Path,
) -> None:
    if ledger.empty:
        _save(pd.DataFrame(), output_dir / "Table3_case_study_comparison.csv", "empty ledger")
        return

    # Filter to panel pairs
    keys = _panel_keys(panel)
    df = ledger.copy()
    if keys:
        df["_key"] = list(zip(df["drug"].str.lower().str.strip(), df["disease"].str.lower().str.strip()))
        df = df[df["_key"].isin(keys)].drop(columns=["_key"])

    if df.empty:
        df = ledger.head(30)

    # Merge panel metadata
    if not panel.empty and {"drug", "disease"}.issubset(panel.columns):
        panel_meta = panel[["drug", "disease"] + [c for c in ["case_type", "expected_quality", "narrative"] if c in panel.columns]].copy()
        panel_meta["drug"] = panel_meta["drug"].str.lower().str.strip()
        panel_meta["disease"] = panel_meta["disease"].str.lower().str.strip()
        df_merge = df.copy()
        df_merge["drug_lc"] = df_merge["drug"].str.lower().str.strip()
        df_merge["disease_lc"] = df_merge["disease"].str.lower().str.strip()
        df = df_merge.merge(panel_meta, left_on=["drug_lc", "disease_lc"], right_on=["drug", "disease"], how="left", suffixes=("", "_panel")).drop(columns=["drug_lc", "disease_lc"])
        for col in ["drug_panel", "disease_panel"]:
            if col in df.columns:
                df.drop(columns=[col], inplace=True)

    cols = [
        "drug", "disease",
        "case_type", "expected_quality",
        "coverage_tier",
        "drug_mapping_method", "disease_mapping_method",
        "drug_mapping_score", "disease_mapping_score",
        "articles_retrieved", "usable_articles",
        "therapeutic_ratio", "adverse_burden", "irrelevant_noise_rate",
        "safety_overlap_gamma",
        "structural_consistency_score",
        "prior_mean", "posterior_mean",
        "credible_interval_lower", "credible_interval_upper", "credible_interval_width",
        "kl_divergence", "mean_shift",
        "evidence_readiness_score", "uncertainty_level",
        "quality_flag",
        "final_interpretation",
    ]
    available = [c for c in cols if c in df.columns]
    _save(df[available], output_dir / "Table3_case_study_comparison.csv")


# ── Table 4 — top evidence-ready candidates for review ──────────────────────

def table4_top_candidates(ledger: pd.DataFrame, output_dir: Path) -> None:
    if ledger.empty:
        _save(pd.DataFrame(), output_dir / "Table4_top_evidence_ready_candidates_for_review.csv", "empty ledger")
        return

    # Apply minimum quality gates
    df = ledger.copy()
    for col in ["evidence_readiness_score", "credible_interval_width"]:
        if col not in df.columns:
            df[col] = float("nan")

    df["_score"] = pd.to_numeric(df["evidence_readiness_score"], errors="coerce")
    df["_ci_w"] = pd.to_numeric(df["credible_interval_width"], errors="coerce")

    # Gates: score ≥ 40 and CI width ≤ 0.50 (some uncertainty is acceptable)
    thresholded = df[(df["_score"] >= 40) & (df["_ci_w"] <= 0.50)]
    if thresholded.empty:
        thresholded = df.nlargest(25, "_score")

    # Exclude safety-conflicted or failed flags
    exclude_flags = {
        "Safety-conflicted evidence", "Safety-concerning",
        "Insufficient evidence", "Literature noise dominated",
    }
    if "quality_flag" in thresholded.columns:
        filtered = thresholded[~thresholded["quality_flag"].isin(exclude_flags)]
        if filtered.empty:
            filtered = thresholded
    else:
        filtered = thresholded

    top = filtered.sort_values("_score", ascending=False).drop(columns=["_score", "_ci_w"]).head(25)

    cols = [
        "drug", "disease", "coverage_tier",
        "drug_mapping_method", "disease_mapping_method",
        "articles_retrieved", "therapeutic_ratio",
        "safety_overlap_gamma", "structural_consistency_score",
        "posterior_mean", "credible_interval_width",
        "kl_divergence", "evidence_readiness_score", "uncertainty_level",
        "quality_flag", "final_interpretation",
    ]
    available = [c for c in cols if c in top.columns]
    _save(top[available], output_dir / "Table4_top_evidence_ready_candidates_for_review.csv")


# ── Table 5 — failed or conflicted examples ──────────────────────────────────

def table5_failed_conflicted(ledger: pd.DataFrame, output_dir: Path) -> None:
    if ledger.empty:
        _save(pd.DataFrame(), output_dir / "Table5_failed_or_conflicted_repurposing_examples.csv", "empty ledger")
        return

    conflicted_flags = {
        "Safety-conflicted evidence",
        "Safety-concerning",
        "Literature-conflicted",
        "Literature noise dominated",
        "Terminology uncertainty",
        "Insufficient evidence",
        "failed_repurposing",
    }
    df = ledger.copy()
    if "quality_flag" in df.columns:
        mask = df["quality_flag"].isin(conflicted_flags)
    elif "case_type" in df.columns:
        mask = df["case_type"].isin({"failed_repurposing", "safety_conflicted", "noisy_evidence"})
    else:
        mask = pd.Series([False] * len(df))

    conflicted = df[mask].copy()
    if conflicted.empty:
        # Fall back to highest CI width (most uncertain)
        if "credible_interval_width" in df.columns:
            conflicted = df.nlargest(20, "credible_interval_width")
        else:
            conflicted = df.tail(20)

    conflicted = conflicted.sort_values(
        ["quality_flag", "credible_interval_width"] if "credible_interval_width" in conflicted.columns else ["quality_flag"],
        ascending=[True, False] if "credible_interval_width" in conflicted.columns else [True],
    )

    cols = [
        "drug", "disease", "case_type", "coverage_tier",
        "drug_mapping_status", "disease_mapping_status",
        "articles_retrieved", "therapeutic_ratio", "adverse_burden",
        "irrelevant_noise_rate", "safety_overlap_gamma",
        "posterior_mean", "credible_interval_width",
        "kl_divergence", "evidence_readiness_score",
        "quality_flag", "final_interpretation",
    ]
    available = [c for c in cols if c in conflicted.columns]
    _save(conflicted[available], output_dir / "Table5_failed_or_conflicted_repurposing_examples.csv")


# ── Supplementary tables ─────────────────────────────────────────────────────

def supp_full_pair_audit(ledger: pd.DataFrame, supp_dir: Path) -> None:
    if ledger.empty:
        return
    _save(ledger, supp_dir / "SuppTable_full_pair_audit.csv", "all pairs, all columns")


def supp_coverage_tier_summary(ledger: pd.DataFrame, supp_dir: Path) -> None:
    if ledger.empty or "coverage_tier" not in ledger.columns:
        return
    summary = (
        ledger["coverage_tier"].fillna("missing").value_counts()
        .rename_axis("coverage_tier").reset_index(name="pair_count")
    )
    if "drug" in ledger.columns:
        summary["unique_drugs"] = ledger.groupby("coverage_tier")["drug"].nunique().reindex(summary["coverage_tier"]).values
    if "disease" in ledger.columns:
        summary["unique_diseases"] = ledger.groupby("coverage_tier")["disease"].nunique().reindex(summary["coverage_tier"]).values
    _save(summary, supp_dir / "SuppTable_coverage_tier_summary.csv")


def supp_quality_flag_definitions(supp_dir: Path) -> None:
    rows = [
        {"quality_flag": "High evidence quality",
         "definition": "Full Bayesian audit passed; low noise, low safety conflict, posterior well-constrained."},
        {"quality_flag": "Moderate evidence quality",
         "definition": "Good evidence coverage with some uncertainty; further review warranted."},
        {"quality_flag": "Insufficient evidence",
         "definition": "Too few articles or trials to make a reliable assessment."},
        {"quality_flag": "Safety-conflicted evidence",
         "definition": "Safety overlap gamma is high and adverse burden is notable."},
        {"quality_flag": "Safety-concerning",
         "definition": "Safety signal present but not dominant; should be monitored."},
        {"quality_flag": "Literature-conflicted",
         "definition": "Both therapeutic and adverse signals present in literature at similar rates."},
        {"quality_flag": "Literature noise dominated",
         "definition": "Irrelevant noise rate is high; retrieved literature is poorly specific."},
        {"quality_flag": "Terminology uncertainty",
         "definition": "Drug or disease mapping is low-confidence or unmapped."},
    ]
    _save(pd.DataFrame(rows), supp_dir / "SuppTable_quality_flag_definitions.csv")


# ── main entry point ─────────────────────────────────────────────────────────

def generate_all_tables(
    ledger_path: Path,
    audit_dir: Path,
    output_dir: Path,
    supp_dir: Path,
    panel_csv: Optional[Path] = None,
) -> List[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    ledger = _read_csv(ledger_path)
    clinical_audit = _read_csv(audit_dir / "01_clinical_trial_extraction_audit.csv")
    terminology_audit = _read_csv(audit_dir / "02_terminology_mapping_audit.csv")
    graph_audit = _read_csv(audit_dir / "03_graph_construction_audit.csv")
    panel = _read_csv(panel_csv) if panel_csv else pd.DataFrame()

    print("Generating manuscript tables…")
    table1_preprocessing_readiness(clinical_audit, terminology_audit, graph_audit, output_dir)
    table2_metric_definitions(output_dir)
    table3_case_study_comparison(ledger, panel, output_dir)
    table4_top_candidates(ledger, output_dir)
    table5_failed_conflicted(ledger, output_dir)

    print("Generating supplementary tables…")
    supp_full_pair_audit(ledger, supp_dir)
    supp_coverage_tier_summary(ledger, supp_dir)
    supp_quality_flag_definitions(supp_dir)

    generated = sorted(str(p.name) for p in output_dir.glob("Table*.csv"))
    supp_generated = sorted(str(p.name) for p in supp_dir.glob("SuppTable*.csv"))
    print(f"  {len(generated)} main tables, {len(supp_generated)} supplementary tables written.")
    return generated


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate all publication tables.")
    p.add_argument("--ledger_path", required=True)
    p.add_argument("--audit_dir", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--supp_dir", default=None)
    p.add_argument("--panel_csv", default=None)
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    out = Path(args.output_dir)
    supp = Path(args.supp_dir) if args.supp_dir else out.parent / "supplementary_tables"
    generate_all_tables(
        ledger_path=Path(args.ledger_path),
        audit_dir=Path(args.audit_dir),
        output_dir=out,
        supp_dir=supp,
        panel_csv=Path(args.panel_csv) if args.panel_csv else None,
    )
