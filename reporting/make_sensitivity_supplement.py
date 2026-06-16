"""
Generate offline supplementary sensitivity tables for the Bayesian model.

This script does not rerun ClinicalTrials.gov, MeSH mapping, PubMed, OpenAI, or
graph construction. It reuses the retained publication-run ledger and varies
only the modelling constants requested by reviewers:

  - cmax
  - tau
  - likelihood strength lambda
  - adverse-evidence weight in p_penalised
  - safety coefficient in p_final

Outputs:
  - SuppTable_sensitivity_parameter_summary.csv
  - SuppTable_sensitivity_pair_level.csv
  - SuppText_sensitivity_analysis.md
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from scipy.stats import beta as beta_dist

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from data_quality.quality_flags import quality_flag  # noqa: E402


FEATURE_COLUMNS = {
    "GraphDistanceToIndication": "graph_distance",
    "RandomWalkScore": "random_walk_score",
    "StructuralLikelihood": "structural_likelihood",
    "PreferentialAttachment": "preferential_attachment",
    "KatzSimilarity": "katz_similarity",
}


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


def clamp_prob(value: Any, default: float = 0.5, eps: float = 1e-6) -> float:
    out = safe_float(value, default)
    if out is None:
        out = default
    out = max(0.0, min(1.0, float(out)))
    return min(max(out, eps), 1.0 - eps)


def sigmoid(value: float) -> float:
    if value >= 0:
        exp_neg = math.exp(-value)
        return 1.0 / (1.0 + exp_neg)
    exp_pos = math.exp(value)
    return exp_pos / (1.0 + exp_pos)


def beta_params_from_prob(probability: float, strength: float) -> tuple[float, float]:
    p = clamp_prob(probability)
    s = max(float(strength), 0.0)
    return 1.0 + s * p, 1.0 + s * (1.0 - p)


def concentration_c(article_count: int, cmax: float, tau: float) -> float:
    m = int(max(article_count, 0))
    cmax = float(max(cmax, 0.0))
    tau = float(max(tau, 1e-6))
    return cmax * (1.0 - math.exp(-m / tau))


def load_run_config(run_dir: Path) -> Dict[str, Any]:
    config_path = run_dir / "run_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    bayes = config.get("bayesian", {})
    return {
        "cmax": float(bayes.get("cmax", 200.0)),
        "tau": float(bayes.get("tau", 25.0)),
        "likelihood_strength": float(bayes.get("likelihood_strength", 50.0)),
        "likelihood_intercept": float(bayes.get("likelihood_intercept", 0.0)),
        "weights": dict(bayes.get("weights", {})),
        "adverse_weight": 2.0,
        "safety_coefficient": float(
            config.get("safety_penalty_settings", {}).get("penalty_scale", 0.5)
        ),
    }


def pair_key(drug: Any, disease: Any) -> tuple[str, str]:
    return str(drug or "").strip().lower(), str(disease or "").strip().lower()


def load_run_components(run_dir: Path) -> Dict[tuple[str, str], Dict[str, Any]]:
    runs_dir = run_dir / "runs"
    components: Dict[tuple[str, str], Dict[str, Any]] = {}
    if not runs_dir.exists():
        return components

    for path in sorted(runs_dir.glob("run_*.json")):
        if path.stat().st_size == 0:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        drug = payload.get("drug")
        disease = payload.get("disease")
        comp = payload.get("components")
        if not drug or not disease or not isinstance(comp, dict):
            continue
        components[pair_key(drug, disease)] = {
            "components": comp,
            "case_type": payload.get("case_type", ""),
            "path": str(path),
        }
    return components


def graph_probability(row: pd.Series, config: Dict[str, Any]) -> tuple[float, float]:
    score = float(config["likelihood_intercept"])
    for feature_name, ledger_column in FEATURE_COLUMNS.items():
        weight = float(config["weights"].get(feature_name, 0.0))
        value = safe_float(row.get(ledger_column), 0.0) or 0.0
        score += weight * value
    return clamp_prob(sigmoid(score)), score


def adjusted_readiness_score(
    row: pd.Series,
    scenario_ci_width: float,
    baseline_ci_width: float,
) -> float:
    baseline_score = safe_float(row.get("evidence_readiness_score"), 0.0) or 0.0
    adjusted = baseline_score + 15.0 * (baseline_ci_width - scenario_ci_width)
    return round(max(0.0, min(100.0, adjusted)), 3)


def recompute_pair(
    row: pd.Series,
    config: Dict[str, Any],
    run_record: Optional[Dict[str, Any]],
    *,
    cmax: float,
    tau: float,
    likelihood_strength: float,
    adverse_weight: float,
    safety_coefficient: float,
) -> Dict[str, Any]:
    comp = (run_record or {}).get("components", {})
    counts = comp.get("counts") if isinstance(comp, dict) else None
    if isinstance(counts, dict):
        therapeutic = int(safe_float(counts.get("therapeutic"), 0.0) or 0)
        adverse = int(safe_float(counts.get("adverse"), 0.0) or 0)
        irrelevant = int(safe_float(counts.get("irrelevant"), 0.0) or 0)
        article_count = int(safe_float(comp.get("M"), therapeutic + adverse + irrelevant) or 0)
    else:
        therapeutic = int(safe_float(row.get("therapeutic_count"), 0.0) or 0)
        adverse = int(safe_float(row.get("adverse_count"), 0.0) or 0)
        irrelevant = int(safe_float(row.get("irrelevant_count"), 0.0) or 0)
        article_count = therapeutic + adverse + irrelevant

    if article_count > 0:
        p_raw = therapeutic / article_count
        p_penalised = max((therapeutic - adverse_weight * adverse) / article_count, 0.0)
        adverse_rate = adverse / article_count
        irrelevant_rate = irrelevant / article_count
    else:
        p_raw = 0.5
        p_penalised = 0.5
        adverse_rate = 0.0
        irrelevant_rate = 0.0

    gamma = safe_float(comp.get("gamma") if isinstance(comp, dict) else None, None)
    if gamma is None:
        gamma = safe_float(row.get("safety_overlap_gamma"), None)

    if isinstance(comp, dict) and "safety_relation" in comp:
        safety_relation = bool(comp.get("safety_relation"))
    else:
        ledger_prior = safe_float(row.get("prior_mean"), None)
        safety_relation = bool(
            gamma is not None
            and ledger_prior is not None
            and ledger_prior < p_penalised - 1e-6
        )

    if gamma is None or not safety_relation:
        p_final = p_penalised
    else:
        p_final = p_penalised * (1.0 - safety_coefficient * gamma)
    p_final = clamp_prob(p_final)

    c_m = concentration_c(article_count, cmax=cmax, tau=tau)
    prior_a, prior_b = beta_params_from_prob(p_final, c_m)

    p_likelihood = safe_float(comp.get("p_likelihood") if isinstance(comp, dict) else None, None)
    graph_score = safe_float(comp.get("graph_score") if isinstance(comp, dict) else None, None)
    if p_likelihood is None or graph_score is None:
        p_likelihood, graph_score = graph_probability(row, config)

    like_a, like_b = beta_params_from_prob(p_likelihood, likelihood_strength)
    post_a = prior_a + (like_a - 1.0)
    post_b = prior_b + (like_b - 1.0)
    post_mean = post_a / (post_a + post_b)
    ci_low = float(beta_dist.ppf(0.025, post_a, post_b))
    ci_high = float(beta_dist.ppf(0.975, post_a, post_b))
    ci_width = ci_high - ci_low

    baseline_ci_width = safe_float(row.get("credible_interval_width"), ci_width) or ci_width
    readiness = adjusted_readiness_score(row, ci_width, baseline_ci_width)
    entity_score = min(
        safe_float(row.get("drug_mapping_score"), 0.5) or 0.5,
        safe_float(row.get("disease_mapping_score"), 0.5) or 0.5,
    )
    flag = quality_flag(
        readiness_score=readiness,
        entity_mapping_score=entity_score,
        literature_count=article_count,
        adverse_rate=adverse_rate,
        irrelevant_rate=irrelevant_rate,
        safety_gamma=gamma,
        structural_score=safe_float(row.get("structural_consistency_score"), None),
        credible_interval_width=ci_width,
    )

    return {
        "p_raw": p_raw,
        "p_penalised": p_penalised,
        "p_final": p_final,
        "gamma": gamma,
        "safety_relation": safety_relation,
        "article_count": article_count,
        "graph_score": graph_score,
        "p_likelihood": p_likelihood,
        "cM": c_m,
        "post_mean": post_mean,
        "credible_interval_width": ci_width,
        "evidence_readiness_score": readiness,
        "quality_flag": flag,
    }


def scenario_grid(config: Dict[str, Any]) -> List[Dict[str, Any]]:
    base = {
        "cmax": config["cmax"],
        "tau": config["tau"],
        "likelihood_strength": config["likelihood_strength"],
        "adverse_weight": config["adverse_weight"],
        "safety_coefficient": config["safety_coefficient"],
    }
    scenarios = [
        {
            "scenario": "baseline",
            "parameter": "baseline",
            "multiplier": 1.0,
            "parameter_value": "",
            **base,
        }
    ]
    labels = {
        "cmax": "cmax",
        "tau": "tau",
        "likelihood_strength": "lambda_likelihood_strength",
        "adverse_weight": "adverse_evidence_weight",
        "safety_coefficient": "safety_coefficient",
    }
    for parameter, label in labels.items():
        for multiplier in (0.5, 1.5):
            values = dict(base)
            values[parameter] = base[parameter] * multiplier
            scenarios.append(
                {
                    "scenario": f"{label}_{multiplier:.1f}x",
                    "parameter": label,
                    "multiplier": multiplier,
                    "parameter_value": values[parameter],
                    **values,
                }
            )
    return scenarios


def rank_pairs(df: pd.DataFrame, value_column: str, rank_column: str) -> pd.DataFrame:
    out = df.copy()
    out[rank_column] = out[value_column].rank(
        method="min",
        ascending=False,
    ).astype(int)
    return out


def build_pair_level_table(
    ledger: pd.DataFrame,
    config: Dict[str, Any],
    run_components: Dict[tuple[str, str], Dict[str, Any]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for scenario in scenario_grid(config):
        scenario_rows = []
        for _, pair in ledger.iterrows():
            run_record = run_components.get(pair_key(pair.get("drug"), pair.get("disease")))
            result = recompute_pair(
                pair,
                config,
                run_record,
                cmax=scenario["cmax"],
                tau=scenario["tau"],
                likelihood_strength=scenario["likelihood_strength"],
                adverse_weight=scenario["adverse_weight"],
                safety_coefficient=scenario["safety_coefficient"],
            )
            scenario_rows.append(
                {
                    "scenario": scenario["scenario"],
                    "parameter": scenario["parameter"],
                    "multiplier": scenario["multiplier"],
                    "parameter_value": scenario["parameter_value"],
                    "drug": pair.get("drug"),
                    "disease": pair.get("disease"),
                    "case_type": (run_record or {}).get("case_type", ""),
                    "baseline_quality_flag_from_ledger": pair.get("quality_flag"),
                    "baseline_posterior_mean_from_ledger": safe_float(
                        pair.get("posterior_mean"),
                        None,
                    ),
                    **result,
                }
            )
        scenario_df = pd.DataFrame(scenario_rows)
        scenario_df = rank_pairs(scenario_df, "post_mean", "rank")
        rows.extend(scenario_df.to_dict("records"))

    all_rows = pd.DataFrame(rows)
    baseline = (
        all_rows[all_rows["scenario"] == "baseline"][
            ["drug", "disease", "post_mean", "rank", "quality_flag"]
        ]
        .rename(
            columns={
                "post_mean": "baseline_post_mean",
                "rank": "baseline_rank",
                "quality_flag": "baseline_quality_flag",
            }
        )
    )
    out = all_rows.merge(baseline, on=["drug", "disease"], how="left")
    out["posterior_mean_delta"] = out["post_mean"] - out["baseline_post_mean"]
    out["abs_posterior_mean_delta"] = out["posterior_mean_delta"].abs()
    out["rank_delta"] = out["rank"] - out["baseline_rank"]
    out["abs_rank_delta"] = out["rank_delta"].abs()
    out["quality_flag_changed"] = out["quality_flag"] != out["baseline_quality_flag"]
    return out.sort_values(["scenario", "rank", "drug", "disease"])


def build_summary_table(pair_level: pd.DataFrame) -> pd.DataFrame:
    baseline = pair_level[pair_level["scenario"] == "baseline"][
        ["drug", "disease", "baseline_post_mean", "baseline_rank", "baseline_quality_flag"]
    ].copy()
    baseline_top10 = set(
        baseline.nsmallest(min(10, len(baseline)), "baseline_rank")
        .apply(lambda row: (row["drug"], row["disease"]), axis=1)
        .tolist()
    )

    rows = []
    for scenario, df in pair_level.groupby("scenario", sort=False):
        if scenario == "baseline":
            continue
        scenario_top10 = set(
            df.nsmallest(min(10, len(df)), "rank")
            .apply(lambda row: (row["drug"], row["disease"]), axis=1)
            .tolist()
        )
        spearman = df["baseline_post_mean"].corr(df["post_mean"], method="spearman")
        n_pairs = len(df)
        flag_changed = int(df["quality_flag_changed"].sum())
        stable_flag_pct = 100.0 * (1.0 - flag_changed / max(n_pairs, 1))
        top10_overlap = len(baseline_top10 & scenario_top10)
        top10_overlap_pct = 100.0 * top10_overlap / max(len(baseline_top10), 1)
        stable = (
            float(spearman) >= 0.95
            and float(df["abs_posterior_mean_delta"].median()) <= 0.05
            and stable_flag_pct >= 90.0
        )
        rows.append(
            {
                "scenario": scenario,
                "parameter": df["parameter"].iloc[0],
                "multiplier": df["multiplier"].iloc[0],
                "parameter_value": df["parameter_value"].iloc[0],
                "n_pairs": n_pairs,
                "spearman_rank_correlation": round(float(spearman), 4),
                "top10_overlap_count": top10_overlap,
                "top10_overlap_percent": round(top10_overlap_pct, 1),
                "median_abs_rank_change": round(float(df["abs_rank_delta"].median()), 2),
                "max_abs_rank_change": int(df["abs_rank_delta"].max()),
                "median_abs_posterior_mean_delta": round(
                    float(df["abs_posterior_mean_delta"].median()),
                    4,
                ),
                "max_abs_posterior_mean_delta": round(
                    float(df["abs_posterior_mean_delta"].max()),
                    4,
                ),
                "quality_flag_changed_pairs": flag_changed,
                "quality_flag_stable_percent": round(stable_flag_pct, 1),
                "stable_by_predefined_rule": "yes" if stable else "no",
            }
        )
    return pd.DataFrame(rows)


def write_methods_note(output_dir: Path, config: Dict[str, Any], summary: pd.DataFrame) -> None:
    stable_count = int((summary["stable_by_predefined_rule"] == "yes").sum())
    total = len(summary)
    lines = [
        "# Supplementary Sensitivity Analysis",
        "",
        "A one-at-a-time perturbation analysis was run offline from the retained "
        "publication ledger. Upstream evidence was held fixed: ClinicalTrials.gov "
        "records, MeSH mappings, PubMed/PMC classifications, safety-overlap gamma, "
        "and graph features were not regenerated.",
        "",
        "Baseline constants:",
        f"- cmax = {config['cmax']}",
        f"- tau = {config['tau']}",
        f"- lambda / likelihood_strength = {config['likelihood_strength']}",
        f"- adverse-evidence weight = {config['adverse_weight']}",
        f"- safety coefficient = {config['safety_coefficient']}",
        "",
        "Each constant was varied by 50% below and 50% above the baseline value, "
        "while all other constants were held at baseline. Stability was summarised "
        "using Spearman rank correlation, top-10 overlap, absolute rank changes, "
        "absolute posterior-mean shifts, and quality-flag changes.",
        "",
        "A scenario was labelled stable when Spearman rank correlation was at least "
        "0.95, median absolute posterior-mean shift was at most 0.05, and at least "
        "90% of quality flags were unchanged.",
        "",
        f"Stable scenarios by this rule: {stable_count}/{total}.",
        "",
        "Primary files:",
        "- SuppTable_sensitivity_parameter_summary.csv",
        "- SuppTable_sensitivity_pair_level.csv",
        "",
    ]
    (output_dir / "SuppText_sensitivity_analysis.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def generate_sensitivity_tables(run_dir: Path, output_dir: Optional[Path] = None) -> Dict[str, Path]:
    output = output_dir or run_dir / "supplementary_tables"
    output.mkdir(parents=True, exist_ok=True)

    ledger_path = run_dir / "ledgers" / "full_evidence_quality_ledger.csv"
    ledger = pd.read_csv(ledger_path)
    config = load_run_config(run_dir)
    run_components = load_run_components(run_dir)

    pair_level = build_pair_level_table(ledger, config, run_components)
    summary = build_summary_table(pair_level)

    summary_path = output / "SuppTable_sensitivity_parameter_summary.csv"
    pair_path = output / "SuppTable_sensitivity_pair_level.csv"
    pair_level.to_csv(pair_path, index=False)
    summary.to_csv(summary_path, index=False)
    write_methods_note(output, config, summary)
    return {
        "summary": summary_path,
        "pair_level": pair_path,
        "note": output / "SuppText_sensitivity_analysis.md",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate supplementary one-at-a-time sensitivity tables.",
    )
    parser.add_argument("--run_dir", default="outputs/20260610_bayesian")
    parser.add_argument("--output_dir", default=None)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    paths = generate_sensitivity_tables(
        run_dir=Path(args.run_dir),
        output_dir=Path(args.output_dir) if args.output_dir else None,
    )
    print("Generated sensitivity supplementary files:")
    for label, path in paths.items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
