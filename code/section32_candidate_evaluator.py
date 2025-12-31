import argparse
import pandas as pd
import numpy as np
from scipy.special import expit  # sigmoid


DEFAULT_WEIGHTS = {
    # You can edit these if you have the finalized weights elsewhere.
    # NOTE: Using GraphDistance with a negative weight assumes smaller distance => higher likelihood.
    "GraphDistanceToIndication": -23.32,
    "RandomWalkScore": 0.25,
    "StructuralLikelihood": 0.11,
    "PreferentialAttachment": -0.04,
    "KatzSimilarity": 1.65,
}


def ensure_columns(df: pd.DataFrame, required: list[str], df_name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"{df_name} is missing columns: {missing}. Available: {list(df.columns)}"
        )


def weighted_graph_score(row: pd.Series, weights: dict[str, float]) -> float:
    return float(sum(weights[f] * row[f] for f in weights.keys()))


def graph_probability(score: float) -> float:
    return float(expit(score))


def concentration_c(M: float, cmax: float = 200.0, tau: float = 25.0) -> float:
    """
    Evidence-based prior concentration. Keep consistent with your manuscript.
    """
    return float(cmax * (1.0 - np.exp(-M / tau)))


def beta_prior_params(p_final: float, M: float) -> tuple[float, float]:
    """
    Construct Beta prior parameters from final prior mean and literature count M.
    """
    cM = concentration_c(M)
    alpha = 1.0 + cM * p_final
    beta = 1.0 + cM * (1.0 - p_final)
    return alpha, beta


def posterior_mean_from_soft_likelihood(p_prior: float, M: float, p_like: float) -> float:
    """
    Soft update: combine Beta prior with a Bernoulli-like soft likelihood probability.

    This is a pragmatic, monotone fusion that preserves:
      - more articles => stronger prior
      - stronger p_like => larger posterior shift upward (and vice versa)
    """
    alpha, beta = beta_prior_params(p_prior, M)
    return float((alpha * p_like) / (alpha * p_like + beta * (1.0 - p_like)))


def shortlist_for_drug(
    df_unknown: pd.DataFrame,
    drug: str,
    drug_col: str,
    disease_col: str,
    weights: dict[str, float],
    distance_col: str,
    distance_max: float | None,
    min_p_like: float,
    top_n: int,
) -> pd.DataFrame:
    """
    Shortlist structurally-informative non-edge candidates for a given drug.
    """
    df = df_unknown[df_unknown[drug_col].astype(str).str.lower() == drug.lower()].copy()

    if df.empty:
        return df

    # Optional distance filter
    if distance_max is not None:
        df = df[df[distance_col] <= distance_max].copy()

    if df.empty:
        return df

    df["graph_score"] = df.apply(weighted_graph_score, axis=1, weights=weights)
    df["p_like"] = df["graph_score"].map(graph_probability)

    # Filter to avoid the non-informative regime
    df = df[df["p_like"] >= min_p_like].copy()

    df.sort_values("p_like", ascending=False, inplace=True)
    return df.head(top_n)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Section 3.2: shortlist structurally-informative non-edge drug–disease candidates."
    )

    parser.add_argument("--unknown", required=True, help="Path to graph_features_unknown.csv")
    parser.add_argument("--drug", action="append", required=True, help="Drug name (repeatable)")
    parser.add_argument("--out", default="section32_shortlist.csv", help="Output CSV path")

    # Column mapping (matches your CSV)
    parser.add_argument("--drug-col", default="Drug")
    parser.add_argument("--disease-col", default="Disease")
    parser.add_argument("--distance-col", default="GraphDistanceToIndication")

    # Selection controls
    parser.add_argument("--distance-max", type=float, default=3.0,
                        help="Keep candidates with GraphDistanceToIndication <= this. Use a float; set to -1 to disable.")
    parser.add_argument("--min-p-like", type=float, default=0.55,
                        help="Minimum structural likelihood probability for inclusion.")
    parser.add_argument("--candidates", type=int, default=30,
                        help="Max shortlisted candidates per drug.")

    # Optional priors merge (for running posterior only on shortlisted pairs)
    parser.add_argument("--priors", default=None,
                        help="Optional CSV with columns: Drug, Disease, p_prior_final, pubmed_count")
    parser.add_argument("--prior-col", default="p_prior_final")
    parser.add_argument("--pubmed-col", default="pubmed_count")
    parser.add_argument("--min-pubmed", type=float, default=0.0,
                        help="If priors provided, require pubmed_count >= this threshold before computing posterior.")

    args = parser.parse_args()

    df_unknown = pd.read_csv(args.unknown)

    required_unknown = [
        args.drug_col, args.disease_col, args.distance_col,
        "RandomWalkScore", "StructuralLikelihood", "PreferentialAttachment", "KatzSimilarity"
    ]
    ensure_columns(df_unknown, required_unknown, "Unknown feature file")

    weights = DEFAULT_WEIGHTS.copy()

    # Distance max toggle
    distance_max = None if args.distance_max is not None and args.distance_max < 0 else args.distance_max

    all_rows = []
    for d in args.drug:
        df_s = shortlist_for_drug(
            df_unknown=df_unknown,
            drug=d,
            drug_col=args.drug_col,
            disease_col=args.disease_col,
            weights=weights,
            distance_col=args.distance_col,
            distance_max=distance_max,
            min_p_like=args.min_p_like,
            top_n=args.candidates,
        )
        if df_s.empty:
            continue

        df_s["SelectedDrug"] = d
        all_rows.append(df_s)

    if not all_rows:
        print("[WARN] No candidates matched your filters. Try increasing --distance-max or lowering --min-p-like.")
        pd.DataFrame().to_csv(args.out, index=False)
        return

    df_out = pd.concat(all_rows, ignore_index=True)

    # Optional: merge priors and compute posterior on shortlist only
    if args.priors:
        df_priors = pd.read_csv(args.priors)
        ensure_columns(df_priors, [args.drug_col, args.disease_col, args.prior_col, args.pubmed_col], "Priors file")

        df_out = df_out.merge(
            df_priors[[args.drug_col, args.disease_col, args.prior_col, args.pubmed_col]],
            on=[args.drug_col, args.disease_col],
            how="left",
        )

        # Compute posterior only where prior info exists and passes pubmed threshold
        has_prior = df_out[args.prior_col].notna() & df_out[args.pubmed_col].notna()
        enough_pubmed = has_prior & (df_out[args.pubmed_col] >= args.min_pubmed)

        df_out["posterior_mean"] = np.nan
        df_out["delta_mu"] = np.nan

        idx = df_out.index[enough_pubmed]
        for i in idx:
            p_prior = float(df_out.loc[i, args.prior_col])
            M = float(df_out.loc[i, args.pubmed_col])
            p_like = float(df_out.loc[i, "p_like"])
            post = posterior_mean_from_soft_likelihood(p_prior, M, p_like)
            df_out.loc[i, "posterior_mean"] = post
            df_out.loc[i, "delta_mu"] = post - p_prior

        # If posterior exists, sort by delta_mu then p_like
        if df_out["delta_mu"].notna().any():
            df_out.sort_values(["delta_mu", "p_like"], ascending=[False, False], inplace=True)
        else:
            df_out.sort_values("p_like", ascending=False, inplace=True)
    else:
        df_out.sort_values("p_like", ascending=False, inplace=True)

    df_out.to_csv(args.out, index=False)
    print(f"[OK] Wrote shortlist to {args.out}")


if __name__ == "__main__":
    main()
