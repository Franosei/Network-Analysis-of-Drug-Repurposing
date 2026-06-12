"""
sensitivity_cmax_tau.py

Prior → posterior sensitivity analysis for (cmax, tau) in the concentration function:

    c(M) = cmax * (1 - exp(-M / tau))

This script runs NO LLM calls and NO PubMed calls.
It treats semantic outputs as fixed constants (from a validated run):
  - M
  - p_final (literature + safety adjusted prior mean)

It also treats likelihood as fixed constants (graph evidence):
  - p_likelihood
  - likelihood_strength

Outputs:
  - CSV of the full grid with prior/posterior summaries
  - Heatmaps (publication-ready) for posterior mean, KL, delta_mu, and c(M)
"""

import os
import math
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid


# ---------------------------------------------------------------------
# Fixed inputs from your validated run (NO LLM calls here)
# ---------------------------------------------------------------------
@dataclass(frozen=True)
class FixedRunInputs:
    # Semantic outputs (fixed)
    M: int = 150
    p_final: float = 0.203167  # enhanced prior mean after safety penalty
    # Likelihood outputs (fixed, static graph)
    p_likelihood: float = 0.47849860958775603
    likelihood_strength: float = 50.0

    # Optional: keep for reporting (not used in computation)
    p_raw: float = 0.646667
    p_penalised: float = 0.353333
    gamma: float = 0.85
    counts: Dict[str, int] = None


# ---------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------
def clamp01(x: Any, eps: float = 1e-6) -> float:
    """
    Clamp a value into (0,1) with small epsilon padding for Beta stability.
    """
    try:
        v = float(x)
    except Exception:
        v = 0.5
    v = max(min(v, 1.0), 0.0)
    v = min(max(v, eps), 1.0 - eps)
    return v


def concentration_c(M: int, cmax: float, tau: float) -> float:
    """
    c(M) = cmax * (1 - exp(-M/tau))
    """
    M = int(max(M, 0))
    cmax = float(max(cmax, 0.0))
    tau = float(max(tau, 1e-6))
    return cmax * (1.0 - math.exp(-M / tau))


def beta_params_from_prob(p: float, strength: float) -> Tuple[float, float]:
    """
    α = 1 + strength * p
    β = 1 + strength * (1 - p)
    """
    p = clamp01(p)
    strength = max(float(strength), 0.0)
    a = 1.0 + strength * p
    b = 1.0 + strength * (1.0 - p)
    return a, b


def beta_mean(a: float, b: float) -> float:
    return float(a / (a + b))


def compute_kl_and_delta(prior_a: float, prior_b: float, post_a: float, post_b: float, n_grid: int = 3000) -> Tuple[float, float, float]:
    """
    Compute KL(post || prior), prior mean, posterior mean, delta_mu.

    We approximate KL using numerical integration over (0,1) because it is robust,
    consistent with your existing plotting/diagnostics approach.
    """
    x = np.linspace(0.001, 0.999, int(max(n_grid, 500)))
    prior_pdf = beta.pdf(x, a=prior_a, b=prior_b)
    post_pdf = beta.pdf(x, a=post_a, b=post_b)

    prior_pdf = np.clip(prior_pdf, 1e-12, None)
    post_pdf = np.clip(post_pdf, 1e-12, None)

    prior_pdf /= trapezoid(prior_pdf, x)
    post_pdf /= trapezoid(post_pdf, x)

    kl = trapezoid(rel_entr(post_pdf, prior_pdf), x)
    mu_prior = trapezoid(x * prior_pdf, x)
    mu_post = trapezoid(x * post_pdf, x)
    delta = mu_post - mu_prior

    return float(kl), float(mu_prior), float(mu_post), float(delta)


# ---------------------------------------------------------------------
# Sensitivity analysis
# ---------------------------------------------------------------------
def run_sensitivity(
    fixed: FixedRunInputs,
    cmax_values: List[float],
    tau_values: List[float],
) -> pd.DataFrame:
    """
    Compute prior → posterior outputs for each (cmax, tau) pair.
    """
    M = int(fixed.M)
    p_final = clamp01(fixed.p_final)
    p_like = clamp01(fixed.p_likelihood)
    like_strength = float(fixed.likelihood_strength)

    # Likelihood Beta is fixed
    like_a, like_b = beta_params_from_prob(p_like, like_strength)

    rows = []
    for cmax in cmax_values:
        for tau in tau_values:
            cM = concentration_c(M, cmax=cmax, tau=tau)
            prior_a, prior_b = beta_params_from_prob(p_final, strength=cM)

            # Pseudo-count fusion (subtract baseline +1)
            post_a = prior_a + (like_a - 1.0)
            post_b = prior_b + (like_b - 1.0)

            kl, mu_prior, mu_post, delta = compute_kl_and_delta(prior_a, prior_b, post_a, post_b)

            rows.append(
                {
                    "cmax": float(cmax),
                    "tau": float(tau),
                    "M": M,
                    "cM": float(cM),
                    "p_final": float(p_final),
                    "p_likelihood": float(p_like),
                    "likelihood_strength": float(like_strength),
                    "prior_a": float(prior_a),
                    "prior_b": float(prior_b),
                    "prior_mean": float(mu_prior),
                    "post_a": float(post_a),
                    "post_b": float(post_b),
                    "post_mean": float(mu_post),
                    "delta_mu": float(delta),
                    "KL_post_prior": float(kl),
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Plotting (publication-friendly, no clutter)
# ---------------------------------------------------------------------
def _heatmap(
    df: pd.DataFrame,
    value_col: str,
    cmax_values: List[float],
    tau_values: List[float],
    title: str,
    outpath: str,
    xlabel: str = r"$\tau$",
    ylabel: str = r"$c_{\max}$",
) -> None:
    """
    Create a heatmap with rows=cmax, cols=tau.
    """
    pivot = df.pivot(index="cmax", columns="tau", values=value_col).reindex(index=cmax_values, columns=tau_values)

    plt.figure(figsize=(9, 5.5))
    im = plt.imshow(pivot.values, aspect="auto", origin="lower")

    plt.xticks(range(len(tau_values)), [str(t) for t in tau_values])
    plt.yticks(range(len(cmax_values)), [str(c) for c in cmax_values])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    cbar = plt.colorbar(im)
    cbar.set_label(value_col)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def save_outputs(df: pd.DataFrame, out_dir: str, tag: str, cmax_values: List[float], tau_values: List[float]) -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)

    csv_path = os.path.join(out_dir, f"sensitivity_cmax_tau_{tag}.csv")
    df.to_csv(csv_path, index=False)

    # Heatmaps
    fig_post = os.path.join(out_dir, f"heatmap_post_mean_{tag}.png")
    _heatmap(
        df,
        value_col="post_mean",
        cmax_values=cmax_values,
        tau_values=tau_values,
        title="Posterior Mean Sensitivity to $c_{max}$ and $\\tau$",
        outpath=fig_post,
    )

    fig_kl = os.path.join(out_dir, f"heatmap_KL_{tag}.png")
    _heatmap(
        df,
        value_col="KL_post_prior",
        cmax_values=cmax_values,
        tau_values=tau_values,
        title="KL Divergence Sensitivity (KL(post || prior))",
        outpath=fig_kl,
    )

    fig_delta = os.path.join(out_dir, f"heatmap_delta_mu_{tag}.png")
    _heatmap(
        df,
        value_col="delta_mu",
        cmax_values=cmax_values,
        tau_values=tau_values,
        title="Mean-Shift Sensitivity ($\\Delta \\mu = E[post]-E[prior]$)",
        outpath=fig_delta,
    )

    fig_cM = os.path.join(out_dir, f"heatmap_cM_{tag}.png")
    _heatmap(
        df,
        value_col="cM",
        cmax_values=cmax_values,
        tau_values=tau_values,
        title="Prior Concentration $c(M)$ Across Hyperparameters",
        outpath=fig_cM,
    )

    return {
        "csv": csv_path,
        "heatmap_post_mean": fig_post,
        "heatmap_KL": fig_kl,
        "heatmap_delta_mu": fig_delta,
        "heatmap_cM": fig_cM,
    }


def summarize_grid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compact stats to report in text / supplement.
    """
    cols = ["cM", "prior_mean", "post_mean", "delta_mu", "KL_post_prior"]
    summary_rows = []
    for c in cols:
        s = df[c].dropna()
        summary_rows.append(
            {
                "metric": c,
                "mean": float(s.mean()),
                "sd": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
                "min": float(s.min()),
                "max": float(s.max()),
            }
        )
    return pd.DataFrame(summary_rows)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
if __name__ == "__main__":
    # Fixed constants from your validated run
    fixed = FixedRunInputs(
        M=150,
        p_final=0.203167,
        p_likelihood=0.47849860958775603,
        likelihood_strength=50.0,
        p_raw=0.646667,
        p_penalised=0.353333,
        gamma=0.85,
        counts={"therapeutic": 97, "irrelevant": 31, "adverse": 22},
    )

    # 5×5 grid (centered around your manuscript defaults cmax=200, tau=25)
    cmax_values = [50, 100, 200, 300, 400]
    tau_values = [10, 15, 25, 40, 60]

    df = run_sensitivity(fixed, cmax_values=cmax_values, tau_values=tau_values)

    tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = "sensitivity_cmax_tau_outputs"

    paths = save_outputs(df, out_dir=out_dir, tag=tag, cmax_values=cmax_values, tau_values=tau_values)

    summ = summarize_grid(df)
    summ_path = os.path.join(out_dir, f"sensitivity_summary_{tag}.csv")
    summ.to_csv(summ_path, index=False)

    print("\n" + "=" * 80)
    print("SENSITIVITY ANALYSIS COMPLETE (no LLM calls)")
    print("=" * 80)
    print(f"Saved CSV grid: {paths['csv']}")
    print(f"Saved summary:  {summ_path}")
    print(f"Saved heatmap (post mean): {paths['heatmap_post_mean']}")
    print(f"Saved heatmap (KL):        {paths['heatmap_KL']}")
    print(f"Saved heatmap (delta_mu):  {paths['heatmap_delta_mu']}")
    print(f"Saved heatmap (cM):        {paths['heatmap_cM']}")
    print("\nSummary stats:")
    print(summ.to_string(index=False))
