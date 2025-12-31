"""
bayesian_predictor_stability.py

Extension of  existing bayesian_predictor.py that ADDS:
  - 10-run repeatability evaluation for a fixed drug–disease pair
  - Publication-ready stability plots:
      (1) Prior stability: p_raw, p_penalised, p_final, gamma vs run index
      (2) Posterior stability: posterior mean vs run index
      (3) KL stability: KL divergence vs run index
  - Saves:
      - per-run JSON logs (existing behavior)
      - a single CSV summary for the 10 runs
      - the three stability figures

NOTE:
- Core inference code path (priors/likelihood/posterior) is unchanged.
- This only adds a stability routine and plotting/saving.
"""

import os
import json
import math
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid

from pubmed_utils import LLMClassifier, PubMedSearchConfig, LLMConfig


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


def sigmoid(z: float) -> float:
    """
    Numerically stable sigmoid.
    """
    if z >= 0:
        ez = math.exp(-z)
        return 1.0 / (1.0 + ez)
    ez = math.exp(z)
    return ez / (1.0 + ez)


def beta_params_from_prob(p: float, strength: float) -> Tuple[float, float]:
    """
    Convert probability p and pseudo-count strength into Beta(alpha, beta).
    α = 1 + strength*p
    β = 1 + strength*(1-p)
    """
    p = clamp01(p)
    strength = max(float(strength), 0.0)
    alpha = 1.0 + strength * p
    beta_ = 1.0 + strength * (1.0 - p)
    return alpha, beta_


def concentration_c(M: int, cmax: float = 200.0, tau: float = 25.0) -> float:
    """
    Evidence-scaled concentration function c(M), bounded and increasing in M.

    Default choice (smooth saturation):
      c(M) = cmax * (1 - exp(-M / tau))
    """
    M = int(max(M, 0))
    cmax = float(max(cmax, 0.0))
    tau = float(max(tau, 1e-6))
    return cmax * (1.0 - math.exp(-M / tau))


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
@dataclass(frozen=True)
class PredictorConfig:
    # PubMed/LLM retrieval
    pubmed_max_articles: int = 150
    pubmed_filter_level: str = "high"
    pubmed_years_back: int = 10
    pubmed_use_cache: bool = True

    llm_model: str = "gpt-4o-mini"
    llm_batch_size: int = 5
    llm_delay_s: float = 2.0
    llm_max_retries: int = 2

    # Prior concentration scaling (must match manuscript)
    cmax: float = 200.0
    tau: float = 25.0

    # Likelihood (graph evidence) pseudo-count strength
    likelihood_strength: float = 50.0

    # Likelihood score -> probability mapping
    likelihood_intercept: float = 0.0

    # Output & reproducibility
    plots_dir: str = "plots"
    logs_dir: str = "runs"
    save_run_logs: bool = True

    # Stability evaluation (NEW)
    stability_dir: str = "stability"
    stability_n_runs: int = 10


# ---------------------------------------------------------------------
# PubMed semantic prior wrapper (aligned with pubmed_utils.py)
# ---------------------------------------------------------------------
def make_classifier(cfg: PredictorConfig) -> LLMClassifier:
    search_cfg = PubMedSearchConfig(
        max_results=int(cfg.pubmed_max_articles),
        years_back=int(cfg.pubmed_years_back),
    )
    llm_cfg = LLMConfig(
        model=str(cfg.llm_model),
        delay_s=float(cfg.llm_delay_s),
        batch_size=int(cfg.llm_batch_size),
        max_retries=int(cfg.llm_max_retries),
    )
    return LLMClassifier(search_cfg=search_cfg, llm_cfg=llm_cfg)


# ---------------------------------------------------------------------
# Predictor
# ---------------------------------------------------------------------
class BayesianRepurposingPredictor:
    """
    Bayesian engine for drug repurposing using:
      - PubMed semantic priors (Praw, Ppen, Pfinal) + evidence scaling c(M)
      - Graph feature evidence mapped to probability via sigmoid
      - Beta pseudo-count fusion for posterior
      - Transparent printing and optional JSON logging for reproducibility
    """

    def __init__(self, known_path: str, unknown_path: str, weights_dict: Dict[str, float], cfg: Optional[PredictorConfig] = None):
        self.cfg = cfg or PredictorConfig()
        self.feature_df = self.load_combined_features(known_path, unknown_path)
        self.weights_dict = dict(weights_dict)
        self.existing_pairs = set(
            zip(
                self.feature_df[self.feature_df["Label"] == 1]["Drug"],
                self.feature_df[self.feature_df["Label"] == 1]["Disease"],
            )
        )
        self.classifier = make_classifier(self.cfg)

        os.makedirs(self.cfg.plots_dir, exist_ok=True)
        os.makedirs(self.cfg.logs_dir, exist_ok=True)
        os.makedirs(self.cfg.stability_dir, exist_ok=True)

    @staticmethod
    def load_combined_features(known_csv: str, unknown_csv: str) -> pd.DataFrame:
        df_known = pd.read_csv(known_csv)
        df_unknown = pd.read_csv(unknown_csv)
        df = pd.concat([df_known, df_unknown], ignore_index=True)
        df["Drug"] = df["Drug"].astype(str).str.strip().str.lower()
        df["Disease"] = df["Disease"].astype(str).str.strip().str.lower()
        return df

    def _pubmed_prior(self, drug: str, disease: str) -> Dict[str, Any]:
        return self.classifier.build_semantic_prior(
            drug=drug,
            disease=disease,
            max_articles=int(self.cfg.pubmed_max_articles),
            filter_level=str(self.cfg.pubmed_filter_level),
            use_cache=bool(self.cfg.pubmed_use_cache),
        )

    def _graph_score(self, drug: str, disease: str) -> Tuple[float, Dict[str, float]]:
        drug_l, disease_l = drug.lower(), disease.lower()
        row = self.feature_df[(self.feature_df["Drug"] == drug_l) & (self.feature_df["Disease"] == disease_l)]
        if row.empty:
            return float(self.cfg.likelihood_intercept), {}

        feature_vals: Dict[str, float] = {}
        s = float(self.cfg.likelihood_intercept)
        r0 = row.iloc[0]
        for feature, weight in self.weights_dict.items():
            val = float(r0.get(feature, 0.0))
            feature_vals[feature] = val
            s += float(weight) * val
        return s, feature_vals

    def _graph_probability(self, score: float) -> float:
        return clamp01(sigmoid(score))

    def compute_components(self, drug: str, disease: str) -> Dict[str, Any]:
        prior_result: Dict[str, Any] = {}
        try:
            prior_result = self._pubmed_prior(drug, disease)
        except Exception as e:
            prior_result = {
                "prior": 0.5,
                "penalised_prior": 0.5,
                "enhanced_prior": 0.5,
                "gamma": None,
                "raw_counts": {},
                "total_articles": 0,
            }
            print(f"[WARN] PubMed prior failed for {drug} → {disease}: {e}")

        p_raw = clamp01(prior_result.get("prior", 0.5))
        p_pen = clamp01(prior_result.get("penalised_prior", 0.5))
        p_final = clamp01(prior_result.get("enhanced_prior", 0.5))
        gamma = prior_result.get("gamma", None)
        gamma = clamp01(gamma) if gamma is not None else None

        counts = prior_result.get("raw_counts", {})
        try:
            counts = dict(counts)
        except Exception:
            counts = {}

        M = int(prior_result.get("total_articles", 0))

        cM = concentration_c(M, cmax=self.cfg.cmax, tau=self.cfg.tau)
        prior_a, prior_b = beta_params_from_prob(p_final, strength=cM)

        score, feature_vals = self._graph_score(drug, disease)
        p_like = self._graph_probability(score)

        like_a, like_b = beta_params_from_prob(p_like, strength=self.cfg.likelihood_strength)

        post_a = prior_a + (like_a - 1.0)
        post_b = prior_b + (like_b - 1.0)
        post_mean = post_a / (post_a + post_b)

        return {
            "p_raw": p_raw,
            "p_penalised": p_pen,
            "p_final": p_final,
            "gamma": gamma,
            "counts": counts,
            "M": M,
            "cM": cM,
            "prior_a": prior_a,
            "prior_b": prior_b,
            "graph_score": score,
            "p_likelihood": p_like,
            "likelihood_strength": float(self.cfg.likelihood_strength),
            "like_a": like_a,
            "like_b": like_b,
            "feature_values": feature_vals,
            "post_a": post_a,
            "post_b": post_b,
            "post_mean": float(post_mean),
        }

    @staticmethod
    def compute_kl_and_mean_shift(prior_pdf: np.ndarray, post_pdf: np.ndarray, x: np.ndarray) -> Tuple[float, float, float, float]:
        prior_pdf = np.clip(prior_pdf, 1e-12, None)
        post_pdf = np.clip(post_pdf, 1e-12, None)

        prior_pdf /= trapezoid(prior_pdf, x)
        post_pdf /= trapezoid(post_pdf, x)

        kl = trapezoid(rel_entr(post_pdf, prior_pdf), x)
        mu_prior = trapezoid(x * prior_pdf, x)
        mu_post = trapezoid(x * post_pdf, x)
        delta = mu_post - mu_prior
        return float(round(kl, 6)), float(round(mu_prior, 6)), float(round(mu_post, 6)), float(round(delta, 6))

    def plot_distributions(self, comp: Dict[str, Any], drug: str, disease: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        x = np.linspace(0.001, 0.999, 1000)

        prior_pdf = beta.pdf(x, a=comp["prior_a"], b=comp["prior_b"])
        like_pdf = beta.pdf(x, a=comp["like_a"], b=comp["like_b"])
        post_pdf = beta.pdf(x, a=comp["post_a"], b=comp["post_b"])

        prior_ci = beta.interval(0.95, comp["prior_a"], comp["prior_b"])
        like_ci = beta.interval(0.95, comp["like_a"], comp["like_b"])
        post_ci = beta.interval(0.95, comp["post_a"], comp["post_b"])

        plt.figure(figsize=(10, 6))
        plt.plot(x, prior_pdf / prior_pdf.max(), label="Prior")
        plt.plot(x, like_pdf / like_pdf.max(), label="Likelihood")
        plt.plot(x, post_pdf / post_pdf.max(), label="Posterior")

        plt.axvspan(*prior_ci, alpha=0.1)
        plt.axvspan(*like_ci, alpha=0.1)
        plt.axvspan(*post_ci, alpha=0.1)

        plt.title(f"Bayesian Update: {drug} → {disease}")
        plt.xlabel("Association Strength")
        plt.ylabel("Normalized Density")
        plt.legend()
        plt.tight_layout()

        fname = f"{drug.strip().lower().replace(' ', '_')}_{disease.strip().lower().replace(' ', '_')}.png"
        fpath = os.path.join(self.cfg.plots_dir, fname)
        plt.savefig(fpath, dpi=300)
        plt.close()

        return x, prior_pdf, post_pdf

    def _print_diagnostics(self, drug: str, disease: str, comp: Dict[str, Any]) -> None:
        print("\n" + "-" * 70)
        print(f"PAIR: {drug} → {disease}")
        print("-" * 70)

        print("[Semantic prior / PubMed]")
        print(f"  Filter level: {self.cfg.pubmed_filter_level}")
        print(f"  Years back:   {self.cfg.pubmed_years_back}")
        print(f"  Max articles: {self.cfg.pubmed_max_articles}")
        print(f"  Use cache:    {self.cfg.pubmed_use_cache}")
        print(f"  Articles used (M): {comp['M']}")
        print(f"  Counts (raw): {comp['counts']}")
        print(f"  P_raw:        {comp['p_raw']:.6f}")
        print(f"  P_penalised:  {comp['p_penalised']:.6f}")
        print(f"  P_final:      {comp['p_final']:.6f}")
        print(f"  gamma:        {comp['gamma'] if comp['gamma'] is not None else 'None'}")
        print(f"  c(M):         {comp['cM']:.6f}")
        print(f"  Prior Beta:   alpha={comp['prior_a']:.6f}, beta={comp['prior_b']:.6f}")

        print("\n[Graph evidence / Likelihood]")
        print(f"  Intercept:             {self.cfg.likelihood_intercept:.6f}")
        print(f"  Linear score:          {comp['graph_score']:.6f}")
        print(f"  p_likelihood(sigmoid): {comp['p_likelihood']:.6f}")
        print(f"  Likelihood strength:   {comp['likelihood_strength']:.6f}")
        print(f"  Likelihood Beta:       alpha={comp['like_a']:.6f}, beta={comp['like_b']:.6f}")

        print("\n[Posterior]")
        print(f"  Posterior Beta: alpha={comp['post_a']:.6f}, beta={comp['post_b']:.6f}")
        print(f"  Posterior mean: {comp['post_mean']:.6f}")

    def _save_run_log(self, drug: str, disease: str, comp: Dict[str, Any], extra: Optional[Dict[str, Any]] = None) -> Optional[str]:
        if not self.cfg.save_run_logs:
            return None

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_drug = drug.strip().lower().replace(" ", "_").replace("/", "_")
        safe_dis = disease.strip().lower().replace(" ", "_").replace("/", "_")
        path = os.path.join(self.cfg.logs_dir, f"run_{safe_drug}_{safe_dis}_{ts}.json")

        payload = {
            "timestamp": ts,
            "drug": drug,
            "disease": disease,
            "config": {
                "pubmed_max_articles": self.cfg.pubmed_max_articles,
                "pubmed_filter_level": self.cfg.pubmed_filter_level,
                "pubmed_years_back": self.cfg.pubmed_years_back,
                "pubmed_use_cache": self.cfg.pubmed_use_cache,
                "llm_model": self.cfg.llm_model,
                "llm_batch_size": self.cfg.llm_batch_size,
                "llm_delay_s": self.cfg.llm_delay_s,
                "llm_max_retries": self.cfg.llm_max_retries,
                "cmax": self.cfg.cmax,
                "tau": self.cfg.tau,
                "likelihood_strength": self.cfg.likelihood_strength,
                "likelihood_intercept": self.cfg.likelihood_intercept,
            },
            "components": comp,
        }
        if extra:
            payload["extra"] = extra

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
            return path
        except Exception:
            return None

    # -----------------------------------------------------------------
    # NEW: Stability evaluation for a fixed drug–disease pair (n runs)
    # -----------------------------------------------------------------
    def evaluate_pair_stability(self, drug: str, disease: str, n_runs: Optional[int] = None) -> Dict[str, Any]:
        """
        Repeat the full pipeline n times with identical configuration to assess
        LLM-driven stability of priors (p_raw, p_penalised, p_final, gamma) and
        downstream posterior summaries.

        Saves:
          - per-run JSON logs (existing mechanism)
          - one CSV with all run-level metrics
          - 3 publication-ready plots in cfg.stability_dir
        """
        n = int(n_runs if n_runs is not None else self.cfg.stability_n_runs)
        n = max(1, n)

        drug_clean = drug.strip()
        disease_clean = disease.strip()

        rows: List[Dict[str, Any]] = []

        print("\n" + "=" * 80)
        print(f"STABILITY EVALUATION: {drug_clean} → {disease_clean} (n={n})")
        print("=" * 80)

        # Fixed x-grid for KL computations
        x = np.linspace(0.001, 0.999, 2000)

        for run_idx in range(1, n + 1):
            print(f"\n[RUN {run_idx}/{n}] {drug_clean} → {disease_clean}")
            comp = self.compute_components(drug_clean, disease_clean)
            self._print_diagnostics(drug_clean, disease_clean, comp)

            # KL + mean shift (prior→posterior) computed deterministically from Beta params
            prior_pdf = beta.pdf(x, a=comp["prior_a"], b=comp["prior_b"])
            post_pdf = beta.pdf(x, a=comp["post_a"], b=comp["post_b"])
            kl, mu_prior, mu_post, delta = self.compute_kl_and_mean_shift(prior_pdf, post_pdf, x)

            counts = comp.get("counts", {}) or {}
            T = int(counts.get("therapeutic", 0))
            A = int(counts.get("adverse", 0))
            N = int(counts.get("irrelevant", 0))

            row = {
                "run": run_idx,
                "p_raw": comp["p_raw"],
                "p_penalised": comp["p_penalised"],
                "p_final": comp["p_final"],
                "gamma": comp["gamma"] if comp["gamma"] is not None else np.nan,
                "T": T,
                "A": A,
                "N": N,
                "M": comp["M"],
                "cM": comp["cM"],
                "p_likelihood": comp["p_likelihood"],
                "post_mean": comp["post_mean"],
                "KL": kl,
                "E_prior": mu_prior,
                "E_post": mu_post,
                "delta_mu": delta,
            }
            rows.append(row)

            # Save per-run log (unchanged behavior)
            self._save_run_log(drug_clean, disease_clean, comp, extra={"run_index": run_idx, "KL": kl, "delta_mu": delta})

        df = pd.DataFrame(rows)

        # Save CSV summary
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_drug = drug_clean.lower().replace(" ", "_").replace("/", "_")
        safe_dis = disease_clean.lower().replace(" ", "_").replace("/", "_")
        csv_path = os.path.join(self.cfg.stability_dir, f"stability_{safe_drug}_{safe_dis}_{ts}.csv")
        df.to_csv(csv_path, index=False)

        # Plot 1: Prior stability (probabilities)
        fig1_path = os.path.join(self.cfg.stability_dir, f"stability_priors_{safe_drug}_{safe_dis}_{ts}.png")
        self._plot_prior_stability(df, drug_clean, disease_clean, fig1_path)

        # Plot 2: Posterior mean stability
        fig2_path = os.path.join(self.cfg.stability_dir, f"stability_posterior_{safe_drug}_{safe_dis}_{ts}.png")
        self._plot_posterior_stability(df, drug_clean, disease_clean, fig2_path)

        # Plot 3: KL stability
        fig3_path = os.path.join(self.cfg.stability_dir, f"stability_kl_{safe_drug}_{safe_dis}_{ts}.png")
        self._plot_kl_stability(df, drug_clean, disease_clean, fig3_path)

        # Summary stats for quick reporting
        summary = self._stability_summary(df)

        print("\n" + "-" * 80)
        print("STABILITY SUMMARY (mean ± SD | min..max)")
        print("-" * 80)
        for k, v in summary.items():
            print(f"{k:14s}: {v}")
        print("-" * 80)
        print(f"[SAVED] CSV: {csv_path}")
        print(f"[SAVED] Figure: {fig1_path}")
        print(f"[SAVED] Figure: {fig2_path}")
        print(f"[SAVED] Figure: {fig3_path}")

        return {
            "csv_path": csv_path,
            "fig_prior_path": fig1_path,
            "fig_posterior_path": fig2_path,
            "fig_kl_path": fig3_path,
            "summary": summary,
            "dataframe": df,
        }

    # ---------------------------
    # NEW: Publication-ready plots
    # ---------------------------
    @staticmethod
    def _plot_prior_stability(df: pd.DataFrame, drug: str, disease: str, outpath: str) -> None:
        """
        Prior stability plot:
          x = run index
          y = p_raw, p_penalised, p_final, gamma
        """
        plt.figure(figsize=(10, 5))
        x = df["run"].values

        plt.plot(x, df["p_raw"].values, marker="o", linewidth=1.5, label="p_raw")
        plt.plot(x, df["p_penalised"].values, marker="o", linewidth=1.5, label="p_penalised")
        plt.plot(x, df["p_final"].values, marker="o", linewidth=1.5, label="p_final")
        plt.plot(x, df["gamma"].values, marker="o", linewidth=1.5, label="gamma")

        plt.ylim(0, 1)
        plt.xlabel("Run index")
        plt.ylabel("Value (0–1)")
        plt.title(f"Prior Stability Across Repeated Runs: {drug} → {disease}")
        plt.legend(frameon=False, ncol=2)
        plt.tight_layout()
        plt.savefig(outpath, dpi=300)
        plt.close()

    @staticmethod
    def _plot_posterior_stability(df: pd.DataFrame, drug: str, disease: str, outpath: str) -> None:
        """
        Posterior mean stability plot:
          x = run index
          y = posterior mean
        """
        plt.figure(figsize=(10, 4.5))
        x = df["run"].values

        plt.plot(x, df["post_mean"].values, marker="o", linewidth=1.5, label="posterior mean")
        plt.ylim(0, 1)
        plt.xlabel("Run index")
        plt.ylabel("Posterior mean")
        plt.title(f"Posterior Stability Across Repeated Runs: {drug} → {disease}")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(outpath, dpi=300)
        plt.close()

    @staticmethod
    def _plot_kl_stability(df: pd.DataFrame, drug: str, disease: str, outpath: str) -> None:
        """
        KL divergence stability plot:
          x = run index
          y = KL(prior || posterior) measured as KL(post || prior) in compute_kl_and_mean_shift
        """
        plt.figure(figsize=(10, 4.5))
        x = df["run"].values

        plt.plot(x, df["KL"].values, marker="o", linewidth=1.5, label="KL divergence")
        plt.xlabel("Run index")
        plt.ylabel("KL divergence")
        plt.title(f"Belief-Update Stability Across Repeated Runs: {drug} → {disease}")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(outpath, dpi=300)
        plt.close()

    @staticmethod
    def _stability_summary(df: pd.DataFrame) -> Dict[str, str]:
        """
        Compact stats you can paste into the manuscript or supplement.
        """
        def fmt(col: str) -> str:
            s = df[col].dropna()
            if s.empty:
                return "NA"
            mu = float(s.mean())
            sd = float(s.std(ddof=1)) if len(s) > 1 else 0.0
            mn = float(s.min())
            mx = float(s.max())
            return f"{mu:.6f} ± {sd:.6f} | {mn:.6f}..{mx:.6f}"

        return {
            "p_raw": fmt("p_raw"),
            "p_penalised": fmt("p_penalised"),
            "p_final": fmt("p_final"),
            "gamma": fmt("gamma"),
            "post_mean": fmt("post_mean"),
            "KL": fmt("KL"),
            "M": fmt("M"),
            "T": fmt("T"),
            "A": fmt("A"),
            "N": fmt("N"),
        }

    # -----------------------------------------------------------------
    # Existing drug evaluation (unchanged)
    # -----------------------------------------------------------------
    def evaluate_drug(self, drug: str, diseases: List[str], top_k: int = 5) -> None:
        print(f"\n=== Evaluating repurposing for: {drug} ===")
        results: List[Dict[str, Any]] = []

        for disease in diseases:
            drug_l, disease_l = drug.lower(), disease.lower()
            if (drug_l, disease_l) in self.existing_pairs:
                print(f"[INFO] Already known: {drug} → {disease}")
                continue

            comp = self.compute_components(drug, disease)
            self._print_diagnostics(drug, disease, comp)

            results.append(
                {
                    "disease": disease,
                    "post_mean": comp["post_mean"],
                    "components": comp,
                }
            )

            self._save_run_log(drug, disease, comp)

        if not results:
            print("No new repurposing candidates.")
            return

        results.sort(key=lambda r: r["post_mean"], reverse=True)
        top = results[: int(max(1, top_k))]

        print("\n" + "=" * 70)
        print("TOP CANDIDATES (by posterior mean)")
        print("=" * 70)
        for r in top:
            comp = r["components"]
            print(
                f"{drug} → {r['disease']}: "
                f"PosteriorMean={comp['post_mean']:.6f} | "
                f"P_final={comp['p_final']:.6f} | "
                f"M={comp['M']} | "
                f"gamma={comp['gamma'] if comp['gamma'] is not None else 'None'} | "
                f"p_like={comp['p_likelihood']:.6f}"
            )

        # Plot + KL diagnostics for top-k
        kl_rows = []
        for r in top:
            dis = r["disease"]
            comp = r["components"]
            print(f"\n[INFO] Plotting: {drug} → {dis}")
            x, prior_pdf, post_pdf = self.plot_distributions(comp, drug, dis)

            like_pdf = beta.pdf(x, a=comp["like_a"], b=comp["like_b"])
            like_pdf = like_pdf / trapezoid(like_pdf, x)
            e_like = trapezoid(x * like_pdf, x)

            kl, mu_prior, mu_post, delta = self.compute_kl_and_mean_shift(prior_pdf, post_pdf, x)

            kl_rows.append(
                {
                    "Disease": dis,
                    "KL Divergence": kl,
                    "E[Prior]": mu_prior,
                    "E[Likelihood]": float(round(e_like, 6)),
                    "E[Posterior]": mu_post,
                    "Δμ": delta,
                    "P_raw": float(round(comp["p_raw"], 6)),
                    "P_pen": float(round(comp["p_penalised"], 6)),
                    "P_final": float(round(comp["p_final"], 6)),
                    "gamma": comp["gamma"],
                    "M": comp["M"],
                    "c(M)": float(round(comp["cM"], 6)),
                    "p_like": float(round(comp["p_likelihood"], 6)),
                }
            )

        print("\n" + "=" * 70)
        print("KL DIVERGENCE AND MEAN SHIFT (TOP-K)")
        print("=" * 70)
        print(pd.DataFrame(kl_rows).to_string(index=False))


# ---------------------------------------------------------------------
# Script usage
# ---------------------------------------------------------------------
if __name__ == "__main__":
    weights_dict = {
        "GraphDistanceToIndication": 0.1148,
        "RandomWalkScore": 0.247,
        "StructuralLikelihood": -0.1154,
        "PreferentialAttachment": -0.041,
        "KatzSimilarity": 1.6515,
    }

    cfg = PredictorConfig(
        pubmed_max_articles=150,
        pubmed_filter_level="high",
        pubmed_years_back=10,
        pubmed_use_cache=True,
        llm_model="gpt-4o-mini",
        llm_batch_size=5,
        llm_delay_s=2.0,
        llm_max_retries=2,
        cmax=200.0,
        tau=25.0,
        likelihood_strength=50.0,
        likelihood_intercept=0.0,
        save_run_logs=True,
        stability_dir="stability",
        stability_n_runs=10,
    )

    predictor = BayesianRepurposingPredictor(
        known_path="graph/graph_features_known.csv",
        unknown_path="graph/graph_features_unknown.csv",
        weights_dict=weights_dict,
        cfg=cfg,
    )

    drug = "dexamethasone"
    disease = "COVID-19"


    # 10-run stability evaluation (saves CSV + 3 plots)
    predictor.evaluate_pair_stability(drug, disease, n_runs=10)
