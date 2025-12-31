"""
bayesian_predictor.py

Bayesian drug repurposing predictor that integrates:
  1) Semantic priors from PubMed abstracts (LLMClassifier in pubmed_utils.py)
  2) Graph-feature evidence (weighted score mapped to probability)
  3) Beta pseudo-count fusion (prior + evidence) for posterior inference
  4) Transparent diagnostics: prints and logs priors, counts, gamma, evidence strength, α/β, etc.

Key fixes vs. your previous version:
- Aligns with updated pubmed_utils.py API (search_cfg/llm_cfg; max_articles; filter_level; use_cache).
- Uses evidence-scaled prior concentration c(M) based on total_articles (M), not a hard-coded 100.
- Uses a bounded mapping (sigmoid) to convert graph feature score to probability before Beta construction.
- Uses one consistent likelihood_strength in compute + plotting (no 50/100 mismatch).
- Propagates and prints gamma (side-effect overlap confidence) and M (articles used).
- Consistent Beta parameterization: α = 1 + strength * p ; β = 1 + strength * (1-p).
- Missing features -> neutral likelihood probability (0.5) rather than 1.0.
- Saves per-run detailed diagnostics JSON (optional) for reproducibility.

Assumptions:
- pubmed_utils.LLMClassifier.build_semantic_prior returns:
    {
      "prior": float,
      "penalised_prior": float,
      "enhanced_prior": float,   # treated as P_final
      "gamma": float|None,
      "raw_counts": Counter,
      "labelled_abstracts": list,
      "total_articles": int
    }
- side_effect_updater.update_prior returns dict with keys including "p_final" and "gamma"
  (already handled inside pubmed_utils.py).
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

    - cmax controls maximum pseudo-count weight you allow PubMed priors to have.
    - tau controls how quickly you approach cmax.

    You should set (cmax, tau) to match your manuscript values.
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
    likelihood_intercept: float = 0.0  # can be tuned/calibrated

    # Output & reproducibility
    plots_dir: str = "plots"
    logs_dir: str = "runs"
    save_run_logs: bool = True


# ---------------------------------------------------------------------
# PubMed semantic prior wrapper (aligned with new pubmed_utils.py)
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

    @staticmethod
    def load_combined_features(known_csv: str, unknown_csv: str) -> pd.DataFrame:
        df_known = pd.read_csv(known_csv)
        df_unknown = pd.read_csv(unknown_csv)
        df = pd.concat([df_known, df_unknown], ignore_index=True)
        df["Drug"] = df["Drug"].astype(str).str.strip().str.lower()
        df["Disease"] = df["Disease"].astype(str).str.strip().str.lower()
        return df

    def _pubmed_prior(self, drug: str, disease: str) -> Dict[str, Any]:
        """
        Retrieve semantic prior with fully specified parameters.
        """
        return self.classifier.build_semantic_prior(
            drug=drug,
            disease=disease,
            max_articles=int(self.cfg.pubmed_max_articles),
            filter_level=str(self.cfg.pubmed_filter_level),
            use_cache=bool(self.cfg.pubmed_use_cache),
        )

    def _graph_score(self, drug: str, disease: str) -> Tuple[float, Dict[str, float]]:
        """
        Compute linear feature score s = intercept + Σ(w_i * x_i).

        Returns:
          (score, feature_values_used)
        """
        drug_l, disease_l = drug.lower(), disease.lower()
        row = self.feature_df[(self.feature_df["Drug"] == drug_l) & (self.feature_df["Disease"] == disease_l)]
        if row.empty:
            # Neutral baseline if features are missing
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
        """
        Map unbounded score to probability via sigmoid.
        """
        return clamp01(sigmoid(score))

    def compute_components(self, drug: str, disease: str) -> Dict[str, Any]:
        """
        Compute prior, likelihood, posterior parameters and diagnostics.

        Returns a dict with everything needed for reporting and plotting.
        """
        # ---- Semantic prior from PubMed
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
        p_final = clamp01(prior_result.get("enhanced_prior", 0.5))  # manuscript P_final
        gamma = prior_result.get("gamma", None)
        gamma = clamp01(gamma) if gamma is not None else None

        counts = prior_result.get("raw_counts", {})
        try:
            counts = dict(counts)
        except Exception:
            counts = {}

        M = int(prior_result.get("total_articles", 0))

        # Evidence-scaled concentration for the prior
        cM = concentration_c(M, cmax=self.cfg.cmax, tau=self.cfg.tau)

        prior_a, prior_b = beta_params_from_prob(p_final, strength=cM)

        # ---- Graph evidence likelihood as probability
        score, feature_vals = self._graph_score(drug, disease)
        p_like = self._graph_probability(score)

        like_a, like_b = beta_params_from_prob(p_like, strength=self.cfg.likelihood_strength)

        # ---- Pseudo-count fusion (subtract the +1 baselines so we only add evidence)
        post_a = prior_a + (like_a - 1.0)
        post_b = prior_b + (like_b - 1.0)

        post_mean = post_a / (post_a + post_b)

        return {
            # priors
            "p_raw": p_raw,
            "p_penalised": p_pen,
            "p_final": p_final,
            "gamma": gamma,
            "counts": counts,
            "M": M,
            "cM": cM,
            "prior_a": prior_a,
            "prior_b": prior_b,
            # likelihood
            "graph_score": score,
            "p_likelihood": p_like,
            "likelihood_strength": float(self.cfg.likelihood_strength),
            "like_a": like_a,
            "like_b": like_b,
            "feature_values": feature_vals,
            # posterior
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
        """
        Plot normalized prior, likelihood, posterior Beta PDFs, with consistent parameters.
        """
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
        plt.grid(alpha=0.3)
        plt.tight_layout()

        fname = f"{drug.strip().lower().replace(' ', '_')}_{disease.strip().lower().replace(' ', '_')}.png"
        fpath = os.path.join(self.cfg.plots_dir, fname)
        plt.savefig(fpath)
        plt.close()

        return x, prior_pdf, post_pdf

    def _print_diagnostics(self, drug: str, disease: str, comp: Dict[str, Any]) -> None:
        """
        Print all key calculation components for transparency.
        """
        print("\n" + "-" * 70)
        print(f"PAIR: {drug} → {disease}")
        print("-" * 70)

        # PubMed prior diagnostics
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

        # Graph likelihood diagnostics
        print("\n[Graph evidence / Likelihood]")
        print(f"  Intercept:          {self.cfg.likelihood_intercept:.6f}")
        print(f"  Linear score:       {comp['graph_score']:.6f}")
        print(f"  p_likelihood(sigmoid): {comp['p_likelihood']:.6f}")
        print(f"  Likelihood strength:   {comp['likelihood_strength']:.6f}")
        print(f"  Likelihood Beta:    alpha={comp['like_a']:.6f}, beta={comp['like_b']:.6f}")
        if comp["feature_values"]:
            print("  Feature values used:")
            for k, v in comp["feature_values"].items():
                w = self.weights_dict.get(k, 0.0)
                print(f"    - {k}: value={v:.6f}, weight={float(w):.6f}, contrib={(float(w)*float(v)):.6f}")
        else:
            print("  Feature values used: (missing) -> neutral baseline")

        # Posterior diagnostics
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
        # IMPORTANT: set these two to match your manuscript
        cmax=200.0,
        tau=25.0,
        likelihood_strength=50.0,
        likelihood_intercept=0.0,
        save_run_logs=True,
    )

    predictor = BayesianRepurposingPredictor(
        known_path="graph/graph_features_known.csv",
        unknown_path="graph/graph_features_unknown.csv",
        weights_dict=weights_dict,
        cfg=cfg,
    )

    drugs = ["Azithromycin"]
    diseases = [
        "ureterolithiasis",
        "squamous cell carcinoma of head and neck",
        "frailty",
        "venous thrombosis",
        "constriction, pathologic", 
        "pregnancy in diabetics",
        "gout",
        "body weight",
        "adenomatous polyps",
        "albuminuria"     
    ]

    for drug in drugs:
        predictor.evaluate_drug(drug, diseases, top_k=10)
