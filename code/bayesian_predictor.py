import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid  # ✅ Use this instead of deprecated np.trapz
from pubmed_utils import build_semantic_prior

# ========= Load Prior Condition–Drug Pairs ========= #
def load_existing_trial_pairs(json_path="processed_data/condition_drug_pairs.json"):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return set((entry["intervention"], entry["condition"]) for entry in data)

def load_centrality_scores(csv_path="graph/drug_centrality.csv"):
    df = pd.read_csv(csv_path)

    # Select only numeric centrality columns to normalize
    centrality_cols = [
        "DegreeCentrality", "EigenvectorCentrality", "BetweennessCentrality",
        "ClusteringCoefficient", "RandomWalkCentrality"
    ]

    # Normalize each column to [0, 1] (Min-Max Scaling)
    for col in centrality_cols:
        min_val = df[col].min()
        max_val = df[col].max()
        if max_val > min_val:
            df[col] = (df[col] - min_val) / (max_val - min_val)
        else:
            df[col] = 0.0  # If constant, set to 0

    return df.set_index("Drug").to_dict("index")

# ========== Bayesian Components ========== #
def compute_components(drug, disease, centrality_scores):
    try:
        prior_result = build_semantic_prior(drug, disease, max_count=5)
        raw_prior = prior_result["prior"] if prior_result else 0.0
        penalised_prior = prior_result["penalised_prior"] if prior_result else 0.0
    except Exception as e:
        print(f"[!] Error computing semantic prior for {drug} → {disease}:", e)
        raw_prior = 0.0
        penalised_prior = 0.0

    scores = centrality_scores.get(drug, {})
    degree = scores.get("DegreeCentrality", 0)
    eigen = scores.get("EigenvectorCentrality", 0)
    between = scores.get("BetweennessCentrality", 0)
    clustering = scores.get("ClusteringCoefficient", 0)
    random_walk = scores.get("RandomWalkCentrality", 0)
    
    # Randon forest Feature Importances:
    # DegreeCentrality: 0.1925
    # EigenvectorCentrality: 0.2283
    # BetweennessCentrality: 0.1281
    # ClusteringCoefficient: 0.1573
    # RandomWalkCentrality: 0.2937

    likelihood = 1 + 0.1925 * degree + 0.2283 * eigen + 0.1281 * between + 0.1573 * clustering + 0.2937 * random_walk
    
    # Convert likelihood to pseudo-counts
    likelihood_strength = 10  # Emperically chosen
    likelihood_a = likelihood * likelihood_strength
    likelihood_b = (1 - (likelihood - 1)) * likelihood_strength

    # Use penalised prior for calculations
    prior_a = penalised_prior * 100 + 1
    prior_b = (1 - penalised_prior) * 100 + 1

    # Add pseudo-counts from likelihood
    posterior_a = prior_a + likelihood_a
    posterior_b = prior_b + likelihood_b

    posterior_mean = posterior_a / (posterior_a + posterior_b)

    return raw_prior, penalised_prior, likelihood, posterior_mean, prior_a, prior_b, posterior_a, posterior_b

def compute_kl_and_mean_shift(prior_dist, posterior_dist, x):
    prior_dist = np.clip(prior_dist, 1e-10, None)
    posterior_dist = np.clip(posterior_dist, 1e-10, None)

    prior_dist /= trapezoid(prior_dist, x)
    posterior_dist /= trapezoid(posterior_dist, x)

    kl = trapezoid(rel_entr(posterior_dist, prior_dist), x)
    mu_prior = trapezoid(x * prior_dist, x)
    mu_post = trapezoid(x * posterior_dist, x)
    delta_mu = mu_post - mu_prior

    return round(kl, 4), round(mu_prior, 4), round(mu_post, 4), round(delta_mu, 4)

# ========== Plotting ========== #
def plot_distributions(prior_a, prior_b, likelihood, posterior_a, posterior_b, drug, disease):
    x = np.linspace(0.001, 0.999, 1000)  # Better resolution and avoids 0s

    prior_dist = beta.pdf(x, a=prior_a, b=prior_b)

    likelihood_center = min(max(likelihood / 10, 0.01), 0.99)
    likelihood_a = likelihood_center * 100
    likelihood_b = (1 - likelihood_center) * 100
    likelihood_dist = beta.pdf(x, a=likelihood_a + 1, b=likelihood_b + 1)

    posterior_dist = beta.pdf(x, a=posterior_a, b=posterior_b)

    prior_plot = prior_dist / prior_dist.max()
    likelihood_plot = likelihood_dist / likelihood_dist.max()
    posterior_plot = posterior_dist / posterior_dist.max()

    plt.figure(figsize=(8, 5))
    plt.plot(x, prior_plot, label="Prior (Penalised)", color="blue", linewidth=2)
    plt.plot(x, likelihood_plot, label="Likelihood (Network)", color="red", linewidth=2)
    plt.plot(x, posterior_plot, label="Posterior", color="purple", linewidth=2)
    plt.title(f"Bayesian Inference for {drug} → {disease}", fontsize=13)
    plt.xlabel("θ (Latent Association Strength)")
    plt.ylabel("Relative Density")
    plt.legend(loc="upper right", frameon=False)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()

    return x, prior_dist, posterior_dist

# ========== Main ========== #
if __name__ == "__main__":
    print("Loading prior trial pairs and centrality scores...")
    existing_pairs = load_existing_trial_pairs("processed_data/condition_drug_pairs.json")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    # Custom user input
    drugs_to_evaluate = ["Thalidomide"]
    diseases_to_evaluate = ["Multiple Myeloma", "Arthritis, Rheumatoid", "COVID-19","Glioblastoma","Lupus Erythematosus, Systemic"]

    for drug in drugs_to_evaluate:
        print(f"\n=== Predicting repurposing candidates for: {drug} ===")
        results = []
        for disease in diseases_to_evaluate:
            if (drug, disease) in existing_pairs:
                print(f"Trial has already been done on {drug} → {disease}")
                continue

            raw_prior, penalised_prior, likelihood, posterior_mean, prior_a, prior_b, post_a, post_b = compute_components(
                drug, disease, centrality_scores
            )

            results.append({
                "disease": disease,
                "raw_prior": raw_prior,
                "penalised_prior": penalised_prior,
                "likelihood": likelihood,
                "posterior_mean": posterior_mean,
                "prior_a": prior_a,
                "prior_b": prior_b,
                "posterior_a": post_a,
                "posterior_b": post_b
            })

        if not results:
            print("No new repurposing candidates found.")
            continue

        results = sorted(results, key=lambda r: r["posterior_mean"], reverse=True)
        top_diseases = results[:5]

        print("\nTop Repurposing Candidates (Posterior Mean):")
        for r in top_diseases:
            print(f"{drug} → {r['disease']} — Posterior Mean: {round(r['posterior_mean'], 6)} | Raw Prior: {round(r['raw_prior'], 4)}")

        kl_outputs = []
        for row in top_diseases:
            disease = row["disease"]
            print(f"\nPlotting and computing KL for: {drug} → {disease}")
            x, prior_dist, posterior_dist = plot_distributions(
                row["prior_a"], row["prior_b"], row["likelihood"],
                row["posterior_a"], row["posterior_b"],
                drug, disease
            )
            kl, mu_prior, mu_post, delta_mu = compute_kl_and_mean_shift(prior_dist, posterior_dist, x)
            kl_outputs.append({
                "Disease": disease,
                "KL Divergence": kl,
                "E[Prior]": mu_prior,
                "E[Posterior]": mu_post,
                "Δμ": delta_mu,
                "Raw Prior": round(row["raw_prior"], 3),
                "Penalised Prior": round(row["penalised_prior"], 3)
            })

        print("\nKL Divergence and Mean Shift (Top Candidates):")
        df = pd.DataFrame(kl_outputs)
        print(df.to_string(index=False))

