import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid
from pubmed_utils import pubmed_count_requests


# ========== Data Loaders ========== #
def load_trial_links(data_dir="data"):
    trial_links = defaultdict(set)
    all_diseases = set()
    all_drugs = set()
    for file in os.listdir(data_dir):
        if file.endswith(".json"):
            disease = file.replace(".json", "").replace("_", " ").lower()
            all_diseases.add(disease)
            with open(os.path.join(data_dir, file), "r", encoding="utf-8") as f:
                trials = json.load(f)
                for trial in trials:
                    for drug in trial.get("interventions", []):
                        drug = drug.lower().strip()
                        trial_links[(drug, disease)].add(trial["nctId"])
                        all_drugs.add(drug)
    return trial_links, sorted(all_diseases), sorted(all_drugs)


def load_centrality_scores(csv_path="graph/drug_centrality.csv"):
    df = pd.read_csv(csv_path)
    return df.set_index("Drug").to_dict("index")


# ========== Bayesian Components ========== #
def compute_components(drug, disease, centrality_scores):
    try:
        pubmed_joint = pubmed_count_requests(drug, disease)
        pubmed_d = pubmed_count_requests(drug, "")
        pubmed_c = pubmed_count_requests(disease, "")
        prior = pubmed_joint / (pubmed_d + pubmed_c + 1e-5)
    except Exception:
        prior = 0.0

    scores = centrality_scores.get(drug.lower(), {})
    degree = scores.get("DegreeCentrality", 0)
    eigen = scores.get("EigenvectorCentrality", 0)
    between = scores.get("BetweennessCentrality", 0)

    likelihood = 1 + 0.2 * degree + 0.3 * eigen + 0.1 * between
    posterior_strength = prior * likelihood

    return prior, likelihood, posterior_strength


def compute_kl_and_mean_shift(prior_dist, posterior_dist, x):
    prior_dist = np.maximum(prior_dist, 1e-10)
    posterior_dist = np.maximum(posterior_dist, 1e-10)
    prior_dist /= trapezoid(prior_dist, x)
    posterior_dist /= trapezoid(posterior_dist, x)
    kl = trapezoid(rel_entr(posterior_dist, prior_dist), x)
    mu_prior = trapezoid(x * prior_dist, x)
    mu_post = trapezoid(x * posterior_dist, x)
    delta_mu = mu_post - mu_prior
    return round(kl, 4), round(mu_prior, 4), round(mu_post, 4), round(delta_mu, 4)


# ========== Plotting ========== #
def plot_distributions(prior, likelihood, posterior, drug, disease):
    x = np.linspace(0, 1, 500)

    prior_a = prior * 100 + 1
    prior_b = (1 - prior) * 100 + 1
    prior_dist = beta.pdf(x, a=prior_a, b=prior_b)

    likelihood_center = min(max(likelihood / 10, 0.01), 0.99)
    likelihood_a = likelihood_center * 100
    likelihood_b = (1 - likelihood_center) * 100
    likelihood_dist = beta.pdf(x, a=likelihood_a + 1, b=likelihood_b + 1)

    posterior_a = prior_a * likelihood
    posterior_b = prior_b
    posterior_dist = beta.pdf(x, a=posterior_a, b=posterior_b)

    prior_plot = prior_dist / prior_dist.max()
    likelihood_plot = likelihood_dist / likelihood_dist.max()
    posterior_plot = posterior_dist / posterior_dist.max()

    plt.figure(figsize=(8, 5))
    plt.plot(x, prior_plot, label="Prior (Literature)", color="blue", linewidth=2)
    plt.plot(x, likelihood_plot, label="Likelihood (Network)", color="red", linewidth=2)
    plt.plot(x, posterior_plot, label="Posterior", color="purple", linewidth=2)
    plt.title(f"Bayesian Inference for {drug.title()} → {disease.title()}", fontsize=13)
    plt.xlabel("θ (latent association strength)")
    plt.ylabel("Relative Density")
    plt.legend(loc="upper right", frameon=False)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()

    return x, prior_dist, posterior_dist


# ========== Main: Repurposing Prediction ========== #
if __name__ == "__main__":
    print("Loading data...")
    trial_links, diseases, drugs = load_trial_links("data")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    drug = "dexamethasone"
    results = []

    print(f"\nPredicting repurposing candidates for: {drug.title()}")
    for disease in diseases:
        # Focus on diseases not yet trialed for this drug
        if (drug, disease) not in trial_links:
            prior, likelihood, posterior = compute_components(drug, disease, centrality_scores)
            results.append({
                "disease": disease,
                "prior": prior,
                "likelihood": likelihood,
                "posterior": posterior
            })

    if not results:
        print("No new repurposing candidates found.")
        exit()

    # Normalize posterior scores
    total = sum(r["posterior"] for r in results)
    for r in results:
        r["probability"] = round(r["posterior"] / total, 6) if total > 0 else 0.0

    results = sorted(results, key=lambda r: r["probability"], reverse=True)
    top_diseases = results[:4]

    print("\nTop Repurposing Candidates (Normalized Posterior):")
    for r in top_diseases:
        print(f"{drug.title()} → {r['disease'].title()} — Probability: {r['probability']}")

    kl_outputs = []
    for row in top_diseases:
        disease = row["disease"]
        prior = row["prior"]
        likelihood = row["likelihood"]
        posterior = row["posterior"]

        print(f"\nPlotting and computing KL for: {drug.title()} → {disease.title()}")
        x, prior_dist, posterior_dist = plot_distributions(prior, likelihood, posterior, drug, disease)
        kl, mu_prior, mu_post, delta_mu = compute_kl_and_mean_shift(prior_dist, posterior_dist, x)
        kl_outputs.append({
            "Disease": disease.title(),
            "KL Divergence": kl,
            "E[Prior]": mu_prior,
            "E[Posterior]": mu_post,
            "Δμ": delta_mu
        })

    print("\nKL Divergence and Mean Shift (Top Candidates):")
    df = pd.DataFrame(kl_outputs)
    print(df.to_string(index=False))
