import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, norm
from scipy.special import rel_entr
from scipy.integrate import trapezoid
from pubmed_utils import pubmed_count_requests
from bayesian_predictor import compute_posterior, load_centrality_scores, load_trial_links


def compute_components(drug, disease, centrality_scores):
    """Compute prior, likelihood, and posterior for visualisation."""
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
    likelihood = 1 + 2 * degree + 3 * eigen + 2 * between

    posterior = prior * likelihood
    return prior, likelihood, posterior


def compute_kl_and_mean_shift(prior_dist, posterior_dist, x):
    """Compute KL divergence and mean shift between prior and posterior."""
    prior_dist = np.maximum(prior_dist, 1e-10)
    posterior_dist = np.maximum(posterior_dist, 1e-10)

    prior_dist /= trapezoid(prior_dist, x)
    posterior_dist /= trapezoid(posterior_dist, x)

    kl = trapezoid(rel_entr(posterior_dist, prior_dist), x)
    mu_prior = trapezoid(x * prior_dist, x)
    mu_post = trapezoid(x * posterior_dist, x)
    delta_mu = mu_post - mu_prior

    return round(kl, 4), round(mu_prior, 4), round(mu_post, 4), round(delta_mu, 4)


def plot_distributions(prior, likelihood, posterior, drug, disease):
    """Plot Prior, Likelihood, and Posterior curves, return distributions for analysis."""
    x = np.linspace(0, 1, 500)

    prior_dist = beta.pdf(x, a=prior * 50 + 1, b=50)
    likelihood_dist = norm.pdf(x, loc=min(likelihood / 20, 1), scale=0.1)
    posterior_dist = norm.pdf(x, loc=min(posterior / 20, 1), scale=0.1)

    # Normalise for plotting
    prior_plot = prior_dist / prior_dist.max()
    likelihood_plot = likelihood_dist / likelihood_dist.max()
    posterior_plot = posterior_dist / posterior_dist.max()

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(x, prior_plot, label="Prior (Literature)", color="blue", linewidth=2)
    plt.plot(x, likelihood_plot, label="Likelihood (Network)", color="red", linewidth=2)
    plt.plot(x, posterior_plot, label="Posterior", color="purple", linewidth=2)

    plt.title(f"Bayesian Inference Components for {drug.title()} and {disease.title()}", fontsize=13)
    plt.xlabel("θ (latent association strength)")
    plt.ylabel("Relative Density")
    plt.legend(loc="upper right", frameon=False)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()

    return x, prior_dist, posterior_dist


if __name__ == "__main__":
    trial_links, diseases, drugs = load_trial_links("data")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    drug = "dexamethasone"
    disease = "obesity"

    if (drug, disease) in trial_links:
        print(f"{drug.title()} has already been tested for {disease.title()} — please choose a different pair.")
    else:
        print(f"\nAnalysing: {drug.title()} → {disease.title()}")
        prior, likelihood, posterior = compute_components(drug, disease, centrality_scores)
        x, prior_dist, posterior_dist = plot_distributions(prior, likelihood, posterior, drug, disease)

        kl, mu_prior, mu_post, delta_mu = compute_kl_and_mean_shift(prior_dist, posterior_dist, x)
        print(f"\n KL Divergence (Posterior || Prior): {kl}")
        print(f"Expected Value of Prior: {mu_prior}")
        print(f"Expected Value of Posterior: {mu_post}")
        print(f"Δμ (Shift in Mean): {delta_mu}")
