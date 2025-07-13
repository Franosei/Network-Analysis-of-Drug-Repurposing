import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid
from pubmed_utils import LLMClassifier

classifier = LLMClassifier(delay=4)

def build_semantic_prior(drug, disease, max_count=30):
    """
    Builds a semantic prior based on LLM classification of PubMed abstracts.

    Parameters:
        drug (str): Drug name.
        disease (str): Disease name.
        max_count (int): Maximum number of PubMed abstracts to use.

    Returns:
        dict: Contains 'prior', 'penalised_prior', and 'enhanced_prior' fields.
    """
    return classifier.build_semantic_prior(drug, disease, max_count=max_count)


class BayesianRepurposingPredictor:
    """
    A Bayesian engine for drug repurposing that integrates:
      - Semantic priors from LLM-based PubMed analysis
      - Graph-based topological features from known and unknown associations
      - Weighted likelihoods using learned coefficients
      - Posterior inference using Beta distributions
      - KL divergence and mean shift for visual and statistical insight

    Attributes:
        feature_df (pd.DataFrame): Combined graph features for drug-disease pairs.
        weights_dict (dict): Feature weightings learned from training data.
        existing_pairs (set): Set of known drug-disease tuples with label = 1.
    """

    def __init__(self, known_path, unknown_path, weights_dict):
        """
        Initializes the predictor with graph features and feature weights.

        Parameters:
            known_path (str): Path to known pairs CSV file.
            unknown_path (str): Path to unknown pairs CSV file.
            weights_dict (dict): Dictionary of feature weights.
        """
        self.feature_df = self.load_combined_features(known_path, unknown_path)
        self.weights_dict = weights_dict
        self.existing_pairs = set(
            zip(
                self.feature_df[self.feature_df["Label"] == 1]["Drug"],
                self.feature_df[self.feature_df["Label"] == 1]["Disease"]
            )
        )

    def load_combined_features(self, known_csv, unknown_csv):
        """
        Loads and combines known and unknown graph feature files.

        Parameters:
            known_csv (str): Path to CSV with known (label=1) pairs.
            unknown_csv (str): Path to CSV with unknown (label=0) pairs.

        Returns:
            pd.DataFrame: Merged and cleaned feature DataFrame.
        """
        df_known = pd.read_csv(known_csv)
        df_unknown = pd.read_csv(unknown_csv)
        df = pd.concat([df_known, df_unknown], ignore_index=True)
        df["Drug"] = df["Drug"].str.strip().str.lower()
        df["Disease"] = df["Disease"].str.strip().str.lower()
        return df

    def get_likelihood_from_features(self, drug, disease):
        """
        Computes the likelihood score based on weighted graph features.

        Parameters:
            drug (str): Drug name.
            disease (str): Disease name.

        Returns:
            float: Likelihood value.
        """
        drug, disease = drug.lower(), disease.lower()
        row = self.feature_df[
            (self.feature_df["Drug"] == drug) & 
            (self.feature_df["Disease"] == disease)
        ]
        if row.empty:
            print(f"[!] Missing features for {drug} → {disease}")
            return 1.0

        likelihood = 1.0
        for feature, weight in self.weights_dict.items():
            value = row.iloc[0].get(feature, 0)
            likelihood += weight * value
        return likelihood

    def compute_components(self, drug, disease):
        """
        Calculates prior, likelihood, and posterior Beta distribution parameters.

        Parameters:
            drug (str): Drug name.
            disease (str): Disease name.

        Returns:
            tuple: Contains raw prior, penalised prior, enhanced prior, likelihood,
                   posterior mean, prior_a, prior_b, posterior_a, posterior_b
        """
        try:
            prior_result = build_semantic_prior(drug, disease)
            raw_prior = prior_result["prior"]
            penalised_prior = prior_result["penalised_prior"]
            enhanced_prior = prior_result["enhanced_prior"]
        except Exception as e:
            print(f"[!] Error computing prior for {drug} → {disease}: {e}")
            raw_prior = penalised_prior = enhanced_prior = 0.0

        likelihood = self.get_likelihood_from_features(drug, disease)

        likelihood_strength = 50
        likelihood_center = min(max(likelihood, 0.01), 0.99)
        likelihood_a = likelihood_center * likelihood_strength
        likelihood_b = (1 - likelihood_center) * likelihood_strength

        prior_a = enhanced_prior * 100 + 1
        prior_b = (1 - enhanced_prior) * 100 + 1

        posterior_a = prior_a + likelihood_a
        posterior_b = prior_b + likelihood_b

        posterior_mean = posterior_a / (posterior_a + posterior_b)

        return raw_prior, penalised_prior, enhanced_prior, likelihood, posterior_mean, prior_a, prior_b, posterior_a, posterior_b

    def compute_kl_and_mean_shift(self, prior_dist, posterior_dist, x):
        """
        Computes KL divergence and mean shift between two distributions.

        Parameters:
            prior_dist (np.ndarray): Beta PDF of prior.
            posterior_dist (np.ndarray): Beta PDF of posterior.
            x (np.ndarray): Support grid.

        Returns:
            tuple: KL divergence, mean prior, mean posterior, delta mean.
        """
        prior_dist = np.clip(prior_dist, 1e-10, None)
        posterior_dist = np.clip(posterior_dist, 1e-10, None)

        prior_dist /= trapezoid(prior_dist, x)
        posterior_dist /= trapezoid(posterior_dist, x)

        kl = trapezoid(rel_entr(posterior_dist, prior_dist), x)
        mu_prior = trapezoid(x * prior_dist, x)
        mu_post = trapezoid(x * posterior_dist, x)
        delta_mu = mu_post - mu_prior

        return round(kl, 4), round(mu_prior, 4), round(mu_post, 4), round(delta_mu, 4)

    def plot_distributions(self, prior_a, prior_b, likelihood, posterior_a, posterior_b, drug, disease):
        """
        Plots the normalized prior, likelihood, and posterior Beta distributions
        with 95% confidence intervals and saves the figure to disk.

        Parameters:
            prior_a (float): Alpha of prior Beta.
            prior_b (float): Beta of prior Beta.
            likelihood (float): Likelihood score.
            posterior_a (float): Alpha of posterior Beta.
            posterior_b (float): Beta of posterior Beta.
            drug (str): Drug name.
            disease (str): Disease name.

        Returns:
            tuple: x-axis array, prior PDF, posterior PDF
        """
        x = np.linspace(0.001, 0.999, 1000)
        prior_dist = beta.pdf(x, a=prior_a, b=prior_b)
        likelihood_center = min(max(likelihood, 0.01), 0.99)
        la = likelihood_center * 100
        lb = (1 - likelihood_center) * 100
        likelihood_dist = beta.pdf(x, a=la + 1, b=lb + 1)
        posterior_dist = beta.pdf(x, a=posterior_a, b=posterior_b)

        # 95% credible intervals
        prior_ci = beta.interval(0.95, prior_a, prior_b)
        likelihood_ci = beta.interval(0.95, la + 1, lb + 1)
        posterior_ci = beta.interval(0.95, posterior_a, posterior_b)

        plt.figure(figsize=(10, 6))
        plt.plot(x, prior_dist / prior_dist.max(), label="Prior", color="blue")
        plt.plot(x, likelihood_dist / likelihood_dist.max(), label="Likelihood", color="red")
        plt.plot(x, posterior_dist / posterior_dist.max(), label="Posterior", color="purple")

        # Shade 95% CIs
        plt.axvspan(*prior_ci, color="blue", alpha=0.1)
        plt.axvspan(*likelihood_ci, color="red", alpha=0.1)
        plt.axvspan(*posterior_ci, color="purple", alpha=0.1)

        plt.title(f"Bayesian Update: {drug} → {disease}")
        plt.xlabel("Association Strength")
        plt.ylabel("Normalized Density")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()

        os.makedirs("plots", exist_ok=True)
        filename = f"{drug.strip().lower().replace(' ', '_')}_{disease.strip().lower().replace(' ', '_')}.png"
        filepath = os.path.join("plots", filename)
        plt.savefig(filepath)
        plt.close()

        return x, prior_dist, posterior_dist

    def evaluate_drug(self, drug, diseases):
        """
        Evaluates candidate diseases for a given drug using Bayesian inference.

        Parameters:
            drug (str): Drug name.
            diseases (list[str]): List of diseases to evaluate.

        Prints:
            Top 5 candidates, posterior means, and KL divergence stats.
        """
        print(f"\n=== Evaluating repurposing for: {drug} ===")
        results = []

        for disease in diseases:
            drug_l, disease_l = drug.lower(), disease.lower()
            if (drug_l, disease_l) in self.existing_pairs:
                print(f"Already known: {drug} → {disease}")
                continue

            raw_prior, penalised_prior, enhanced_prior, likelihood, post_mean, pa, pb, post_a, post_b = self.compute_components(drug, disease)
            results.append({
                "disease": disease,
                "raw_prior": raw_prior,
                "penalised_prior": penalised_prior,
                "enhanced_prior": enhanced_prior,
                "likelihood": likelihood,
                "posterior_mean": post_mean,
                "prior_a": pa,
                "prior_b": pb,
                "posterior_a": post_a,
                "posterior_b": post_b
            })

        if not results:
            print("No new repurposing candidates.")
            return

        results = sorted(results, key=lambda r: r["posterior_mean"], reverse=True)
        top = results[:5]

        print("\nTop Candidates:")
        for r in top:
            print(f"{drug} → {r['disease']}: Posterior Mean={round(r['posterior_mean'], 4)}, Raw Prior={round(r['raw_prior'], 4)}")

        kl_metrics = []
        for r in top:
            print(f"\nPlotting: {drug} → {r['disease']}")
            x, prior_d, post_d = self.plot_distributions(
                r["prior_a"], r["prior_b"], r["likelihood"],
                r["posterior_a"], r["posterior_b"], drug, r["disease"]
            )

            likelihood_center = min(max(r["likelihood"], 0.01), 0.99)
            la = likelihood_center * 100
            lb = (1 - likelihood_center) * 100
            likelihood_pdf = beta.pdf(x, a=la + 1, b=lb + 1)
            likelihood_pdf /= trapezoid(likelihood_pdf, x)
            e_likelihood = trapezoid(x * likelihood_pdf, x)

            kl, mu_prior, mu_post, delta = self.compute_kl_and_mean_shift(prior_d, post_d, x)
            kl_metrics.append({
                "Disease": r["disease"],
                "KL Divergence": kl,
                "E[Prior]": mu_prior,
                "E[Likelihood]": round(e_likelihood, 4),
                "E[Posterior]": mu_post,
                "Δμ": delta,
                "Likelihood Score": round(r["likelihood"], 4),
                "Raw Prior": round(r["raw_prior"], 4),
                "Penalised Prior": round(r["penalised_prior"], 4),
                "Enhanced Prior": round(r["enhanced_prior"], 4)
            })

        print("\nKL Divergence and Mean Shift Summary:")
        print(pd.DataFrame(kl_metrics).to_string(index=False))


# ==== Run Predictor ====
if __name__ == "__main__":
    weights_dict = {
        "GraphDistanceToIndication": 0.1148,
        "RandomWalkScore": 0.247,
        "StructuralLikelihood": -0.1154,
        "PreferentialAttachment": -0.041,
        "KatzSimilarity": 1.6515
    }

    predictor = BayesianRepurposingPredictor(
        known_path="graph/graph_features_known.csv",
        unknown_path="graph/graph_features_unknown.csv",
        weights_dict=weights_dict
    )

    drugs = ["metformin"]
    diseases = [
                    "Chronic Kidney Disease",           
                    "Acute Kidney Injury",             
                    "renal insufficiency, chronic",         
                    "Polycystic Kidney Diseases",              
                    "Nephrotic Syndrome",              
                    "Glomerulonephritis",              
                    "glomerulonephritis, iga",                   
                    "glomerulonephritis, membranous",         
                    "carcinoma, renal cell"             
                ]


    for drug in drugs:
        predictor.evaluate_drug(drug, diseases)
