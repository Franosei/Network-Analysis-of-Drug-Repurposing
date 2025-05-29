import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pubmed_utils import pubmed_count_requests
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.integrate import trapezoid

# --- Disease and Drug Class Hierarchies ---
disease_classes = {
    "alzheimer's disease": "neurological",
    "asthma": "respiratory",
    "cancer": "oncological",
    "depression": "psychiatric",
    "diabetes": "metabolic",
    "epilepsy": "neurological",
    "hiv aids": "infectious",
    "hypertension": "cardiovascular",
    "obesity": "metabolic",
    "rheumatoid arthritis": "autoimmune"
}

drug_classes = {
    "dexamethasone": "steroid",
    "rituximab": "monoclonal_antibody"
}

# --- Data Loading ---
def load_trial_links(data_dir="data"):
    """Load all drug-disease trial relationships from files."""
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
    """Load precomputed centrality scores for drugs."""
    df = pd.read_csv(csv_path)
    return df.set_index("Drug").to_dict("index")

# --- Prior Computations ---
def compute_group_prior(drug_class, disease_class):
    """Compute class-level (group) prior using literature counts."""
    try:
        joint = pubmed_count_requests(drug_class, disease_class)
        d = pubmed_count_requests(drug_class, "")
        c = pubmed_count_requests(disease_class, "")
        return joint / (d + c + 1e-5)
    except Exception:
        return 0.0

def compute_local_prior(drug, disease):
    """Compute specific drug-disease local prior from PubMed."""
    try:
        joint = pubmed_count_requests(drug, disease)
        d = pubmed_count_requests(drug, "")
        c = pubmed_count_requests(disease, "")
        return joint / (d + c + 1e-5)
    except Exception:
        return 0.0

def hierarchical_prior(drug, disease, lam=0.7):
    """Blend local and group prior based on lambda weight."""
    local = compute_local_prior(drug, disease)
    d_class = drug_classes.get(drug, drug)
    c_class = disease_classes.get(disease, disease)
    group = compute_group_prior(d_class, c_class)
    return lam * local + (1 - lam) * group

# --- Bayesian Posterior ---
def compute_posterior(drug, disease, centrality_scores, lam=0.7):
    """Compute final posterior score using prior and graph likelihood."""
    prior = hierarchical_prior(drug, disease, lam)
    scores = centrality_scores.get(drug.lower(), {})
    degree = scores.get("DegreeCentrality", 0)
    eigen = scores.get("EigenvectorCentrality", 0)
    between = scores.get("BetweennessCentrality", 0)
    likelihood = 1 + 2 * degree + 3 * eigen + 2 * between
    return prior * likelihood, prior

# --- KL Divergence ---
def compute_kl_divergence(prior, posterior):
    """Compute KL divergence and mean shift between two Beta distributions."""
    x = np.linspace(0.001, 0.999, 500)
    a1, b1 = prior * 100 + 1, (1 - prior) * 100 + 1
    a2, b2 = posterior * 100 + 1, (1 - posterior) * 100 + 1

    p1 = beta.pdf(x, a=a1, b=b1)
    p2 = beta.pdf(x, a=a2, b=b2)
    p1 = p1 / trapezoid(p1, x)
    p2 = p2 / trapezoid(p2, x)

    kl = trapezoid(rel_entr(p2, p1), x)
    mu1 = trapezoid(x * p1, x)
    mu2 = trapezoid(x * p2, x)
    return round(kl, 4), round(mu1, 4), round(mu2, 4), round(mu2 - mu1, 4)

# --- Prediction Loop ---
def predict_with_hierarchy(drug, diseases, trial_links, centrality_scores, lam=0.7):
    """Run prediction for a drug against all diseases not seen in trials."""
    results = []
    for disease in diseases:
        if (drug, disease) not in trial_links:
            posterior, prior = compute_posterior(drug, disease, centrality_scores, lam)
            results.append({
                "disease": disease,
                "posterior": posterior,
                "prior": prior
            })

    total = sum(r["posterior"] for r in results)
    for r in results:
        r["probability"] = round(r["posterior"] / total, 6) if total > 0 else 0.0

    return sorted(results, key=lambda x: x["probability"], reverse=True)

# --- Visualisation ---
def plot_hierarchy_prior_distributions(disease_classes):
    """Plot prior distributions grouped by disease class."""
    class_counts = defaultdict(list)

    for disease, cls in disease_classes.items():
        try:
            joint = pubmed_count_requests(cls, disease)
            d = pubmed_count_requests(cls, "")
            c = pubmed_count_requests(disease, "")
            prior = joint / (d + c + 1e-5)
            class_counts[cls].append(prior)
        except Exception:
            continue

    fig, ax = plt.subplots(figsize=(10, 6))
    for cls, priors in class_counts.items():
        x = np.linspace(0, 1, 500)
        alpha = np.mean(priors) * 50 + 1
        beta_param = 50
        y = beta.pdf(x, a=alpha, b=beta_param)
        y /= y.max()
        ax.plot(x, y, label=cls.title(), linewidth=2)

    ax.set_title("Distribution of Hierarchical Priors by Disease Class", fontsize=14)
    ax.set_xlabel("Prior Probability θ")
    ax.set_ylabel("Normalised Density")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

# --- Main ---
if __name__ == "__main__":
    print("Running Hierarchical Bayesian Inference and Visualisation...")
    trial_links, diseases, drugs = load_trial_links("data")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    selected_drug = "dexamethasone"
    predictions = predict_with_hierarchy(
        selected_drug.lower(),
        diseases,
        trial_links,
        centrality_scores,
        lam=0.6
    )

    print(f"\nHierarchical Inference for: {selected_drug.title()}")
    for item in predictions[:10]:
        print(f"{selected_drug.title()} → {item['disease'].title()} — P: {item['probability']}")

    print("\nGenerating Hierarchical Prior Distribution Plot...")
    plot_hierarchy_prior_distributions(disease_classes)

    # --- KL Divergence Report ---
    print("\nKL Divergence and Mean Shift (Top Predictions):")
    kl_records = []
    for item in predictions[:4]:
        kl, mu1, mu2, delta = compute_kl_divergence(item["prior"], item["posterior"])
        kl_records.append({
            "Disease": item["disease"].title(),
            "KL Divergence": kl,
            "E[Prior]": mu1,
            "E[Posterior]": mu2,
            "Δμ": delta
        })

    df = pd.DataFrame(kl_records)
    print(df.to_string(index=False))
