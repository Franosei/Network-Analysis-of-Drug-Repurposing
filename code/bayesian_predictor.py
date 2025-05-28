import os
import json
import pandas as pd
from collections import defaultdict
from pubmed_utils import pubmed_count_requests


def load_trial_links(data_dir="data"):
    """Load all drug–disease pairs from trial files."""
    trial_links = defaultdict(set)
    all_diseases = set()
    all_drugs = set()

    for file in os.listdir(data_dir):
        if file.endswith(".json"):
            disease = file.replace(".json", "").replace("_", " ")
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
    """Load precomputed centrality scores."""
    df = pd.read_csv(csv_path)
    return df.set_index("Drug").to_dict("index")


def compute_posterior(drug, disease, centrality_scores):
    """Compute unnormalized posterior probability for drug–disease pair."""
    try:
        pubmed_joint = pubmed_count_requests(drug, disease)
        pubmed_d = pubmed_count_requests(drug, "")
        pubmed_c = pubmed_count_requests(disease, "")
        pubmed_prior = pubmed_joint / (pubmed_d + pubmed_c + 1e-5)
    except Exception:
        pubmed_prior = 0.0

    scores = centrality_scores.get(drug.lower(), {})
    degree = scores.get("DegreeCentrality", 0)
    eigen = scores.get("EigenvectorCentrality", 0)
    between = scores.get("BetweennessCentrality", 0)
    likelihood = 1 + 2 * degree + 3 * eigen + 2 * between

    return pubmed_prior * likelihood


def predict_unexplored_diseases_for_drug(drug, diseases, trial_links, centrality_scores, normalise=True):
    """Predict unexplored diseases for a given drug with optional posterior normalisation."""
    results = []
    
    for disease in diseases:
        if (drug, disease) not in trial_links:
            score = compute_posterior(drug, disease, centrality_scores)
            results.append({
                "disease": disease,
                "posterior": score
            })

    if normalise:
        total = sum(item["posterior"] for item in results)
        for item in results:
            item["probability"] = round(item["posterior"] / total, 6) if total > 0 else 0.0

    # Sort by probability if normalised, else by raw posterior
    sort_key = "probability" if normalise else "posterior"
    return sorted(results, key=lambda x: x[sort_key], reverse=True)


if __name__ == "__main__":
    print("Loading trial data and centrality scores...")
    trial_links, diseases, drugs = load_trial_links("data")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    selected_drug = "dexamethasone"  # Replace as needed

    print(f"\n Predicting unexplored diseases for: {selected_drug.title()}")
    predictions = predict_unexplored_diseases_for_drug(
        selected_drug.lower(),
        diseases,
        trial_links,
        centrality_scores,
        normalise=True
    )

    print(f"\n Normalised Bayesian Probabilities for {selected_drug.title()}:")
    for item in predictions[:10]:
        print(f"{selected_drug.title()} → {item['disease'].title()} — Probability: {item['probability']}")

