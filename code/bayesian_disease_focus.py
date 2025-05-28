import os
import json
import pandas as pd
from collections import defaultdict
from pubmed_utils import pubmed_count_requests
from bayesian_predictor import compute_posterior  # Assumes compute_posterior is defined in bayesian_predictor.py


def load_trial_links(data_dir="data"):
    """Load all drug–disease pairs from trial files."""
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
                        drug = drug.lower()
                        all_drugs.add(drug)
                        trial_links[(drug, disease)].add(trial["nctId"])

    return trial_links, list(all_diseases), list(all_drugs)


def load_centrality_scores(csv_path="graph/drug_centrality.csv"):
    """Load precomputed centrality scores."""
    df = pd.read_csv(csv_path)
    return df.set_index("Drug").to_dict("index")


def predict_unexplored_drugs_for_disease(disease, drugs, trial_links, centrality_scores, normalise=True):
    """Compare multiple drugs for repurposing against a single disease."""
    results = []

    for drug in drugs:
        if (drug, disease) not in trial_links:
            score = compute_posterior(drug, disease, centrality_scores)
            results.append({
                "drug": drug,
                "posterior": score
            })

    if normalise:
        total = sum(item["posterior"] for item in results)
        for item in results:
            item["probability"] = round(item["posterior"] / total, 6) if total > 0 else 0.0

    sort_key = "probability" if normalise else "posterior"
    return sorted(results, key=lambda x: x[sort_key], reverse=True)


if __name__ == "__main__":
    print("Loading trial data and centrality scores...")
    trial_links, diseases, drugs = load_trial_links("data")
    centrality_scores = load_centrality_scores("graph/drug_centrality.csv")

    selected_disease = "obesity"

    print(f"\nComparing candidate drugs for repurposing in: {selected_disease.title()}")
    predictions = predict_unexplored_drugs_for_disease(
        selected_disease.lower(),
        drugs,
        trial_links,
        centrality_scores,
        normalise=True
    )

    print(f"\nNormalised Bayesian Probabilities for Drug Repurposing in {selected_disease.title()}:")
    for item in predictions[:10]:
        print(f"{item['drug'].title()} → {selected_disease.title()} — Probability: {item['probability']}")
