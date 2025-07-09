import os
import json
import pandas as pd

# === File paths ===
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
centrality_path = os.path.join(BASE_DIR, "graph", "drug_centrality.csv")
pairs_path = os.path.join(BASE_DIR, "processed_data", "condition_drug_pairs.json")
output_path = os.path.join(BASE_DIR, "graph", "training_dataset.csv")

# === Load centrality CSV ===
centrality_df = pd.read_csv(centrality_path)
centrality_df["Drug"] = centrality_df["Drug"].str.lower()
centrality_map = centrality_df.set_index("Drug").to_dict("index")

# === Load matched condition–drug pairs with phases ===
with open(pairs_path, "r", encoding="utf-8") as f:
    pair_data = json.load(f)

# === Helper: Extract highest trial phase (only 2, 3, or 4) ===
def get_labeled_phase(phases):
    if not phases:
        return None
    cleaned = [p.strip().upper().replace(" ", "") for p in phases]
    if any("PHASE4" in p or "PHASE3" in p for p in cleaned):
        return 1  # considered successful
    if any("PHASE2" in p for p in cleaned):
        return 0  # not yet successful
    return None  # exclude if no phase 2–4

# === Generate dataset ===
rows = []
for entry in pair_data:
    drug = entry["intervention"].lower()
    disease = entry["condition"].lower()
    phases = entry.get("phases", [])
    label = get_labeled_phase(phases)

    if label is None:
        continue  # exclude pairs not in Phase 2–4

    centrality = centrality_map.get(drug)
    if centrality is None:
        continue  # skip if centrality not found

    rows.append({
        "Drug": drug,
        "Disease": disease,
        "DegreeCentrality": centrality["DegreeCentrality"],
        "EigenvectorCentrality": centrality["EigenvectorCentrality"],
        "BetweennessCentrality": centrality["BetweennessCentrality"],
        "ClusteringCoefficient": centrality["ClusteringCoefficient"],
        "RandomWalkCentrality": centrality["RandomWalkCentrality"],
        "PhaseSuccess": label
    })

# === Save to CSV ===
df = pd.DataFrame(rows)
df.to_csv(output_path, index=False)
print(f"\n Training dataset saved to:\n{output_path}")
print("\n Sample rows:\n", df.head(8))
