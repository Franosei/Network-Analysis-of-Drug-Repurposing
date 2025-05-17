import os
import json
import pandas as pd

target_drug = "dexamethasone"
data_dir = "../data"
summary = []

for file in os.listdir(data_dir):
    if file.endswith(".json"):
        path = os.path.join(data_dir, file)
        with open(path, "r", encoding="utf-8") as f:
            trials = json.load(f)
        for trial in trials:
            interventions = [i.lower() for i in trial.get("interventions", [])]
            if target_drug in interventions:
                summary.append({
                    "Disease": file.replace(".json", "").replace("_", " "),
                    "Title": trial.get("briefTitle", ""),
                    "Phases": trial.get("phases", []),
                    "Status": trial.get("status", ""),
                    "NCT ID": trial.get("nctId", "")
                })


df = pd.DataFrame(summary)
print(df.head(10))
