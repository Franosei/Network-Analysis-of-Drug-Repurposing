import os
import json
import re
import requests

class ClinicalTrialFetcher:
    def __init__(self, output_dir="data", min_trials=10):
        self.base_url = "https://clinicaltrials.gov/api/v2/studies"
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", output_dir))
        self.min_trials = min_trials
        os.makedirs(self.output_dir, exist_ok=True)

        self.valid_statuses = {
            "RECRUITING",
            "ACTIVE_NOT_RECRUITING",
            "ENROLLING_BY_INVITATION",
            "COMPLETED"
        }
        self.valid_phases = {"PHASE1", "PHASE2", "PHASE3", "PHASE4"}
        self.valid_intervention_type = "DRUG"

    def fetch_trials(self, condition, max_trials=500):
        """
        Fetch trials from ClinicalTrials.gov for a given condition.
        """
        params = {
            "query.term": condition,
            "pageSize": max_trials
        }
        response = requests.get(self.base_url, params=params)
        response.raise_for_status()
        return response.json().get("studies", [])

    def extract_trial_data(self, study):
        """
        Extract and structure relevant data from a single trial entry.
        """
        protocol = study.get("protocolSection", {})
        interventions_raw = protocol.get("armsInterventionsModule", {}).get("interventions", [])

        interventions = [
            i.get("name", "")
            for i in interventions_raw
            if i.get("type", "").upper() == self.valid_intervention_type
        ]
        intervention_types = [
            i.get("type", "")
            for i in interventions_raw
            if i.get("type", "").upper() == self.valid_intervention_type
        ]

        return {
            "nctId": protocol.get("identificationModule", {}).get("nctId", ""),
            "title": protocol.get("identificationModule", {}).get("officialTitle", ""),
            "briefTitle": protocol.get("identificationModule", {}).get("briefTitle", ""),
            "status": protocol.get("statusModule", {}).get("overallStatus", ""),
            "startDate": protocol.get("statusModule", {}).get("startDateStruct", {}).get("date", ""),
            "completionDate": protocol.get("statusModule", {}).get("completionDateStruct", {}).get("date", ""),
            "phases": protocol.get("designModule", {}).get("phases", []),
            "interventions": interventions,
            "interventionTypes": intervention_types,
            "conditions": protocol.get("conditionsModule", {}).get("conditions", [])
        }

    def sanitize_filename(self, name):
        """
        Sanitize filename to remove or replace unsafe characters.
        """
        return re.sub(r"[^\w\-_\.]", "_", name.lower().replace(" ", "_"))

    def save_to_json(self, data, filename):
        """
        Save structured trial data to a JSON file.
        """
        output_path = os.path.join(self.output_dir, filename)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        print(f"Saved {len(data)} trials to {output_path}")

    def is_eligible_trial(self, study):
        """
        Determine if a study meets status, phase, and intervention type criteria.
        """
        protocol = study.get("protocolSection", {})
        status = protocol.get("statusModule", {}).get("overallStatus", "").upper()
        phases = set(p.upper() for p in protocol.get("designModule", {}).get("phases", []))
        interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])

        has_valid_intervention = any(
            i.get("type", "").upper() == self.valid_intervention_type
            for i in interventions
        )

        return (
            status in self.valid_statuses and
            phases.intersection(self.valid_phases) and
            has_valid_intervention
        )

    def process_area(self, therapeutic_area):
        """
        Fetch, filter, and save clinical trials for a specific therapeutic area.
        """
        print(f"\nProcessing: {therapeutic_area}")
        trials = self.fetch_trials(therapeutic_area)
        eligible_trials = []

        for study in trials:
            try:
                if self.is_eligible_trial(study):
                    trial_info = self.extract_trial_data(study)
                    eligible_trials.append(trial_info)
            except Exception as e:
                print(f"Skipping trial due to error: {e}")

        if not eligible_trials:
            print(f"No eligible trials found for {therapeutic_area}")
            return

        safe_filename = self.sanitize_filename(therapeutic_area) + ".json"
        self.save_to_json(eligible_trials, safe_filename)

    def run_batch(self, therapeutic_areas):
        """
        Process multiple therapeutic areas sequentially.
        """
        for area in therapeutic_areas:
            self.process_area(area)


if __name__ == "__main__":
    therapeutic_list = [
        "cancer",
        "diabetes",
        "asthma",
        "Alzheimer's disease",
        "depression",
        "hypertension",
        "HIV/AIDS",
        "rheumatoid arthritis",
        "epilepsy",
        "obesity"
    ]

    fetcher = ClinicalTrialFetcher()
    fetcher.run_batch(therapeutic_list)
