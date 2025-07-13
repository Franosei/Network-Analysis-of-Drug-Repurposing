import os
import json
import re
import requests
import time
import random


class ClinicalTrialFetcher:
    """
    A class to fetch, filter, and save unique clinical trial data from ClinicalTrials.gov API.
    Ensures no trial is repeated across therapeutic areas.
    """

    def __init__(self, output_dir="data", min_trials=10):
        self.base_url = "https://clinicaltrials.gov/api/v2/studies"
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", output_dir))
        self.min_trials = min_trials
        os.makedirs(self.output_dir, exist_ok=True)

        self.valid_statuses = {
            "RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION",
            "COMPLETED", "WITHDRAWN", "TERMINATED"
        }
        self.valid_phases = {"PHASE2", "PHASE3", "PHASE4"}
        self.valid_intervention_type = "DRUG"
        self.seen_nct_ids = set()  # Global deduplication tracker

    def fetch_trials_paginated(self, condition, max_unique_trials=500):
        """
        Fetch trials page-by-page for a condition until enough unique, eligible trials are found.
        Includes retry logic for rate-limiting.
        """
        print(f"\n[Fetching]: {condition}")
        unique_trials = []
        next_page_token = None
        retries = 0
        max_retries = 5

        while len(unique_trials) < max_unique_trials:
            params = {
                "query.term": condition,
                "pageSize": 100,
            }
            if next_page_token:
                params["pageToken"] = next_page_token

            try:
                response = requests.get(self.base_url, params=params)
                if response.status_code == 429:
                    raise requests.exceptions.HTTPError("Rate limit exceeded", response=response)
                response.raise_for_status()
                result = response.json()
                studies = result.get("studies", [])
                next_page_token = result.get("nextPageToken")

                for study in studies:
                    if self.is_eligible_trial(study):
                        nct_id = study.get("protocolSection", {}).get("identificationModule", {}).get("nctId", "")
                        if nct_id and nct_id not in self.seen_nct_ids:
                            trial_info = self.extract_trial_data(study)
                            unique_trials.append(trial_info)
                            self.seen_nct_ids.add(nct_id)
                    if len(unique_trials) >= max_unique_trials:
                        break

                # Delay between requests to avoid triggering rate limit
                time.sleep(random.uniform(0.3, 0.8))
                retries = 0  # Reset retry count after a successful call

            except requests.exceptions.HTTPError as e:
                if response.status_code == 429:
                    wait_time = 2 ** retries
                    print(f"[Rate limited] Waiting {wait_time}s before retrying... (Attempt {retries+1})")
                    time.sleep(wait_time)
                    retries += 1
                    if retries > max_retries:
                        print("[Aborted] Too many retries due to rate limiting.")
                        break
                else:
                    print(f"[HTTP Error] {e}")
                    break
            except Exception as e:
                print(f"[Error] {e}")
                break

            if not next_page_token:
                break 

        print(f"[{condition}] Found {len(unique_trials)} unique eligible trials.")
        return unique_trials

    def extract_trial_data(self, study):
        protocol = study.get("protocolSection", {})
        interventions_raw = protocol.get("armsInterventionsModule", {}).get("interventions", [])

        interventions = [
            i.get("name", "") for i in interventions_raw
            if i.get("type", "").upper() == self.valid_intervention_type
        ]
        intervention_types = [
            i.get("type", "") for i in interventions_raw
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

    def is_eligible_trial(self, study):
        protocol = study.get("protocolSection", {})
        status = protocol.get("statusModule", {}).get("overallStatus", "").upper()
        phases = set(p.upper() for p in protocol.get("designModule", {}).get("phases", []))
        interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])

        has_valid_intervention = any(
            i.get("type", "").upper() == self.valid_intervention_type for i in interventions
        )

        return (
            status in self.valid_statuses and
            phases.intersection(self.valid_phases) and
            has_valid_intervention
        )

    def sanitize_filename(self, name):
        return re.sub(r"[^\w\-_\.]", "_", name.lower().replace(" ", "_"))

    def save_to_json(self, data, filename):
        output_path = os.path.join(self.output_dir, filename)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        print(f"[Saved] {len(data)} trials → {output_path}")

    def run_batch(self, therapeutic_areas, trials_per_area=1000):
        for area in therapeutic_areas:
            trials = self.fetch_trials_paginated(area, max_unique_trials=trials_per_area)
            if trials:
                filename = self.sanitize_filename(area) + ".json"
                self.save_to_json(trials, filename)

if __name__ == "__main__":
    therapeutic_list = [
        "cancer", "diabetes", "asthma", "Alzheimer's disease", "depression",
        "hypertension", "HIV/AIDS", "rheumatoid arthritis", "epilepsy", "obesity"
    ]

    fetcher = ClinicalTrialFetcher()
    fetcher.run_batch(therapeutic_list, trials_per_area=1000)
