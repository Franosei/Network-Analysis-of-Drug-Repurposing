import os
import json
import re
import time
import random
from typing import Any, Dict, List, Optional, Tuple

import requests


class ClinicalTrialFetcher:
    """
    Fetch, filter, and save unique clinical trial data from ClinicalTrials.gov v2 API.

    Key behavior:
    - Fetches page-by-page until nextPageToken is exhausted (i.e., retrieves all available studies
      that match query.term).
    - Filters to eligible trials:
        (1) overallStatus in valid_statuses
        (2) phases intersect valid_phases (Phase 2-4)
        (3) has DRUG intervention
        (4) trial "date year" <= cutoff_year (default: 2020)
    - Global de-duplication across therapeutic areas (no repeated NCT IDs).
    """

    def __init__(
        self,
        output_dir: str = "data",
        page_size: int = 100,
        max_retries: int = 7,
        base_backoff_s: float = 1.0,
        jitter_s: float = 0.5,
        timeout_s: float = 30.0,
        cutoff_year: int = 2020,
        include_missing_dates: bool = False,
    ):
        self.base_url = "https://clinicaltrials.gov/api/v2/studies"
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", output_dir))
        os.makedirs(self.output_dir, exist_ok=True)

        self.page_size = int(page_size)

        self.max_retries = int(max_retries)
        self.base_backoff_s = float(base_backoff_s)
        self.jitter_s = float(jitter_s)
        self.timeout_s = float(timeout_s)

        self.cutoff_year = int(cutoff_year)
        self.include_missing_dates = bool(include_missing_dates)

        self.valid_statuses = {
            "RECRUITING",
            "ACTIVE_NOT_RECRUITING",
            "ENROLLING_BY_INVITATION",
            "COMPLETED",
            "WITHDRAWN",
            "TERMINATED",
        }
        self.valid_phases = {"PHASE2", "PHASE3", "PHASE4"}
        self.valid_intervention_type = "DRUG"

        self.seen_nct_ids = set()

        self.session = requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": "ClinicalTrialFetcher/1.0",
            }
        )

    # -----------------------------
    # Public API

    def fetch_trials_paginated(
        self,
        condition: str,
        max_unique_trials: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Fetch trials page-by-page for a condition until:
        - nextPageToken is exhausted (default, fetch all), OR
        - max_unique_trials is reached (if provided).

        Applies eligibility filtering + cutoff_year filtering locally.
        """
        print(f"\n[Fetching]: {condition}")

        unique_trials: List[Dict[str, Any]] = []
        next_page_token: Optional[str] = None

        while True:
            params = {
                "query.term": condition,
                "pageSize": self.page_size,
            }
            if next_page_token:
                params["pageToken"] = next_page_token

            result = self._get_with_retries(self.base_url, params=params)
            if result is None:
                break

            studies = result.get("studies", []) or []
            next_page_token = result.get("nextPageToken")

            for study in studies:
                if not self.is_eligible_trial(study):
                    continue

                # Enforce cutoff year (<= 2020) to avoid leakage into 2021–2025 test set
                if not self.is_within_cutoff(study):
                    continue

                nct_id = (
                    study.get("protocolSection", {})
                    .get("identificationModule", {})
                    .get("nctId", "")
                )
                if not nct_id:
                    continue

                if nct_id in self.seen_nct_ids:
                    continue

                unique_trials.append(self.extract_trial_data(study))
                self.seen_nct_ids.add(nct_id)

                if max_unique_trials is not None and len(unique_trials) >= int(max_unique_trials):
                    print(f"[{condition}] Reached cap: {max_unique_trials} unique eligible trials.")
                    return unique_trials

            time.sleep(random.uniform(0.2, 0.6))

            if not next_page_token:
                break

        print(f"[{condition}] Found {len(unique_trials)} unique eligible trials (<= {self.cutoff_year}).")
        return unique_trials

    def run_batch(self, therapeutic_areas: List[str], trials_per_area: Optional[int] = None) -> None:
        """
        If trials_per_area is None => fetch ALL eligible trials per area (subject to cutoff_year).
        If trials_per_area is an int => stop after that many unique eligible trials per area.
        """
        for area in therapeutic_areas:
            trials = self.fetch_trials_paginated(area, max_unique_trials=trials_per_area)
            if trials:
                filename = self.sanitize_filename(area) + f"_upto_{self.cutoff_year}.json"
                self.save_to_json(trials, filename)
            else:
                print(f"[{area}] No eligible trials found within cutoff year <= {self.cutoff_year}.")

    # -----------------------------
    # Filtering logic
    # -----------------------------
    def is_eligible_trial(self, study: Dict[str, Any]) -> bool:
        protocol = study.get("protocolSection", {}) or {}

        status = str(protocol.get("statusModule", {}).get("overallStatus", "")).upper()
        phases = {str(p).upper() for p in (protocol.get("designModule", {}).get("phases", []) or [])}
        interventions = protocol.get("armsInterventionsModule", {}).get("interventions", []) or []

        has_valid_intervention = any(
            str(i.get("type", "")).upper() == self.valid_intervention_type for i in interventions
        )

        return (
            status in self.valid_statuses
            and bool(phases.intersection(self.valid_phases))
            and has_valid_intervention
        )

    def is_within_cutoff(self, study: Dict[str, Any]) -> bool:
        """
        Returns True if the trial's reference year <= cutoff_year.

        Reference year selection (in order):
        1) startDateStruct.date
        2) completionDateStruct.date
        3) firstPostedDateStruct.date (if present)
        If none is parseable:
          - include if include_missing_dates=True
          - else exclude (conservative to avoid leakage)
        """
        year, source = self._extract_reference_year(study)
        if year is None:
            return self.include_missing_dates
        return year <= self.cutoff_year

    def _extract_reference_year(self, study: Dict[str, Any]) -> Tuple[Optional[int], Optional[str]]:
        protocol = study.get("protocolSection", {}) or {}
        status_mod = protocol.get("statusModule", {}) or {}

        candidates = [
            (status_mod.get("startDateStruct", {}) or {}).get("date", ""),
            (status_mod.get("completionDateStruct", {}) or {}).get("date", ""),
            (status_mod.get("firstPostedDateStruct", {}) or {}).get("date", ""),
        ]
        sources = ["startDate", "completionDate", "firstPostedDate"]

        for date_str, src in zip(candidates, sources):
            y = self._parse_year(date_str)
            if y is not None:
                return y, src

        return None, None

    @staticmethod
    def _parse_year(date_str: str) -> Optional[int]:
        """
        Extract a 4-digit year from ClinicalTrials.gov date strings.
        Handles common formats like:
          - "2019-04-12"
          - "April 2019"
          - "2019"
        """
        if not date_str:
            return None

        s = str(date_str).strip()
        # Direct year
        m = re.search(r"\b(19\d{2}|20\d{2})\b", s)
        if not m:
            return None

        try:
            return int(m.group(1))
        except Exception:
            return None

    # -----------------------------
    # Extraction / persistence
    # -----------------------------
    def extract_trial_data(self, study: Dict[str, Any]) -> Dict[str, Any]:
        protocol = study.get("protocolSection", {}) or {}
        interventions_raw = protocol.get("armsInterventionsModule", {}).get("interventions", []) or []

        interventions = [
            i.get("name", "")
            for i in interventions_raw
            if str(i.get("type", "")).upper() == self.valid_intervention_type
        ]
        intervention_types = [
            i.get("type", "")
            for i in interventions_raw
            if str(i.get("type", "")).upper() == self.valid_intervention_type
        ]

        return {
            "nctId": protocol.get("identificationModule", {}).get("nctId", ""),
            "title": protocol.get("identificationModule", {}).get("officialTitle", ""),
            "briefTitle": protocol.get("identificationModule", {}).get("briefTitle", ""),
            "status": protocol.get("statusModule", {}).get("overallStatus", ""),
            "startDate": protocol.get("statusModule", {}).get("startDateStruct", {}).get("date", ""),
            "completionDate": protocol.get("statusModule", {}).get("completionDateStruct", {}).get("date", ""),
            "phases": protocol.get("designModule", {}).get("phases", []) or [],
            "interventions": interventions,
            "interventionTypes": intervention_types,
            "conditions": protocol.get("conditionsModule", {}).get("conditions", []) or [],
        }

    @staticmethod
    def sanitize_filename(name: str) -> str:
        return re.sub(r"[^\w\-_\.]", "_", name.lower().replace(" ", "_"))

    def save_to_json(self, data: List[Dict[str, Any]], filename: str) -> None:
        output_path = os.path.join(self.output_dir, filename)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        print(f"[Saved] {len(data)} trials → {output_path}")

    # -----------------------------
    # HTTP with retries
    # -----------------------------
    def _get_with_retries(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        for attempt in range(1, self.max_retries + 1):
            try:
                resp = self.session.get(url, params=params, timeout=self.timeout_s)

                if resp.status_code == 429 or (500 <= resp.status_code < 600):
                    backoff = self.base_backoff_s * (2 ** (attempt - 1))
                    backoff += random.uniform(0, self.jitter_s)
                    print(
                        f"[Retryable HTTP {resp.status_code}] Backing off {backoff:.2f}s "
                        f"(attempt {attempt}/{self.max_retries})"
                    )
                    time.sleep(backoff)
                    continue

                resp.raise_for_status()
                return resp.json()

            except requests.exceptions.RequestException as e:
                backoff = self.base_backoff_s * (2 ** (attempt - 1))
                backoff += random.uniform(0, self.jitter_s)
                print(
                    f"[Request error] {e} | Backing off {backoff:.2f}s "
                    f"(attempt {attempt}/{self.max_retries})"
                )
                time.sleep(backoff)

        print("[Aborted] Too many retries.")
        return None


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
        "obesity",
    ]

    # Training set: only trials up to and including 2020
    fetcher = ClinicalTrialFetcher(
        output_dir="data",
        page_size=100,
        cutoff_year=2020,
        include_missing_dates=False,  # conservative: avoids leakage
    )

    # Fetch ALL eligible trials per area, subject to cutoff_year
    fetcher.run_batch(therapeutic_list, trials_per_area=None)
