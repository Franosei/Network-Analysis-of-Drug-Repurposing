import os
import json
import re
from difflib import get_close_matches
from typing import Set, Dict, List, Tuple
from xml.etree import ElementTree as ET

class ConditionDrugPairBuilder:
    def __init__(
        self,
        input_dir="data",
        output_dir="processed_data",
        output_filename="condition_drug_pairs.json",
        unmatched_filename="unmatched_pairs.json",
        mesh_path="mesh_data/desc2025.xml"
    ):
        self.input_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", input_dir))
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", output_dir))
        os.makedirs(self.output_dir, exist_ok=True)

        self.output_path = os.path.join(self.output_dir, output_filename)
        self.unmatched_path = os.path.join(self.output_dir, unmatched_filename)
        self.mesh_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", mesh_path))

        self.condition_term_map, self.drug_term_map = self.load_mesh_terms()

    def load_mesh_terms(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        condition_map = {}
        drug_map = {}

        tree = ET.parse(self.mesh_path)
        root = tree.getroot()

        for descriptor in root.findall("DescriptorRecord"):
            canonical = descriptor.findtext("DescriptorName/String").strip()
            tree_numbers = [el.text for el in descriptor.findall("TreeNumberList/TreeNumber")]

            is_condition = any(num.startswith("C") for num in tree_numbers)
            is_drug = any(num.startswith("D") for num in tree_numbers)
            if not (is_condition or is_drug):
                continue

            terms = set()
            terms.add(canonical.lower())
            for entry_term in descriptor.findall(".//TermList/Term/String"):
                terms.add(entry_term.text.strip().lower())

            if is_condition:
                for term in terms:
                    condition_map[term] = canonical
            if is_drug:
                for term in terms:
                    drug_map[term] = canonical

        return condition_map, drug_map

    def clean_text(self, text: str) -> str:
        """
        Clean up final output name by removing content in () or [] and punctuation.
        """
        return re.sub(r"\s*[\[(].*?[\])]\s*", "", text).strip().strip(",. ")

    def normalize_items(self, items: List[str]) -> List[str]:
        """
        Normalize raw input strings and return a flat list of cleaned items.
        """
        if not isinstance(items, list):
            return []

        result = []
        for item in items:
            parts = re.split(r'\s*(?:,|;|\+|\band\b|/)\s*', item)
            for part in parts:
                part = part.strip()
                if part and part.lower() != "placebo":
                    result.append(part)
        return result

    def match(self, term: str, mesh_map: Dict[str, str]) -> str:
        """
        Match using exact and fallback fuzzy logic.
        """
        key = term.lower()
        if key in mesh_map:
            return mesh_map[key]

        close = get_close_matches(key, mesh_map.keys(), n=1, cutoff=0.85)
        if close:
            return mesh_map[close[0]]
        return None

    def extract_pairs(self):
        rows = []
        unmatched = []

        files = [f for f in os.listdir(self.input_dir) if f.endswith(".json")]
        print(f"Processing {len(files)} JSON file(s) in '{self.input_dir}'...")

        for file in files:
            with open(os.path.join(self.input_dir, file), "r", encoding="utf-8") as f:
                trials = json.load(f)

            for trial in trials:
                conditions = self.normalize_items(trial.get("conditions", []))
                interventions = self.normalize_items(trial.get("interventions", []))

                # Safely extract phase info and format as list
                raw_phase = trial.get("phases", trial.get("phase", None))
                if isinstance(raw_phase, list):
                    phases = raw_phase
                elif isinstance(raw_phase, str):
                    phases = [raw_phase]
                else:
                    phases = []

                for cond in conditions:
                    cond_match = self.match(cond, self.condition_term_map)
                    if not cond_match:
                        unmatched.append({"condition": cond, "intervention": None, "reason": "unmatched condition"})
                        continue

                    for drug in interventions:
                        drug_match = self.match(drug, self.drug_term_map)
                        if not drug_match:
                            unmatched.append({"condition": cond_match, "intervention": drug, "reason": "unmatched intervention"})
                            continue

                        rows.append({
                            "condition": self.clean_text(cond_match),
                            "intervention": self.clean_text(drug_match),
                            "phases": phases
                        })

        print(f"Matched: {len(rows)}  Unmatched: {len(unmatched)}")
        return rows, unmatched

    def save_pairs(self, rows: List[Dict[str, str]]):
        with open(self.output_path, "w", encoding="utf-8") as f:
            json.dump(rows, f, indent=2)
        print(f"Saved matched pairs to: {self.output_path}")

    def save_unmatched(self, unmatched: List[Dict[str, str]]):
        with open(self.unmatched_path, "w", encoding="utf-8") as f:
            json.dump(unmatched, f, indent=2)
        print(f"Saved unmatched entries to: {self.unmatched_path}")


if __name__ == "__main__":
    builder = ConditionDrugPairBuilder()
    matched, unmatched = builder.extract_pairs()
    builder.save_pairs(matched)
    builder.save_unmatched(unmatched)
