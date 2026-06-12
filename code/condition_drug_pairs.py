import os
import json
import re
from difflib import SequenceMatcher, get_close_matches
from typing import Dict, List, Tuple, Optional
from xml.etree import ElementTree as ET


class ConditionDrugPairBuilder:
    """
    Process clinical trial JSON files to extract condition–intervention pairs,
    normalize them against MeSH (DescriptorRecord terms + Entry Terms),
    and record phases and status.

    This version is domain-agnostic:
      - no hard-coded disease/drug examples
      - normalization rules are generic patterns that improve MeSH matchability
      - matching uses exact -> fuzzy -> token-guided scoring
    """

    def __init__(
        self,
        input_dir: str = "data",
        output_dir: str = "processed_data",
        output_filename: str = "condition_drug_pairs.json",
        unmatched_filename: str = "unmatched_pairs.json",
        mesh_path: str = "mesh_data/desc2026.xml",
        fuzzy_cutoff: float = 0.80,
        token_jaccard_min: float = 0.60,
        include_placebo: bool = False,
        keep_unmatched_debug_fields: bool = True,
    ):
        self.input_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", input_dir))
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", output_dir))
        os.makedirs(self.output_dir, exist_ok=True)

        self.output_path = os.path.join(self.output_dir, output_filename)
        self.unmatched_path = os.path.join(self.output_dir, unmatched_filename)
        self.mesh_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", mesh_path))

        self.fuzzy_cutoff = float(fuzzy_cutoff)
        self.token_jaccard_min = float(token_jaccard_min)
        self.include_placebo = bool(include_placebo)
        self.keep_unmatched_debug_fields = bool(keep_unmatched_debug_fields)

        # Load MeSH maps keyed by normalized strings
        self.condition_term_map, self.drug_term_map = self.load_mesh_terms()

        # Token indexes speed up token-guided matching
        self._cond_token_index = self._build_token_index(self.condition_term_map)
        self._drug_token_index = self._build_token_index(self.drug_term_map)

    # ------------------------------------------------------------------
    # MeSH loading
    # ------------------------------------------------------------------
    def load_mesh_terms(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Load MeSH terms from XML and build:
          - condition_term_map: normalized_term -> canonical_descriptor
          - drug_term_map:      normalized_term -> canonical_descriptor

        Classification rule (as in your original):
          - Condition: TreeNumber starts with "C"
          - Drug:      TreeNumber starts with "D"
        """
        condition_map: Dict[str, str] = {}
        drug_map: Dict[str, str] = {}

        tree = ET.parse(self.mesh_path)
        root = tree.getroot()

        for descriptor in root.findall("DescriptorRecord"):
            canonical = (descriptor.findtext("DescriptorName/String") or "").strip()
            if not canonical:
                continue

            tree_numbers = [
                el.text for el in descriptor.findall("TreeNumberList/TreeNumber")
                if el is not None and el.text
            ]

            is_condition = any(num.startswith("C") for num in tree_numbers)
            is_drug = any(num.startswith("D") for num in tree_numbers)
            if not (is_condition or is_drug):
                continue

            # Collect canonical + entry terms; store in normalized form
            terms = set()
            terms.add(self.normalize_for_match(canonical))

            for entry_term in descriptor.findall(".//TermList/Term/String"):
                if entry_term is not None and entry_term.text:
                    terms.add(self.normalize_for_match(entry_term.text))

            if is_condition:
                for t in terms:
                    if t:
                        condition_map[t] = canonical
            if is_drug:
                for t in terms:
                    if t:
                        drug_map[t] = canonical

        return condition_map, drug_map

    # ------------------------------------------------------------------
    # Normalization (generic, no examples)
    # ------------------------------------------------------------------
    @staticmethod
    def _strip_brackets(text: str) -> str:
        # Remove bracketed qualifiers: (...) and [...]
        return re.sub(r"\s*[\[(].*?[\])]\s*", " ", text)

    @staticmethod
    def _basic_cleanup(text: str) -> str:
        text = text.lower().strip()
        text = text.replace("&", " and ")
        text = text.replace("’", "'")
        # normalize separators to spaces
        text = re.sub(r"[_\-]+", " ", text)
        text = re.sub(r"\s+", " ", text).strip()
        return text

    @staticmethod
    def _drop_common_qualifiers(text: str) -> str:
        """
        Remove generic qualifiers frequently attached to biomedical entities
        that degrade ontology matching but usually do not define the core entity.
        """
        # staging / grading patterns (generic roman/arabic)
        text = re.sub(r"\b(stage|grade)\s*(?:[0-9]+|[ivx]+)\b", " ", text)

        # common clinical qualifiers (generic list; not disease/drug specific)
        text = re.sub(
            r"\b(advanced|metastatic|recurrent|relapsed|resistant|refractory|"
            r"unresectable|chronic|acute|severe|mild|moderate|progressive|"
            r"localized|locally advanced|disseminated)\b",
            " ",
            text,
        )

        # remove common study-context tokens
        text = re.sub(
            r"\b(randomized|open label|double blind|single blind|placebo controlled|"
            r"controlled|pilot|feasibility|observational)\b",
            " ",
            text,
        )

        # remove measurement-like tokens and symbols often appended
        text = re.sub(r"\b(positive|negative)\b", " ", text)
        text = re.sub(r"\b\d+(\.\d+)?\s*(mg|g|mcg|µg|ug|ml|l|iu|%)\b", " ", text)

        return re.sub(r"\s+", " ", text).strip()

    @staticmethod
    def _normalize_numbers(text: str) -> str:
        """
        Normalize common numeric variants without referencing particular entities.
        """
        # roman numerals for types/classes
        text = re.sub(r"\btype\s+ii\b", "type 2", text)
        text = re.sub(r"\btype\s+iii\b", "type 3", text)
        text = re.sub(r"\btype\s+iv\b", "type 4", text)
        text = re.sub(r"\btype\s+i\b", "type 1", text)

        # normalize whitespace
        return re.sub(r"\s+", " ", text).strip()

    @staticmethod
    def _remove_punctuation(text: str) -> str:
        # Keep alphanumerics/spaces only; drop remaining punctuation
        text = re.sub(r"[^a-z0-9\s]+", " ", text)
        return re.sub(r"\s+", " ", text).strip()

    def normalize_for_match(self, text: str) -> str:
        """
        Canonical normalization used both for:
          - MeSH term keys
          - incoming trial strings
        """
        if not text:
            return ""

        t = self._strip_brackets(text)
        t = self._basic_cleanup(t)
        t = self._drop_common_qualifiers(t)
        t = self._normalize_numbers(t)
        t = self._remove_punctuation(t)

        return t

    @staticmethod
    def clean_output_name(text: str) -> str:
        # For final output display: remove bracketed noise if any remains
        return re.sub(r"\s*[\[(].*?[\])]\s*", "", text).strip().strip(",. ")

    # ------------------------------------------------------------------
    # Trial string splitting/tokenization
    # ------------------------------------------------------------------
    def normalize_items(self, items: List[str]) -> List[str]:
        if not isinstance(items, list):
            return []

        result: List[str] = []
        for item in items:
            if not item:
                continue

            parts = re.split(r"\s*(?:,|;|\+|/|\band\b)\s*", str(item))
            for part in parts:
                part = part.strip()
                if not part:
                    continue
                if not self.include_placebo and part.lower() == "placebo":
                    continue
                result.append(part)

        return result

    @staticmethod
    def _tokenize(text: str) -> List[str]:
        return [tok for tok in text.split() if tok]

    @staticmethod
    def _jaccard(a: List[str], b: List[str]) -> float:
        sa, sb = set(a), set(b)
        if not sa or not sb:
            return 0.0
        return len(sa & sb) / len(sa | sb)

    @staticmethod
    def _seq_ratio(a: str, b: str) -> float:
        return SequenceMatcher(None, a, b).ratio()

    def _build_token_index(self, mesh_map: Dict[str, str]) -> Dict[str, List[str]]:
        """
        token -> list of mesh keys that contain that token
        """
        idx: Dict[str, List[str]] = {}
        for k in mesh_map.keys():
            toks = self._tokenize(k)
            for tok in toks:
                idx.setdefault(tok, []).append(k)
        return idx

    # ------------------------------------------------------------------
    # Matching: exact -> difflib close -> token-guided scoring
    # ------------------------------------------------------------------
    def match_term(
        self,
        raw_term: str,
        mesh_map: Dict[str, str],
        token_index: Dict[str, List[str]],
    ) -> Tuple[Optional[str], Optional[Dict[str, object]]]:
        """
        Returns:
          - canonical MeSH descriptor name (or None)
          - debug dict with match metadata (or None)
        """
        norm = self.normalize_for_match(raw_term)
        if not norm:
            return None, {"method": "empty", "normalized": norm}

        # 1) Exact
        if norm in mesh_map:
            return mesh_map[norm], {
                "method": "exact",
                "normalized": norm,
                "matched_key": norm,
                "score": 1.0,
                "sequence_ratio": 1.0,
                "token_jaccard_score": 1.0,
            }

        # 2) Fuzzy over all keys
        close = get_close_matches(norm, mesh_map.keys(), n=1, cutoff=self.fuzzy_cutoff)
        if close:
            norm_tokens = self._tokenize(norm)
            close_tokens = self._tokenize(close[0])
            seq_ratio = self._seq_ratio(norm, close[0])
            token_jaccard = self._jaccard(norm_tokens, close_tokens)
            return mesh_map[close[0]], {
                "method": "fuzzy",
                "normalized": norm,
                "matched_key": close[0],
                "score": round(seq_ratio, 4),
                "sequence_ratio": round(seq_ratio, 4),
                "token_jaccard_score": round(token_jaccard, 4),
            }

        # 3) Token-guided candidate narrowing + scoring
        toks = self._tokenize(norm)
        if not toks:
            return None, {"method": "no_tokens", "normalized": norm}

        candidates: List[str] = []
        seen = set()
        for tok in toks:
            for cand in token_index.get(tok, []):
                if cand not in seen:
                    candidates.append(cand)
                    seen.add(cand)

        if not candidates:
            return None, {"method": "no_candidates", "normalized": norm}

        best_key = None
        best_score = 0.0
        best_jaccard = 0.0
        best_ratio = 0.0

        for cand in candidates:
            cand_toks = self._tokenize(cand)
            j = self._jaccard(toks, cand_toks)
            if j < self.token_jaccard_min:
                continue

            r = self._seq_ratio(norm, cand)
            score = 0.65 * j + 0.35 * r

            if score > best_score:
                best_score = score
                best_key = cand
                best_jaccard = j
                best_ratio = r

        if best_key:
            return mesh_map[best_key], {
                "method": "token_score",
                "normalized": norm,
                "matched_key": best_key,
                "score": round(best_score, 4),
                "sequence_ratio": round(best_ratio, 4),
                "token_jaccard_score": round(best_jaccard, 4),
            }

        return None, {"method": "token_score_failed", "normalized": norm}

    # ------------------------------------------------------------------
    # Core extraction
    # ------------------------------------------------------------------
    def extract_pairs(self) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
        rows: List[Dict[str, object]] = []
        unmatched: List[Dict[str, object]] = []

        files = [f for f in os.listdir(self.input_dir) if f.endswith(".json")]
        print(f"Processing {len(files)} JSON file(s) in '{self.input_dir}'...")

        for file in files:
            path = os.path.join(self.input_dir, file)
            with open(path, "r", encoding="utf-8") as f:
                trials = json.load(f)

            for trial in trials:
                conditions = self.normalize_items(trial.get("conditions", []))
                interventions = self.normalize_items(trial.get("interventions", []))

                raw_phase = trial.get("phases", trial.get("phase", None))
                if isinstance(raw_phase, list):
                    phases = raw_phase
                elif isinstance(raw_phase, str):
                    phases = [raw_phase]
                else:
                    phases = []

                status = str(trial.get("status", "")).strip()

                for cond in conditions:
                    cond_match, cond_dbg = self.match_term(
                        raw_term=cond,
                        mesh_map=self.condition_term_map,
                        token_index=self._cond_token_index,
                    )
                    if not cond_match:
                        entry = {
                            "raw_condition": cond,
                            "raw_intervention": None,
                            "reason": "unmatched condition",
                        }
                        if self.keep_unmatched_debug_fields and cond_dbg:
                            entry.update({"debug": cond_dbg})
                        unmatched.append(entry)
                        continue

                    for drug in interventions:
                        drug_match, drug_dbg = self.match_term(
                            raw_term=drug,
                            mesh_map=self.drug_term_map,
                            token_index=self._drug_token_index,
                        )
                        if not drug_match:
                            entry = {
                                "raw_condition": cond_match,
                                "raw_intervention": drug,
                                "reason": "unmatched intervention",
                            }
                            if self.keep_unmatched_debug_fields and drug_dbg:
                                entry.update({"debug": drug_dbg})
                            unmatched.append(entry)
                            continue

                        rows.append({
                            "condition": self.clean_output_name(cond_match),
                            "intervention": self.clean_output_name(drug_match),
                            "phases": phases,
                            "status": status,
                        })

        print(f"Matched: {len(rows)}  Unmatched: {len(unmatched)}")
        return rows, unmatched

    def save_pairs(self, rows: List[Dict[str, object]]) -> None:
        with open(self.output_path, "w", encoding="utf-8") as f:
            json.dump(rows, f, indent=2)
        print(f"Saved matched pairs to: {self.output_path}")

    def save_unmatched(self, unmatched: List[Dict[str, object]]) -> None:
        with open(self.unmatched_path, "w", encoding="utf-8") as f:
            json.dump(unmatched, f, indent=2)
        print(f"Saved unmatched entries to: {self.unmatched_path}")


if __name__ == "__main__":
    builder = ConditionDrugPairBuilder(
        input_dir="data",
        output_dir="processed_data",
        mesh_path="mesh_data/desc2026.xml",
        fuzzy_cutoff=0.80,
        token_jaccard_min=0.60,
        include_placebo=False,
        keep_unmatched_debug_fields=True,
    )
    matched, unmatched = builder.extract_pairs()
    builder.save_pairs(matched)
    builder.save_unmatched(unmatched)
