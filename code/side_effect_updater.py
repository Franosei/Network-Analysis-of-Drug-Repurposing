import json
from typing import Any, Dict, Set

import requests
from openai import OpenAI


class SideEffectUpdater:
    """
    Retrieves post-marketing adverse events from openFDA and applies a bounded penalty
    to a literature-derived prior based on semantic overlap between side effects and
    the target disease.

    Contract (stable return schema):
      {
        "p_final": float,            # adjusted prior after penalty (clamped to [0,1])
        "gamma": float,              # overlap score in [0,1]
        "penalty_scale": float,      # multiplier applied to gamma (default 0.5)
        "matching_effects": list,    # side effects judged relevant to disease
        "relation": bool             # whether an overlap relation exists
      }
    """

    def __init__(self, api_key: str, fda_limit: int = 50, timeout: int = 10):
        self.client = OpenAI(api_key=api_key)
        self.fda_url = "https://api.fda.gov/drug/event.json"
        self.fda_limit = int(max(1, fda_limit))
        self.timeout = int(max(1, timeout))

    @staticmethod
    def _clamp01(x: Any) -> float:
        try:
            v = float(x)
        except Exception:
            return 0.0
        return max(min(v, 1.0), 0.0)

    def get_side_effects(self, drug: str) -> Set[str]:
        """
        Get a set of unique adverse event terms from openFDA.

        Notes:
          - openFDA may return errors for rare terms or if no results exist.
          - We intentionally keep a modest limit for speed/cost control.
        """
        drug = (drug or "").strip()
        if not drug:
            return set()

        try:
            resp = requests.get(
                self.fda_url,
                params={
                    "search": f'patient.drug.medicinalproduct:"{drug.lower()}"',
                    "limit": self.fda_limit,
                },
                timeout=self.timeout,
            )
            resp.raise_for_status()
            data = resp.json()

            side_effects = {
                reaction["reactionmeddrapt"].lower()
                for case in data.get("results", [])
                for reaction in case.get("patient", {}).get("reaction", [])
                if reaction.get("reactionmeddrapt")
            }
            return side_effects

        except Exception:
            return set()

    def update_prior(self, drug: str, disease: str, current_prior: float, penalty_scale: float = 0.5) -> Dict[str, Any]:
        """
        Apply a side-effect overlap penalty to the current prior.

        Implements:
          P_final = P_current * (1 - penalty_scale * gamma) if relation=True
          otherwise P_final = P_current

        Where gamma in [0,1] is obtained from an LLM assessment.

        Returns the stable dict schema documented in the class docstring.
        """
        p0 = self._clamp01(current_prior)
        penalty_scale = max(float(penalty_scale), 0.0)

        side_effects = self.get_side_effects(drug)
        if not side_effects:
            return {
                "p_final": p0,
                "gamma": 0.0,
                "penalty_scale": penalty_scale,
                "matching_effects": [],
                "relation": False,
            }

        try:
            response = self.client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {
                        "role": "system",
                        "content": (
                            "Analyze whether the drug's side effects are clinically or semantically "
                            "related to the target disease. Return JSON with keys: "
                            "relation (bool), confidence (0-1), matching_effects (list of strings)."
                        ),
                    },
                    {
                        "role": "user",
                        "content": (
                            f"Drug: {drug}\n"
                            f"Disease: {disease}\n"
                            f"Side effects (MedDRA PT terms): {', '.join(sorted(side_effects))}"
                        ),
                    },
                ],
                temperature=0.3,
                response_format={"type": "json_object"},
            )

            raw = (response.choices[0].message.content or "").strip()
            analysis = json.loads(raw) if raw else {}

            relation = bool(analysis.get("relation", False))
            gamma = self._clamp01(analysis.get("confidence", 0.0))
            matching_effects = analysis.get("matching_effects", [])
            if not isinstance(matching_effects, list):
                matching_effects = []

            if relation:
                penalty = self._clamp01(gamma * penalty_scale)
                p_final = self._clamp01(p0 * (1.0 - penalty))
            else:
                p_final = p0

            return {
                "p_final": p_final,
                "gamma": gamma,
                "penalty_scale": penalty_scale,
                "matching_effects": matching_effects,
                "relation": relation,
            }

        except Exception as e:
            print(f"[SideEffectUpdater Error] {str(e)}")
            return {
                "p_final": p0,
                "gamma": 0.0,
                "penalty_scale": penalty_scale,
                "matching_effects": [],
                "relation": False,
            }
