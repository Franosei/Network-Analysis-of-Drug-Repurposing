import requests
import json
from openai import OpenAI

class SideEffectUpdater:
    def __init__(self, api_key):
        self.client = OpenAI(api_key=api_key)
        self.fda_url = "https://api.fda.gov/drug/event.json"
    
    def get_side_effects(self, drug):
        """Get unique side effects from FDA API"""
        try:
            response = requests.get(
                self.fda_url,
                params={'search': f'patient.drug.medicinalproduct:"{drug.lower()}"', 'limit': 50},
                timeout=10
            )
            data = response.json()
            return {
                reaction['reactionmeddrapt'].lower()
                for case in data.get('results', [])
                for reaction in case.get('patient', {}).get('reaction', [])
                if reaction.get('reactionmeddrapt')
            }
        except Exception:
            return set()
    
    def update_prior(self, drug, disease, current_prior):
        """Update prior based on side effect analysis"""
        side_effects = self.get_side_effects(drug)
        if not side_effects:
            return current_prior
        
        try:
            response = self.client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{
                    "role": "system",
                    "content": "Analyze if drug side effects relate to disease. Return JSON with: relation (bool), confidence (0-1), matching_effects (list)"
                }, {
                    "role": "user",
                    "content": f"Analyze if {drug}'s side effects relate to {disease}. Side effects: {', '.join(sorted(side_effects))}"
                }],
                temperature=0.3,
                response_format={"type": "json_object"}
            )
            
            analysis = json.loads(response.choices[0].message.content)

            if analysis.get('relation', False):
                penalty = analysis.get('confidence', 0) * 0.5
                new_prior = max(current_prior * (1 - penalty), 0)

                print(f"Matching side effects between '{drug}' and '{disease}': {analysis.get('matching_effects', [])}")
                print(f"LLM found a relation (confidence: {analysis.get('confidence', 0)}). Prior reduced by {round(penalty * 100, 1)}%.")

                return new_prior
            
        except Exception as e:
            print(f"[SideEffectUpdater Error] {str(e)}")
        
        return current_prior
