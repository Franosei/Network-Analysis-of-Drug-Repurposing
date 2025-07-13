import os
import json
import time
import requests
from collections import Counter
from datetime import datetime
from xml.etree import ElementTree as ET
from openai import OpenAI
from dotenv import load_dotenv
from difflib import get_close_matches
from side_effect_updater import SideEffectUpdater
       

# Load environment variables
load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
client = OpenAI(api_key=OPENAI_API_KEY)
updater = SideEffectUpdater(os.getenv("OPENAI_API_KEY"))


class LLMClassifier:
    def __init__(self, delay=3, model="gpt-4o-mini"):
        self.model = model
        self.delay = delay

    def clean_title(self, title):
        return " ".join(title.strip().lower().split())

    def match_title(self, target, candidate_titles):
        matches = get_close_matches(self.clean_title(target), [self.clean_title(c) for c in candidate_titles], n=1, cutoff=0.6)
        if matches:
            for original in candidate_titles:
                if self.clean_title(original) == matches[0]:
                    return original
        return None

    def pmc_get_pmids(self, term, max_count=30):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pmc",
            "term": term,
            "retmode": "json",
            "retmax": max_count,
            "email": "oseifrancis633@gmail.com"
        }
        try:
            r = requests.get(url, params=params, timeout=10)
            r.raise_for_status()
            data = r.json()
            return data.get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            print(f"[ERROR] Fetching PMIDs for '{term}':", e)
            return []

    def pmc_fetch_abstracts(self, pmids):
        articles = []
        for i in range(0, len(pmids), 20):
            batch = pmids[i:i+20]
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                "db": "pmc",
                "id": ",".join(batch),
                "retmode": "xml",
                "email": "oseifrancis633@gmail.com"
            }
            try:
                r = requests.get(url, params=params, timeout=10)
                r.raise_for_status()
                root = ET.fromstring(r.content)

                for article in root.findall(".//article"):
                    title_elem = article.find(".//article-title")
                    abstract_elem = article.find(".//abstract")
                    intro_text, conclusion_text = "", ""

                    for sec in article.findall(".//sec"):
                        sec_title = sec.find("title")
                        if sec_title is not None and sec_title.text:
                            sec_lower = sec_title.text.lower()
                            if "introduction" in sec_lower:
                                paras = sec.findall(".//p")
                                intro_text = " ".join("".join(p.itertext()).strip() for p in paras)
                            elif "conclusion" in sec_lower:
                                paras = sec.findall(".//p")
                                conclusion_text = " ".join("".join(p.itertext()).strip() for p in paras)

                    articles.append({
                        "title": "".join(title_elem.itertext()).strip() if title_elem is not None else "No Title",
                        "abstract": "".join(abstract_elem.itertext()).strip() if abstract_elem is not None else "No Abstract",
                        "introduction": intro_text or "No Introduction Found",
                        "conclusion": conclusion_text or "No Conclusion Found"
                    })

            except Exception as e:
                print("[ERROR] Fetching abstracts:", e)
            time.sleep(1)
        return articles

    def classify_batch(self, batch, drug, disease, max_retries=2):
        formatted = "\n\n".join([
            f"Title: {a['title']}\nIntroduction: {a['introduction']}\nAbstract: {a['abstract']}\nConclusion: {a['conclusion']}"
            for a in batch
        ])

        prompt = f"""
You are a biomedical AI assistant. Based on the provided articles, classify each article into one of the following categories:
- therapeutic: the drug helps treat the disease
- adverse: the drug causes or worsens the disease
- irrelevant: no meaningful relation between drug and disease

The drug is: **{drug}**  
The disease is: **{disease}**

Return ONLY a valid JSON list. Each item must include:
- "Title": string (match the article title)
- "category": one of ["therapeutic", "adverse", "irrelevant"]

Example output:
[
  {{"Title": "Example A", "category": "therapeutic"}},
  {{"Title": "Example B", "category": "irrelevant"}}
]

NO commentary. No explanation. No markdown. Just the JSON list.
Articles:
{formatted}
"""

        for attempt in range(max_retries + 1):
            try:
                response = client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "You are a biomedical assistant categorizing drug–disease relationships. This data will be used to build a semantic prior so we required more and accurte classification."},
                        {"role": "user", "content": prompt.strip()}
                    ],
                    temperature=0
                )
                content = response.choices[0].message.content.strip()
                content = content.replace("```json", "").replace("```", "").strip()
                result = json.loads(content)

                if isinstance(result, list) and all("category" in x for x in result):
                    return result

            except Exception as e:
                print(f"[WARN] LLM classification failed on attempt {attempt + 1}:", e)
                time.sleep(2)

        return []

    def classify_abstracts(self, articles, drug, disease, batch_size=3):
        all_labels = []
        for i in range(0, len(articles), batch_size):
            batch = articles[i:i+batch_size]
            labels = self.classify_batch(batch, drug, disease)

            if not labels or len(labels) != len(batch):
                print(f"[!] Batch {i//batch_size+1}: mismatch or empty result (expected {len(batch)}, got {len(labels)})")
                print(f"[!] Titles in batch: {[a['title'][:40] for a in batch]}")
            else:
                all_labels.extend(labels)

            time.sleep(self.delay)
        return all_labels

    def build_semantic_prior(self, drug, disease, max_count=30):
        print(f"Processing: {drug} → {disease}")
        term = f"{drug} AND {disease}"
        pmids = self.pmc_get_pmids(term, max_count)
        if not pmids:
            print("[!] No PMIDs found.")
            return None

        articles = self.pmc_fetch_abstracts(pmids)
        if not articles:
            print("[!] No articles found.")
            return None

        print(f" Retrieved {len(articles)} articles. Classifying...")
        labelled = self.classify_abstracts(articles, drug, disease)
        labelled_map = {item["Title"]: item["category"] for item in labelled}

        all_labelled = []
        all_articles = []
        for article in articles:
            matched_title = self.match_title(article["title"], labelled_map.keys())
            if matched_title:
                category = labelled_map[matched_title]
                article["category"] = category
                all_articles.append(article)
                all_labelled.append({"Title": article["title"], "category": category})

        if not all_labelled:
            print("[!] No successfully classified articles.")
            return None

        os.makedirs("literatures", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outpath = f"literatures/classified_articles_{drug}_{disease}_{timestamp}.json"
        with open(outpath, "w", encoding="utf-8") as f:
            json.dump(all_articles, f, indent=2)

        counts = Counter(x["category"] for x in all_labelled)
        total = sum(counts.values())
        therapeutic = counts.get("therapeutic", 0)
        adverse = counts.get("adverse", 0)

        raw_prior = therapeutic / total if total else 0
        penalised_prior = max((therapeutic - 2 * adverse) / total, 0) if total else 0

        # Call SideEffectUpdater here
        enhanced_prior = updater.update_prior(drug, disease, penalised_prior)

        print(f"Completed: {therapeutic} therapeutic, {adverse} adverse, {counts.get('irrelevant', 0)} irrelevant")

        return {
            "prior": raw_prior,
            "penalised_prior": penalised_prior,
            "enhanced_prior": enhanced_prior,
            "raw_counts": counts,
            "labelled_abstracts": all_labelled
        }