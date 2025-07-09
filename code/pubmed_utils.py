import os
import json
import numpy as np
from collections import Counter
from openai import OpenAI
from dotenv import load_dotenv
from time import sleep
import requests
from xml.etree import ElementTree as ET
from datetime import datetime


load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
client = OpenAI(api_key=OPENAI_API_KEY)

def expand_topics(drug, disease):
    prompt = f"""
Given the drug "{drug}" and the disease "{disease}", generate 5 different PubMed search topic phrasings 
that may yield the same useful literature results even if exact keywords are not matched.
every search should have both the drug and disease in the title or abstract OR RELATED.

Return them as a JSON list of strings.
"""
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a biomedical researcher skilled in PubMed querying."},
            {"role": "user", "content": prompt.strip()},
        ],
        temperature=0
    )

    content = response.choices[0].message.content.strip()
    try:
        if content.startswith("```json"):
            content = content.replace("```json", "").replace("```", "").strip()
        elif content.startswith("```"):
            content = content.replace("```", "").strip()
        return json.loads(content)
    except Exception as e:
        print("Failed to parse topic expansion:", e)
        print("Raw:", content)
        return [f"{drug} and {disease}"]

def pmc_get_pmids(term, max_count=5):
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
        print(f"Error fetching PMIDs for term '{term}':", e)
        return []

def pmc_fetch_abstracts(pmids):
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
                title = "".join(title_elem.itertext()).strip() if title_elem is not None else "No Title"
                abstract = "".join(abstract_elem.itertext()).strip() if abstract_elem is not None else "No Abstract"

                conclusion_text = ""
                for sec in article.findall(".//sec"):
                    sec_title = sec.find("title")
                    if sec_title is not None and "conclusion" in sec_title.text.lower():
                        paras = sec.findall(".//p")
                        conclusion_text = " ".join("".join(p.itertext()).strip() for p in paras)
                        break

                articles.append({
                    "title": title,
                    "abstract": abstract,
                    "conclusion": conclusion_text or "No Conclusion Found"
                })

        except Exception as e:
            print("Error fetching abstracts:", e)
        sleep(0.5)
    return articles

def classify_abstracts(articles,drug, disease):
    formatted = "\n\n".join([
        f"Title: {a['title']}\nAbstract: {a['abstract']}\nConclusion: {a['conclusion']}" 
        for a in articles
    ])
    prompt = f"""
Analyze the following PubMed articles if. For each  drug: {drug} and  diseases:{disease} pair, determine whether the drug {drug} has a therapeutic, adverse, or irrelevant effect on the disease {disease}. Use both the abstract and the conclusion to decide.
If the abstract discusses the drug as a cause of the disease (e.g., side effect), label as 'adverse' even if the word ‘treat’ appears.
Please therapeutic effect means the drug has a positive use for the condition.
irrelevant means it does not mention the drug and disease at all, or the relationship is not clear.
Be strict in the classification, they are going to be use for drug repurposing application.
Return a JSON list. Each item must be:
{{"Title": "...", "category": "therapeutic" | "adverse" | "irrelevant"}}

ONLY return the JSON list. No explanation or extra commentary.

Articles:
{formatted}
"""
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are an expert biomedical assistant that categorizes drug-disease relationships using abstract and conclusion."},
            {"role": "user", "content": prompt.strip()},
        ],
        temperature=0,
    )

    raw_content = response.choices[0].message.content.strip()
    if raw_content.startswith("```json"):
        raw_content = raw_content.replace("```json", "").replace("```", "").strip()
    elif raw_content.startswith("```"):
        raw_content = raw_content.replace("```", "").strip()

    try:
        return json.loads(raw_content)
    except Exception as e:
        print("Failed to parse GPT response:", e)
        print("Raw content:", response.choices[0].message.content)
        return []

def build_semantic_prior(drug, disease, max_count=5):
    print(f"Expanding topics for: {drug} and {disease}")
    topics = expand_topics(drug, disease)
    all_labelled = []
    all_articles = []

    for term in topics:
        print(f"Searching: {term}")
        pmids = pmc_get_pmids(term, max_count=max_count)
        if not pmids:
            continue
        articles = pmc_fetch_abstracts(pmids)
        if not articles:
            continue
        print(f"{len(articles)} articles retrieved.")
        labelled = classify_abstracts(articles,drug, disease)
        if labelled:
            all_labelled.extend(labelled)
            for i in range(len(labelled)):
                articles[i]["category"] = labelled[i]["category"]
            all_articles.extend(articles)

    if not all_labelled:
        print("No articles successfully classified.")
        return None

    # Save full results with titles, abstracts, conclusions, and categories
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outpath = f"literatures/classified_articles_{drug}_{disease}_{timestamp}.json"
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(all_articles, f, indent=2)

    # Compute prior
    counts = Counter([x["category"] for x in all_labelled])
    total = sum(counts.values())
    if total == 0:
        return {"prior": 0.0, "penalised_prior": 0.0, "raw_counts": counts}

    therapeutic = counts.get("therapeutic", 0)
    adverse = counts.get("adverse", 0)

    raw_prior = therapeutic / total
    penalised_prior = max((therapeutic - 2 * adverse) / total, 0)

    return {
        "prior": raw_prior,
        "penalised_prior": penalised_prior,
        "raw_counts": counts,
        "labelled_abstracts": all_labelled
    }
