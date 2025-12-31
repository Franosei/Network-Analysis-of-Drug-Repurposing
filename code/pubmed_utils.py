"""
pubmed_utils.py

Semantic prior construction from PubMed abstracts using an LLM classifier.

Key properties (aligned with your manuscript and the updated bayesian_predictor.py):
- Uses PubMed (db=pubmed), not PMC, to increase coverage.
- Uses stable article IDs in LLM prompting and outputs (no fuzzy title matching).
- Supports configurable filtering levels and dynamic date windows.
- Adds caching for PubMed search results to reduce repeated queries and improve reproducibility.
- Returns raw_counts (T, A, N) and total_articles for downstream Beta prior scaling.
- Integrates side-effect penalty via SideEffectUpdater, and captures gamma when available.

Fixes / Enhancements added in this version:
1) Disease synonym expansion layer (LLM):
   - Expands a disease term into a small set of high-quality synonyms / aliases (incl. MeSH-like headings if known).
   - Builds a single OR-clause for disease to improve recall (e.g., "breast cancer" OR "mammary carcinoma"...).
   - Prevents repetition by de-duplicating PMIDs after retrieval and by de-duping synonym strings.

2) Adds Introduction and Conclusion when available (without fabrication):
   - PubMed EFetch gives abstracts; Introduction/Conclusion are typically not present.
   - This module now *optionally* augments records using PMC full text when a PMCID exists for a PMID.
   - If no PMCID exists, intro/conclusion remain empty (no invented text).
   - Extraction is conservative: looks for titled sections like Introduction / Conclusion / Conclusions / Discussion.

Note:
- This keeps your overall structure and logic intact; changes are additive and guarded.
- SideEffectUpdater schema mismatch fix remains: reads p_final, not enhanced_prior.

"""

import os
import json
import time
import math
import re
import requests
from dataclasses import dataclass
from collections import Counter
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from xml.etree import ElementTree as ET
from openai import OpenAI
from dotenv import load_dotenv

from side_effect_updater import SideEffectUpdater


# ---------------------------------------------------------------------
# Environment and clients
# ---------------------------------------------------------------------
load_dotenv()

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "").strip()
if not OPENAI_API_KEY:
    raise RuntimeError("OPENAI_API_KEY is not set.")

NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()  # optional but strongly recommended
NCBI_EMAIL = os.getenv("NCBI_EMAIL", os.getenv("EMAIL", "oseifrancis633@gmail.com")).strip()

client = OpenAI(api_key=OPENAI_API_KEY)
updater = SideEffectUpdater(OPENAI_API_KEY)


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
@dataclass(frozen=True)
class PubMedSearchConfig:
    """
    Controls PubMed retrieval and filtering behavior.
    """
    # retrieval limits (hard ceiling for max_articles unless you change this config)
    max_results: int = 200
    # date window (dynamic by default)
    years_back: int = 10
    # rate limiting
    esearch_sleep_s: float = 0.34
    efetch_sleep_s: float = 0.34
    elink_sleep_s: float = 0.34
    pmc_fetch_sleep_s: float = 0.34
    # batching
    efetch_batch_size: int = 100  # PubMed efetch supports up to 200; keep safe
    # cache
    cache_dir: str = "cache"
    cache_file: str = "pubmed_esearch_cache.json"
    synonym_cache_file: str = "disease_synonym_cache.json"
    pmcid_cache_file: str = "pmid_to_pmcid_cache.json"
    # optional PMC augmentation
    enable_pmc_augmentation: bool = True
    # max pmc augmentation per run (keeps this safe & fast)
    max_pmc_augment: int = 50


@dataclass(frozen=True)
class LLMConfig:
    """
    Controls LLM classification behavior.
    """
    model: str = "gpt-4o-mini"
    delay_s: float = 2.0
    batch_size: int = 5
    max_retries: int = 2


@dataclass(frozen=True)
class SynonymConfig:
    """
    Controls disease synonym expansion behavior.
    """
    model: str = "gpt-4o-mini"
    max_synonyms: int = 6
    delay_s: float = 1.0
    max_retries: int = 2


# ---------------------------------------------------------------------
# Cache utilities
# ---------------------------------------------------------------------
class JsonCache:
    """
    Simple JSON cache.
    """
    def __init__(self, cache_dir: str, cache_file: str):
        self.cache_dir = cache_dir
        self.cache_path = os.path.join(cache_dir, cache_file)
        os.makedirs(cache_dir, exist_ok=True)
        self._data: Dict[str, Any] = self._load()

    def _load(self) -> Dict[str, Any]:
        if not os.path.exists(self.cache_path):
            return {}
        try:
            with open(self.cache_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}

    def get(self, key: str) -> Optional[Any]:
        return self._data.get(key)

    def set(self, key: str, value: Any) -> None:
        self._data[key] = value
        try:
            with open(self.cache_path, "w", encoding="utf-8") as f:
                json.dump(self._data, f, indent=2, ensure_ascii=False)
        except Exception:
            # cache failures should not crash the pipeline
            pass


# ---------------------------------------------------------------------
# Synonym expansion (LLM)
# ---------------------------------------------------------------------
class DiseaseSynonymExpander:
    """
    Expands disease names into synonyms / aliases to improve PubMed recall.

    Output is de-duplicated and includes the original disease term.
    Cached to avoid repeated LLM calls.
    """
    def __init__(self, cfg: SynonymConfig, cache: JsonCache):
        self.cfg = cfg
        self.cache = cache

    @staticmethod
    def _norm(s: str) -> str:
        return re.sub(r"\s+", " ", (s or "").strip().lower())

    @staticmethod
    def _dedup_keep_order(items: List[str]) -> List[str]:
        seen = set()
        out: List[str] = []
        for x in items:
            k = DiseaseSynonymExpander._norm(x)
            if not k or k in seen:
                continue
            seen.add(k)
            out.append(x.strip())
        return out

    def expand(self, disease: str) -> List[str]:
        disease = (disease or "").strip()
        if not disease:
            return []

        cache_key = f"disease_synonyms::{self._norm(disease)}::k={int(self.cfg.max_synonyms)}::m={self.cfg.model}"
        cached = self.cache.get(cache_key)
        if isinstance(cached, list) and cached:
            # ensure original is included
            items = self._dedup_keep_order([disease] + [str(x) for x in cached if str(x).strip()])
            return items[: max(1, int(self.cfg.max_synonyms) + 1)]

        prompt = f"""
You are assisting with PubMed query expansion.

Given the disease term below, return a JSON array of up to {int(self.cfg.max_synonyms)} alternative names/synonyms
that are commonly used in biomedical literature.

Rules:
- Return ONLY valid JSON (a list of strings).
- Do NOT include duplicates or trivial rephrasings (e.g., just adding/removing "disease").
- Prefer specific clinical/biomedical synonyms and well-known aliases.
- If the term is already a MeSH heading, you may include close aliases used in titles/abstracts.
- Avoid overly broad terms that would explode recall (e.g., for "breast cancer", do not return just "neoplasms").

Disease term: "{disease}"
""".strip()

        for attempt in range(self.cfg.max_retries + 1):
            try:
                resp = client.chat.completions.create(
                    model=self.cfg.model,
                    messages=[
                        {"role": "system", "content": "Return strictly valid JSON. No commentary."},
                        {"role": "user", "content": prompt},
                    ],
                    temperature=0,
                )
                content = (resp.choices[0].message.content or "").strip()
                content = content.replace("```json", "").replace("```", "").strip()
                parsed = json.loads(content)

                if not isinstance(parsed, list):
                    raise ValueError("Synonym output is not a JSON list.")

                syns = [str(x).strip() for x in parsed if isinstance(x, str) and str(x).strip()]
                syns = self._dedup_keep_order(syns)

                # Always include the original disease first
                all_terms = self._dedup_keep_order([disease] + syns)
                all_terms = all_terms[: max(1, int(self.cfg.max_synonyms) + 1)]

                # Cache without the original (so we can re-add consistently)
                self.cache.set(cache_key, syns[: int(self.cfg.max_synonyms)])
                time.sleep(self.cfg.delay_s)
                return all_terms

            except Exception as e:
                print(f"[WARN] Disease synonym expansion failed (attempt {attempt + 1}): {e}")
                time.sleep(1.5)

        # fallback: only original
        return [disease]


# ---------------------------------------------------------------------
# PubMed retrieval
# ---------------------------------------------------------------------
class PubMedClient:
    """
    Thin PubMed E-utilities client (ESearch + EFetch) with optional PMC augmentation.
    """
    def __init__(self, cfg: PubMedSearchConfig):
        self.cfg = cfg
        self.cache = JsonCache(cfg.cache_dir, cfg.cache_file)
        self.syn_cache = JsonCache(cfg.cache_dir, cfg.synonym_cache_file)
        self.pmcid_cache = JsonCache(cfg.cache_dir, cfg.pmcid_cache_file)

        self.syn_expander = DiseaseSynonymExpander(SynonymConfig(), self.syn_cache)

    def _dynamic_date_filter(self) -> str:
        end_year = datetime.now().year
        start_year = end_year - self.cfg.years_back
        return f"{start_year}:{end_year}[dp]"

    @staticmethod
    def _quote_term(t: str) -> str:
        t = (t or "").strip()
        if not t:
            return ""
        # Quote if it contains spaces or special characters
        if re.search(r"\s|[():\[\]\"']", t):
            return f"\"{t.replace('\"', '')}\""
        return t

    def build_query(self, drug: str, disease_terms: List[str], filter_level: str = "high") -> str:
        """
        Build a PubMed query with optional filters.

        disease_terms is an OR-clause across synonyms (Title/Abstract + MeSH Terms).
        """
        filter_level = (filter_level or "high").strip().lower()
        date_filter = self._dynamic_date_filter()

        # De-dup disease terms
        disease_terms = [t.strip() for t in (disease_terms or []) if str(t).strip()]
        disease_terms = DiseaseSynonymExpander._dedup_keep_order(disease_terms)
        if not disease_terms:
            # fall back
            disease_terms = []

        if filter_level == "minimal":
            if disease_terms:
                disease_clause = " OR ".join(
                    [
                        f"({self._quote_term(t)}[Title/Abstract] OR {self._quote_term(t)}[MeSH Terms])"
                        for t in disease_terms
                    ]
                )
                return f"({drug}) AND ({disease_clause})"
            return f"{drug}"

        # Drug clause: keep it in Title/Abstract for focus
        drug_clause = f"({self._quote_term(drug)}[Title/Abstract])" if drug else ""

        # Disease clause: include both Title/Abstract and MeSH Terms for each synonym
        if disease_terms:
            disease_clause = " OR ".join(
                [
                    f"({self._quote_term(t)}[Title/Abstract] OR {self._quote_term(t)}[MeSH Terms])"
                    for t in disease_terms
                ]
            )
            disease_clause = f"({disease_clause})"
        else:
            # if no terms, keep empty
            disease_clause = ""

        # Combine core term
        if drug_clause and disease_clause:
            base = f"({drug_clause} AND {disease_clause})"
        elif drug_clause:
            base = drug_clause
        elif disease_clause:
            base = disease_clause
        else:
            base = ""

        additions = [date_filter, "hasabstract[text]"]

        if filter_level in ("high", "strict"):
            additions.append("english[la]")

        if filter_level == "strict":
            additions.append("humans[mh]")
            additions.append("(Randomized Controlled Trial[pt] OR Meta-Analysis[pt] OR systematic[sb])")

        return f"{base} AND " + " AND ".join(additions) if base else " AND ".join(additions)

    def expand_disease_terms(self, disease: str) -> List[str]:
        return self.syn_expander.expand(disease)

    def esearch_pmids(self, query: str, max_results: int, use_cache: bool = True) -> List[str]:
        """
        Retrieve PMIDs for a query, sorted by relevance.
        Uses caching keyed on (query, max_results).
        """
        max_results = int(max(1, max_results))
        cache_key = f"{query}||max={max_results}"

        if use_cache:
            cached = self.cache.get(cache_key)
            if isinstance(cached, list) and cached:
                return cached[:max_results]

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": 0,
            "sort": "relevance",
            "email": NCBI_EMAIL,
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            r = requests.get(url, params=params, timeout=15)
            r.raise_for_status()
            data = r.json()
            total_count = int(data.get("esearchresult", {}).get("count", 0))
            if total_count <= 0:
                return []

            n = min(total_count, max_results)

            pmids: List[str] = []
            batch_size = 5000
            for start in range(0, n, batch_size):
                params["retstart"] = start
                params["retmax"] = min(batch_size, n - start)
                r = requests.get(url, params=params, timeout=15)
                r.raise_for_status()
                data = r.json()
                pmids.extend(data.get("esearchresult", {}).get("idlist", []))
                time.sleep(self.cfg.esearch_sleep_s)

            pmids = pmids[:n]
            if use_cache:
                self.cache.set(cache_key, pmids)
            return pmids

        except Exception as e:
            print(f"[ERROR] PubMed esearch failed for query='{query[:120]}': {e}")
            return []

    def efetch_abstracts(self, pmids: List[str]) -> List[Dict[str, Any]]:
        """
        Fetch PubMed records (title + abstract) for PMIDs.

        Returns list of dicts:
          {
            "id": int,
            "pmid": str,
            "title": str,
            "abstract": str,
            "abstract_sections": [str],
            "introduction": str,
            "conclusion": str,
            "pmcid": Optional[str]
          }
        """
        if not pmids:
            return []

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        out: List[Dict[str, Any]] = []
        total = len(pmids)

        for i in range(0, total, self.cfg.efetch_batch_size):
            batch = pmids[i : i + self.cfg.efetch_batch_size]
            params = {
                "db": "pubmed",
                "id": ",".join(batch),
                "retmode": "xml",
                "rettype": "abstract",
                "email": NCBI_EMAIL,
            }
            if NCBI_API_KEY:
                params["api_key"] = NCBI_API_KEY

            try:
                r = requests.get(url, params=params, timeout=20)
                r.raise_for_status()
                root = ET.fromstring(r.content)

                for pubmed_article in root.findall(".//PubmedArticle"):
                    pmid_elem = pubmed_article.find(".//PMID")
                    pmid = pmid_elem.text.strip() if pmid_elem is not None and pmid_elem.text else ""

                    title_elem = pubmed_article.find(".//ArticleTitle")
                    title = "".join(title_elem.itertext()).strip() if title_elem is not None else "No Title"

                    sections: List[str] = []
                    abstract_parts: List[str] = []

                    for abstract_elem in pubmed_article.findall(".//Abstract/AbstractText"):
                        label = (abstract_elem.get("Label") or "").strip()
                        txt = "".join(abstract_elem.itertext()).strip()
                        if not txt:
                            continue
                        if label:
                            sections.append(f"{label}: {txt}")
                            abstract_parts.append(f"{label}: {txt}")
                        else:
                            sections.append(txt)
                            abstract_parts.append(txt)

                    abstract = " ".join(abstract_parts).strip() if abstract_parts else ""
                    if not abstract:
                        continue

                    out.append(
                        {
                            "pmid": pmid,
                            "title": title,
                            "abstract": abstract,
                            "abstract_sections": sections,
                            "introduction": "",   # will be augmented from PMC when possible
                            "conclusion": "",     # will be augmented from PMC when possible
                            "pmcid": None,
                        }
                    )

                time.sleep(self.cfg.efetch_sleep_s)

            except Exception as e:
                print(f"[ERROR] PubMed efetch batch failed (batch {i//self.cfg.efetch_batch_size + 1}): {e}")
                time.sleep(self.cfg.efetch_sleep_s)

        # Stable integer IDs for LLM mapping
        for idx, rec in enumerate(out):
            rec["id"] = idx

        # Optional PMC augmentation (no fabrication)
        if self.cfg.enable_pmc_augmentation and out:
            try:
                self._augment_with_pmc(out)
            except Exception as e:
                print(f"[WARN] PMC augmentation failed: {e}")

        return out

    # -----------------------------
    # PMC augmentation helpers
    # -----------------------------
    def _pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """
        Resolve PMID -> PMCID via NCBI elink (cached).
        """
        pmid = (pmid or "").strip()
        if not pmid:
            return None

        cache_key = f"pmid2pmcid::{pmid}"
        cached = self.pmcid_cache.get(cache_key)
        if isinstance(cached, str) and cached.strip():
            return cached.strip()
        if cached is None:
            # proceed to lookup
            pass

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
        params = {
            "dbfrom": "pubmed",
            "db": "pmc",
            "id": pmid,
            "retmode": "xml",
            "email": NCBI_EMAIL,
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            r = requests.get(url, params=params, timeout=20)
            r.raise_for_status()
            root = ET.fromstring(r.content)

            # elink to PMC gives IDs like "PMC1234567" in some contexts or numeric IDs (PMCID is typically PMCxxxxxxx).
            # We'll try to parse LinkSetDb/Link/Id, then later confirm via efetch response if needed.
            link_ids = [x.text.strip() for x in root.findall(".//LinkSetDb/Link/Id") if x is not None and x.text]
            if not link_ids:
                self.pmcid_cache.set(cache_key, "")
                return None

            # Those are PMC internal numeric IDs. Convert numeric PMC id -> PMCID using esummary.
            pmc_numeric = link_ids[0]
            pmcid = self._pmc_numeric_to_pmcid(pmc_numeric)
            self.pmcid_cache.set(cache_key, pmcid or "")
            time.sleep(self.cfg.elink_sleep_s)
            return pmcid

        except Exception:
            self.pmcid_cache.set(cache_key, "")
            return None

    def _pmc_numeric_to_pmcid(self, pmc_numeric_id: str) -> Optional[str]:
        """
        Convert PMC numeric ID (from elink) to PMCID (PMCxxxxxxx) via esummary.
        """
        pmc_numeric_id = (pmc_numeric_id or "").strip()
        if not pmc_numeric_id:
            return None

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            "db": "pmc",
            "id": pmc_numeric_id,
            "retmode": "json",
            "email": NCBI_EMAIL,
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            r = requests.get(url, params=params, timeout=20)
            r.raise_for_status()
            data = r.json()
            result = data.get("result", {})
            rec = result.get(str(pmc_numeric_id), {})
            # "articleids" often includes {"idtype":"pmcid","value":"PMCxxxxxxx"}
            for aid in rec.get("articleids", []) or []:
                if str(aid.get("idtype", "")).lower() == "pmcid":
                    v = str(aid.get("value", "")).strip()
                    return v if v else None
            return None
        except Exception:
            return None

    @staticmethod
    def _norm_title(s: str) -> str:
        return re.sub(r"\s+", " ", (s or "").strip().lower())

    @staticmethod
    def _extract_pmc_sections(pmc_xml: bytes) -> Tuple[str, str]:
        """
        Extract Introduction and Conclusion (or closest equivalents) from PMC NXML.
        Conservative: returns empty strings if not found.
        """
        try:
            root = ET.fromstring(pmc_xml)
        except Exception:
            return "", ""

        def sec_title(sec_elem: ET.Element) -> str:
            t = sec_elem.find("./title")
            return "".join(t.itertext()).strip() if t is not None else ""

        def sec_text(sec_elem: ET.Element) -> str:
            return re.sub(r"\s+", " ", "".join(sec_elem.itertext()).strip())

        intro = ""
        concl = ""

        # Look through body sections
        for sec in root.findall(".//body//sec"):
            title = sec_title(sec).lower()

            if not intro and title in {"introduction", "background"}:
                txt = sec_text(sec)
                if txt:
                    intro = txt

            # conclusion could be Conclusion(s) or Discussion/Concluding remarks
            if not concl and title in {"conclusion", "conclusions", "concluding remarks"}:
                txt = sec_text(sec)
                if txt:
                    concl = txt

        # Fallback: if still no conclusion, try discussion as "closest" (but only if explicitly titled discussion)
        if not concl:
            for sec in root.findall(".//body//sec"):
                title = sec_title(sec).lower()
                if title == "discussion":
                    txt = sec_text(sec)
                    if txt:
                        concl = txt
                        break

        return intro.strip(), concl.strip()

    def _fetch_pmc_fulltext_xml(self, pmcid: str) -> Optional[bytes]:
        """
        Fetch PMC fulltext (NXML) by PMCID.
        """
        pmcid = (pmcid or "").strip()
        if not pmcid:
            return None

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "pmc",
            "id": pmcid,
            "retmode": "xml",
            "rettype": "full",
            "email": NCBI_EMAIL,
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            r = requests.get(url, params=params, timeout=30)
            r.raise_for_status()
            time.sleep(self.cfg.pmc_fetch_sleep_s)
            return r.content
        except Exception:
            return None

    def _augment_with_pmc(self, records: List[Dict[str, Any]]) -> None:
        """
        For records with an available PMCID, fetch PMC full text and extract intro/conclusion.
        Limits the number of augmentations per run to cfg.max_pmc_augment.
        """
        augmented = 0
        for rec in records:
            if augmented >= int(self.cfg.max_pmc_augment):
                break

            pmid = rec.get("pmid", "")
            if not pmid:
                continue

            pmcid = self._pmid_to_pmcid(pmid)
            if not pmcid:
                continue

            pmc_xml = self._fetch_pmc_fulltext_xml(pmcid)
            if not pmc_xml:
                continue

            intro, concl = self._extract_pmc_sections(pmc_xml)
            # Only set if found (no fabrication)
            if intro:
                rec["introduction"] = intro
            if concl:
                rec["conclusion"] = concl
            rec["pmcid"] = pmcid

            if intro or concl:
                augmented += 1


# ---------------------------------------------------------------------
# LLM classification
# ---------------------------------------------------------------------
class LLMClassifier:
    """
    Builds a semantic prior for a drug-disease pair by:
      1) retrieving PubMed abstracts (augmented with intro/conclusion if PMC exists)
      2) classifying each record as therapeutic/adverse/irrelevant
      3) computing raw and penalised priors from counts
      4) applying SideEffectUpdater to produce enhanced_prior (= p_final)
    """

    VALID_CATEGORIES = {"therapeutic", "adverse", "irrelevant"}

    def __init__(
        self,
        search_cfg: Optional[PubMedSearchConfig] = None,
        llm_cfg: Optional[LLMConfig] = None,
    ):
        self.search_cfg = search_cfg or PubMedSearchConfig()
        self.llm_cfg = llm_cfg or LLMConfig()
        self.pubmed = PubMedClient(self.search_cfg)

    @staticmethod
    def _clamp01(x: Any, default: float = 0.0) -> float:
        try:
            v = float(x)
        except Exception:
            v = float(default)
        return max(min(v, 1.0), 0.0)

    def _classify_batch(self, batch: List[Dict[str, Any]], drug: str, disease: str) -> List[Dict[str, Any]]:
        """
        Classify a batch of records. Output must be a JSON list of:
          { "id": <int>, "category": "therapeutic"|"adverse"|"irrelevant" }
        """
        payload = [
            {
                "id": rec["id"],
                "pmid": rec.get("pmid", ""),
                "pmcid": rec.get("pmcid", None),
                "title": rec.get("title", "No Title"),
                "abstract": rec.get("abstract", ""),
                "abstract_sections": rec.get("abstract_sections", []),
                "introduction": rec.get("introduction", ""),
                "conclusion": rec.get("conclusion", ""),
            }
            for rec in batch
        ]

        prompt = f"""
You are a biomedical assistant. Classify each PubMed record for the relationship between the specified drug and disease.

Drug: {drug}
Disease: {disease}

Categories (mutually exclusive):
- therapeutic: evidence suggests the drug benefits or treats the disease
- adverse: evidence suggests the drug induces, worsens, or is associated with harm relevant to the disease
- irrelevant: no meaningful therapeutic or adverse relationship is supported by the provided text

You may use:
- abstract (always present)
- introduction / conclusion (only present when available from PMC full text; may be empty)

Return ONLY valid JSON as a list. Each item must include:
- "id": integer (must match an input id)
- "category": one of ["therapeutic", "adverse", "irrelevant"]

Do not include any additional keys. Do not include any commentary.

Input records:
{json.dumps(payload, ensure_ascii=False)}
""".strip()

        for attempt in range(self.llm_cfg.max_retries + 1):
            try:
                resp = client.chat.completions.create(
                    model=self.llm_cfg.model,
                    messages=[
                        {
                            "role": "system",
                            "content": (
                                "You are a biomedical assistant categorizing drug-disease relationships from PubMed/PMC text. "
                                "Return strictly valid JSON. No explanations."
                            ),
                        },
                        {"role": "user", "content": prompt},
                    ],
                    temperature=0,
                )

                content = (resp.choices[0].message.content or "").strip()
                content = content.replace("```json", "").replace("```", "").strip()
                parsed = json.loads(content)

                if not isinstance(parsed, list):
                    raise ValueError("LLM output is not a JSON list.")

                out: List[Dict[str, Any]] = []
                seen_ids = set()

                for item in parsed:
                    if not isinstance(item, dict):
                        continue
                    if "id" not in item or "category" not in item:
                        continue
                    try:
                        rid = int(item["id"])
                    except Exception:
                        continue
                    cat = str(item["category"]).strip().lower()
                    if cat not in self.VALID_CATEGORIES:
                        continue
                    if rid in seen_ids:
                        continue
                    seen_ids.add(rid)
                    out.append({"id": rid, "category": cat})

                expected_ids = {int(r["id"]) for r in batch}
                got_ids = {int(x["id"]) for x in out}
                if expected_ids == got_ids:
                    return out

                raise ValueError(f"LLM batch mismatch: expected {len(expected_ids)} ids, got {len(got_ids)} ids.")

            except Exception as e:
                print(f"[WARN] LLM classification failed (attempt {attempt + 1}): {e}")
                time.sleep(2)

        return []

    def classify_abstracts(self, records: List[Dict[str, Any]], drug: str, disease: str) -> List[Dict[str, Any]]:
        """
        Classify all retrieved records in batches.
        Returns list of {id, category}.
        """
        total = len(records)
        if total == 0:
            return []

        bs = int(max(1, self.llm_cfg.batch_size))
        total_batches = math.ceil(total / bs)

        all_labels: List[Dict[str, Any]] = []

        for b in range(total_batches):
            batch = records[b * bs : (b + 1) * bs]
            labels = self._classify_batch(batch, drug, disease)

            if labels and len(labels) == len(batch):
                all_labels.extend(labels)
                print(f"[INFO] Classified batch {b + 1}/{total_batches} ({len(all_labels)}/{total})")
            else:
                print(f"[WARN] Batch {b + 1}/{total_batches} classification failed or incomplete.")

            time.sleep(self.llm_cfg.delay_s)

        dedup = {int(x["id"]): x["category"] for x in all_labels}
        return [{"id": k, "category": v} for k, v in dedup.items()]

    def _neutral_result(self) -> Dict[str, Any]:
        """
        Consistent neutral fallback when no evidence is available.
        This avoids contradictory priors (e.g., prior=0.0 but enhanced_prior=0.5).
        """
        return {
            "prior": 0.5,
            "penalised_prior": 0.5,
            "enhanced_prior": 0.5,
            "gamma": None,
            "raw_counts": Counter(),
            "labelled_abstracts": [],
            "total_articles": 0,
        }

    @staticmethod
    def _dedup_pmids_keep_order(pmids: List[str]) -> List[str]:
        seen = set()
        out: List[str] = []
        for p in pmids:
            p = (p or "").strip()
            if not p or p in seen:
                continue
            seen.add(p)
            out.append(p)
        return out

    def build_semantic_prior(
        self,
        drug: str,
        disease: str,
        max_articles: Optional[int] = None,
        filter_level: str = "high",
        save_dir: str = "literatures",
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """
        Main entry point: builds a semantic prior dictionary.

        Returns:
          {
            "prior": raw_prior,
            "penalised_prior": penalised_prior,
            "enhanced_prior": enhanced_prior (= p_final),
            "gamma": gamma_or_None,
            "raw_counts": Counter({"therapeutic": T, "adverse": A, "irrelevant": N}),
            "labelled_abstracts": [{"id","pmid","Title","category"}],
            "total_articles": M
          }
        """
        drug = (drug or "").strip()
        disease = (disease or "").strip()

        if not drug or not disease:
            return self._neutral_result()

        # 1) Expand disease terms to improve retrieval
        disease_terms = self.pubmed.expand_disease_terms(disease)
        if disease_terms:
            print(f"[INFO] Disease expansion terms: {disease_terms}")

        # 2) Build query using OR across disease terms
        query = self.pubmed.build_query(drug=drug, disease_terms=disease_terms, filter_level=filter_level)

        limit = int(max_articles) if max_articles is not None else int(self.search_cfg.max_results)
        limit = max(1, limit)

        # Respect configured hard ceiling (explicit)
        if limit > int(self.search_cfg.max_results):
            print(f"[INFO] max_articles={limit} capped to search_cfg.max_results={self.search_cfg.max_results}")
            limit = int(self.search_cfg.max_results)

        print(f"[INFO] PubMed query (level={filter_level}): {query}")

        pmids = self.pubmed.esearch_pmids(query=query, max_results=limit, use_cache=use_cache)
        pmids = self._dedup_pmids_keep_order(pmids)
        if not pmids:
            return self._neutral_result()

        records = self.pubmed.efetch_abstracts(pmids)
        if not records:
            return self._neutral_result()

        labels = self.classify_abstracts(records, drug, disease)
        label_map = {int(x["id"]): x["category"] for x in labels}

        labelled_records: List[Dict[str, Any]] = []
        for rec in records:
            rid = int(rec["id"])
            if rid not in label_map:
                continue
            cat = label_map[rid]
            labelled_records.append(
                {
                    "id": rid,
                    "pmid": rec.get("pmid", ""),
                    "pmcid": rec.get("pmcid", None),
                    "Title": rec.get("title", "No Title"),
                    "category": cat,
                    "abstract": rec.get("abstract", ""),
                    "introduction": rec.get("introduction", ""),
                    "conclusion": rec.get("conclusion", ""),
                }
            )

        if not labelled_records:
            return self._neutral_result()

        # Save for auditability
        os.makedirs(save_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_drug = drug.replace("/", "_")
        safe_disease = disease.replace("/", "_")
        outpath = os.path.join(save_dir, f"classified_pubmed_{safe_drug}_{safe_disease}_{timestamp}.json")
        try:
            with open(outpath, "w", encoding="utf-8") as f:
                json.dump(labelled_records, f, indent=2, ensure_ascii=False)
            print(f"[INFO] Saved classified abstracts to: {outpath}")
        except Exception:
            pass

        # Compute counts and priors
        counts = Counter(x["category"] for x in labelled_records)
        total = int(sum(counts.values()))
        T = int(counts.get("therapeutic", 0))
        A = int(counts.get("adverse", 0))
        N = int(counts.get("irrelevant", 0))

        raw_prior = (T / total) if total else 0.5
        penalised_prior = max((T - 2 * A) / total, 0.0) if total else 0.5

        raw_prior = self._clamp01(raw_prior, default=0.5)
        penalised_prior = self._clamp01(penalised_prior, default=0.5)

        # Side-effect adjustment (FIX: read p_final from SideEffectUpdater dict)
        gamma: Optional[float] = None
        enhanced_prior = penalised_prior

        try:
            upd = updater.update_prior(drug, disease, penalised_prior)

            # Supported return styles:
            # - dict {"p_final": float, "gamma": float, ...}  [current SideEffectUpdater]
            # - dict {"enhanced_prior": float, "gamma": float, ...} [legacy]
            # - tuple (enhanced_prior, gamma)
            # - float enhanced_prior
            if isinstance(upd, dict):
                enhanced_prior = float(upd.get("p_final", upd.get("enhanced_prior", penalised_prior)))
                g = upd.get("gamma", None)
                gamma = float(g) if g is not None else None
            elif isinstance(upd, tuple) and len(upd) >= 2:
                enhanced_prior = float(upd[0])
                gamma = float(upd[1]) if upd[1] is not None else None
            else:
                enhanced_prior = float(upd)

        except Exception as e:
            print(f"[WARN] Side effect adjustment failed: {e}")
            enhanced_prior = penalised_prior
            gamma = None

        enhanced_prior = self._clamp01(enhanced_prior, default=penalised_prior)
        if gamma is not None:
            gamma = self._clamp01(gamma, default=0.0)

        return {
            "prior": raw_prior,
            "penalised_prior": penalised_prior,
            "enhanced_prior": enhanced_prior,  # corresponds to p_final
            "gamma": gamma,
            "raw_counts": counts,
            "labelled_abstracts": [
                {
                    "Title": x["Title"],
                    "category": x["category"],
                    "pmid": x["pmid"],
                    "pmcid": x.get("pmcid", None),
                    "id": x["id"],
                }
                for x in labelled_records
            ],
            "total_articles": total,
        }


# ---------------------------------------------------------------------
# Script usage
# ---------------------------------------------------------------------
# if __name__ == "__main__":
#     classifier = LLMClassifier(
#         search_cfg=PubMedSearchConfig(max_results=200, years_back=10, enable_pmc_augmentation=True),
#         llm_cfg=LLMConfig(model="gpt-4o-mini", delay_s=2.0, batch_size=5, max_retries=2),
#     )
#
#     drug = "metformin"
#     disease = "chronic kidney disease"
#
#     result = classifier.build_semantic_prior(
#         drug=drug,
#         disease=disease,
#         max_articles=150,
#         filter_level="high",
#         use_cache=True,
#     )
#
#     print("\n" + "=" * 60)
#     print("FINAL SUMMARY")
#     print("=" * 60)
#     print(f"Drug: {drug}")
#     print(f"Disease: {disease}")
#     print(f"Total Articles Processed: {result.get('total_articles', 0)}")
#     print(f"Counts: {dict(result.get('raw_counts', {}))}")
#     print(f"Raw prior: {result.get('prior')}")
#     print(f"Penalised prior: {result.get('penalised_prior')}")
#     print(f"Enhanced prior (p_final): {result.get('enhanced_prior')}")
#     print(f"Gamma: {result.get('gamma')}")
#     print("=" * 60)
