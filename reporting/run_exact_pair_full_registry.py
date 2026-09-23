"""Rebuild an exact drug--condition trial graph and audit one target pair.

This runner deliberately does not reuse the manuscript's sampled/fuzzy graph.
It downloads a current ClinicalTrials.gov registry snapshot, retains only
interventional studies with interventions explicitly typed as DRUG, and makes
study-level edges from punctuation-normalised *exact* field values.  The target
edge is removed before structural features are calculated.

It also retrieves every PubMed identifier returned by a declared title/abstract
query, downloads the records, and verifies that both target phrases occur in the
title/abstract after the same conservative normalisation.  It does not infer
therapeutic direction or safety from publication volume.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sqlite3
import sys
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import requests


CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CTGOV_FIELDS = ",".join(
    (
        "NCTId",
        "BriefTitle",
        "OverallStatus",
        "StartDate",
        "CompletionDate",
        "StudyFirstPostDate",
        "LastUpdatePostDate",
        "Condition",
        "InterventionName",
        "InterventionType",
        "Phase",
        "StudyType",
    )
)
USER_AGENT = "exact-pair-registry-audit/1.0 (research reproducibility script)"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalise_exact(value: Any) -> str:
    """Case/punctuation normalisation only; no stemming, ontology, or synonyms."""
    text = unicodedata.normalize("NFKC", str(value or "")).casefold()
    return re.sub(r"[^a-z0-9]+", " ", text).strip()


def contains_exact_phrase(text: str, phrase: str) -> bool:
    normalised = f" {normalise_exact(text)} "
    needle = f" {normalise_exact(phrase)} "
    return needle in normalised


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def request_with_retries(
    session: requests.Session,
    url: str,
    *,
    params: dict[str, Any],
    timeout: int = 90,
    attempts: int = 6,
) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt + 1 == attempts:
                break
            time.sleep(min(30, 2**attempt))
    raise RuntimeError(f"Request failed after {attempts} attempts: {url}") from last_error


def initialise_database(path: Path, *, refresh: bool) -> sqlite3.Connection:
    if refresh and path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS studies (
            nct_id TEXT PRIMARY KEY,
            brief_title TEXT,
            overall_status TEXT,
            study_type TEXT,
            phases_json TEXT,
            start_date TEXT,
            completion_date TEXT,
            first_posted TEXT,
            last_updated TEXT
        );
        CREATE TABLE IF NOT EXISTS study_pairs (
            drug TEXT NOT NULL,
            disease TEXT NOT NULL,
            nct_id TEXT NOT NULL,
            raw_drug TEXT NOT NULL,
            raw_disease TEXT NOT NULL,
            first_posted TEXT,
            start_date TEXT,
            PRIMARY KEY (drug, disease, nct_id, raw_drug, raw_disease)
        );
        CREATE INDEX IF NOT EXISTS idx_study_pairs_drug
            ON study_pairs(drug);
        CREATE INDEX IF NOT EXISTS idx_study_pairs_disease
            ON study_pairs(disease);
        CREATE INDEX IF NOT EXISTS idx_study_pairs_target
            ON study_pairs(drug, disease);
        """
    )
    return connection


def metadata_get(connection: sqlite3.Connection, key: str) -> str | None:
    row = connection.execute("SELECT value FROM metadata WHERE key = ?", (key,)).fetchone()
    return None if row is None else str(row[0])


def metadata_set(connection: sqlite3.Connection, key: str, value: Any) -> None:
    connection.execute(
        "INSERT INTO metadata(key, value) VALUES (?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
        (key, json.dumps(value, ensure_ascii=False)),
    )


def unpack_metadata(value: str | None, default: Any = None) -> Any:
    if value is None:
        return default
    return json.loads(value)


def date_value(module: dict[str, Any], key: str) -> str:
    return str((module.get(key) or {}).get("date") or "")


def extract_study(study: dict[str, Any]) -> tuple[tuple[Any, ...] | None, list[tuple[Any, ...]], Counter]:
    audit: Counter = Counter()
    protocol = study.get("protocolSection") or {}
    identification = protocol.get("identificationModule") or {}
    status = protocol.get("statusModule") or {}
    design = protocol.get("designModule") or {}
    conditions_module = protocol.get("conditionsModule") or {}
    arms = protocol.get("armsInterventionsModule") or {}

    nct_id = str(identification.get("nctId") or "").strip()
    if not nct_id:
        audit["missing_nct_id"] += 1
        return None, [], audit

    study_type = str(design.get("studyType") or "")
    study_row = (
        nct_id,
        str(identification.get("briefTitle") or ""),
        str(status.get("overallStatus") or ""),
        study_type,
        json.dumps(design.get("phases") or []),
        date_value(status, "startDateStruct"),
        date_value(status, "completionDateStruct"),
        date_value(status, "studyFirstPostDateStruct"),
        date_value(status, "lastUpdatePostDateStruct"),
    )
    if study_type != "INTERVENTIONAL":
        audit["non_interventional"] += 1
        return study_row, [], audit

    raw_conditions = [
        str(value).strip()
        for value in conditions_module.get("conditions") or []
        if str(value).strip()
    ]
    interventions = arms.get("interventions") or []
    raw_drugs = [
        str(item.get("name") or "").strip()
        for item in interventions
        if str(item.get("type") or "").upper() == "DRUG"
        and str(item.get("name") or "").strip()
    ]
    if not raw_conditions:
        audit["interventional_missing_condition"] += 1
    if not raw_drugs:
        audit["interventional_without_drug"] += 1

    # Deduplicate repeated spelling within a study while retaining raw provenance.
    condition_values = sorted({(normalise_exact(v), v) for v in raw_conditions if normalise_exact(v)})
    drug_values = sorted({(normalise_exact(v), v) for v in raw_drugs if normalise_exact(v)})
    pair_rows = [
        (
            drug,
            disease,
            nct_id,
            raw_drug,
            raw_disease,
            study_row[7],
            study_row[5],
        )
        for drug, raw_drug in drug_values
        for disease, raw_disease in condition_values
    ]
    audit["interventional_with_drug"] += bool(raw_drugs)
    audit["study_pair_rows"] += len(pair_rows)
    return study_row, pair_rows, audit


def download_registry(
    connection: sqlite3.Connection,
    session: requests.Session,
    *,
    page_size: int,
    max_pages: int | None,
    pause: float,
) -> dict[str, Any]:
    if unpack_metadata(metadata_get(connection, "registry_complete"), False):
        print("ClinicalTrials.gov snapshot already complete; reusing database.", flush=True)
        return unpack_metadata(metadata_get(connection, "registry_audit"), {})

    page_token = unpack_metadata(metadata_get(connection, "next_page_token"))
    page_number = int(unpack_metadata(metadata_get(connection, "pages_downloaded"), 0))
    audit = Counter(unpack_metadata(metadata_get(connection, "registry_audit"), {}))
    if not metadata_get(connection, "registry_started_at"):
        metadata_set(connection, "registry_started_at", utc_now())
        connection.commit()

    while True:
        if max_pages is not None and page_number >= max_pages:
            print(f"Stopped at --max-pages={max_pages}; snapshot is incomplete.", flush=True)
            break
        params: dict[str, Any] = {
            "format": "json",
            "pageSize": page_size,
            "fields": CTGOV_FIELDS,
        }
        if page_token:
            params["pageToken"] = page_token
        result = request_with_retries(session, CTGOV_URL, params=params).json()
        studies = result.get("studies") or []

        page_audit: Counter = Counter()
        with connection:
            for study in studies:
                study_row, pair_rows, study_audit = extract_study(study)
                page_audit.update(study_audit)
                if study_row is not None:
                    connection.execute(
                        "INSERT OR REPLACE INTO studies VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        study_row,
                    )
                if pair_rows:
                    connection.executemany(
                        "INSERT OR IGNORE INTO study_pairs VALUES (?, ?, ?, ?, ?, ?, ?)",
                        pair_rows,
                    )
            page_number += 1
            audit.update(page_audit)
            audit["studies_returned"] += len(studies)
            page_token = result.get("nextPageToken")
            metadata_set(connection, "pages_downloaded", page_number)
            metadata_set(connection, "next_page_token", page_token)
            metadata_set(connection, "registry_audit", dict(audit))

        if page_number == 1 or page_number % 25 == 0 or not page_token:
            print(
                f"Registry page {page_number}: {audit['studies_returned']:,} studies, "
                f"{audit['study_pair_rows']:,} exact study-pair rows seen",
                flush=True,
            )
        if not page_token:
            with connection:
                metadata_set(connection, "registry_complete", True)
                metadata_set(connection, "registry_completed_at", utc_now())
            break
        if pause:
            time.sleep(pause)

    return dict(audit)


def build_pair_edges(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        DROP TABLE IF EXISTS pair_edges;
        CREATE TABLE pair_edges AS
        SELECT drug, disease,
               COUNT(DISTINCT nct_id) AS trial_count,
               MIN(NULLIF(first_posted, '')) AS first_posted,
               MAX(NULLIF(first_posted, '')) AS last_posted
        FROM study_pairs
        GROUP BY drug, disease;
        CREATE UNIQUE INDEX idx_pair_edges_pair ON pair_edges(drug, disease);
        CREATE INDEX idx_pair_edges_drug ON pair_edges(drug);
        CREATE INDEX idx_pair_edges_disease ON pair_edges(disease);
        """
    )
    connection.commit()


def target_trial_rows(
    connection: sqlite3.Connection, drug: str, disease: str
) -> list[dict[str, Any]]:
    query = """
        SELECT DISTINCT p.nct_id, s.brief_title, s.overall_status, s.start_date,
               s.first_posted, s.last_updated, p.raw_drug, p.raw_disease
        FROM study_pairs p
        JOIN studies s ON s.nct_id = p.nct_id
        WHERE p.drug = ? AND p.disease = ?
        ORDER BY s.first_posted, p.nct_id
    """
    columns = (
        "nct_id",
        "brief_title",
        "overall_status",
        "start_date",
        "first_posted",
        "last_updated",
        "raw_drug",
        "raw_disease",
    )
    return [dict(zip(columns, row)) for row in connection.execute(query, (drug, disease))]


def apparent_pre_2020_disease_rows(
    connection: sqlite3.Connection, disease: str
) -> list[dict[str, Any]]:
    """Find current COVID fields attached to studies originally posted pre-2020.

    These rows demonstrate why the current registry cannot be converted into a
    historical snapshot using a study-level date: conditions/interventions can
    be added by later record amendments.
    """
    query = """
        SELECT DISTINCT p.nct_id, s.brief_title, p.raw_drug, p.raw_disease,
               s.first_posted, s.last_updated
        FROM study_pairs p
        JOIN studies s ON s.nct_id = p.nct_id
        WHERE p.disease = ? AND s.first_posted < '2020-01-01'
        ORDER BY s.first_posted, p.nct_id, p.raw_drug
    """
    columns = (
        "nct_id",
        "brief_title",
        "raw_drug",
        "raw_disease",
        "study_first_posted",
        "record_last_updated",
    )
    return [dict(zip(columns, row)) for row in connection.execute(query, (disease,))]


def cutoff_clause(cutoff_exclusive: str | None, alias: str) -> tuple[str, list[str]]:
    if cutoff_exclusive is None:
        return "", []
    return f" AND {alias}.first_posted < ?", [cutoff_exclusive]


def alternative_paths(
    connection: sqlite3.Connection,
    drug: str,
    disease: str,
    *,
    cutoff_exclusive: str | None = None,
) -> list[tuple[str, str]]:
    # In a bipartite graph, after withholding the direct edge, the first possible
    # target path has length three: drug -> other disease -> other drug -> disease.
    cutoff1, params1 = cutoff_clause(cutoff_exclusive, "e1")
    cutoff2, params2 = cutoff_clause(cutoff_exclusive, "e2")
    cutoff3, params3 = cutoff_clause(cutoff_exclusive, "e3")
    query = (
        """
            SELECT e1.disease, e2.drug
            FROM pair_edges e1
            JOIN pair_edges e2 ON e2.disease = e1.disease
            JOIN pair_edges e3 ON e3.drug = e2.drug
            WHERE e1.drug = ? AND e3.disease = ?
              AND e1.disease <> ? AND e2.drug <> ?
            """
        + cutoff1
        + cutoff2
        + cutoff3
        + """
            GROUP BY e1.disease, e2.drug
            ORDER BY e1.disease, e2.drug
            """
    )
    params = [drug, disease, disease, drug, *params1, *params2, *params3]
    return list(connection.execute(query, params))


def node_degree(
    connection: sqlite3.Connection,
    *,
    node_type: str,
    value: str,
    withheld_drug: str,
    withheld_disease: str,
    cutoff_exclusive: str | None = None,
) -> int:
    cutoff, cutoff_params = cutoff_clause(cutoff_exclusive, "pair_edges")
    if node_type == "drug":
        row = connection.execute(
            "SELECT COUNT(*) FROM pair_edges WHERE drug = ? AND disease <> ?" + cutoff,
            [value, withheld_disease if value == withheld_drug else "\0", *cutoff_params],
        ).fetchone()
    else:
        row = connection.execute(
            "SELECT COUNT(*) FROM pair_edges WHERE disease = ? AND drug <> ?" + cutoff,
            [value, withheld_drug if value == withheld_disease else "\0", *cutoff_params],
        ).fetchone()
    return int(row[0])


def structural_snapshot(
    connection: sqlite3.Connection,
    *,
    drug: str,
    disease: str,
    cutoff_exclusive: str | None,
    include_path_rows: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    paths = alternative_paths(
        connection, drug, disease, cutoff_exclusive=cutoff_exclusive
    )
    drug_degree = node_degree(
        connection,
        node_type="drug",
        value=drug,
        withheld_drug=drug,
        withheld_disease=disease,
        cutoff_exclusive=cutoff_exclusive,
    )
    disease_degree = node_degree(
        connection,
        node_type="disease",
        value=disease,
        withheld_drug=drug,
        withheld_disease=disease,
        cutoff_exclusive=cutoff_exclusive,
    )
    degree_cache: dict[tuple[str, str], int] = {}

    def degree(kind: str, value: str) -> int:
        key = (kind, value)
        if key not in degree_cache:
            degree_cache[key] = node_degree(
                connection,
                node_type=kind,
                value=value,
                withheld_drug=drug,
                withheld_disease=disease,
                cutoff_exclusive=cutoff_exclusive,
            )
        return degree_cache[key]

    random_walk_3 = 0.0
    path_rows: list[dict[str, Any]] = []
    for middle_disease, middle_drug in paths:
        denominator = drug_degree * degree("disease", middle_disease) * degree("drug", middle_drug)
        contribution = 0.0 if denominator == 0 else 1.0 / denominator
        random_walk_3 += contribution
        if include_path_rows:
            path_rows.append(
                {
                    "source_drug": drug,
                    "middle_disease": middle_disease,
                    "middle_drug": middle_drug,
                    "target_disease": disease,
                    "path_length": 3,
                    "random_walk_contribution": contribution,
                }
            )
    snapshot = {
        "cutoff_exclusive": cutoff_exclusive,
        "target_edge_withheld": True,
        "drug_degree": drug_degree,
        "disease_degree": disease_degree,
        "alternative_path_count_length_3": len(paths),
        "shortest_path": 3 if paths else None,
        "inverse_shortest_path": (1.0 / 3.0) if paths else 0.0,
        "three_step_random_walk_probability": random_walk_3,
        "preferential_attachment": drug_degree * disease_degree,
    }
    return snapshot, path_rows


def compute_network_audit(
    connection: sqlite3.Connection, *, drug_label: str, disease_label: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    drug = normalise_exact(drug_label)
    disease = normalise_exact(disease_label)
    direct = connection.execute(
        "SELECT trial_count, first_posted, last_posted FROM pair_edges "
        "WHERE drug = ? AND disease = ?",
        (drug, disease),
    ).fetchone()
    current_snapshot, path_rows = structural_snapshot(
        connection,
        drug=drug,
        disease=disease,
        cutoff_exclusive=None,
        include_path_rows=True,
    )
    temporal_snapshots = {}
    for label, cutoff in (
        ("pre_2020", "2020-01-01"),
        ("before_2020_03_01", "2020-03-01"),
        ("before_2020_06_01", "2020-06-01"),
    ):
        temporal_snapshots[label], _ = structural_snapshot(
            connection,
            drug=drug,
            disease=disease,
            cutoff_exclusive=cutoff,
            include_path_rows=False,
        )

    counts = {
        "registry_studies": connection.execute("SELECT COUNT(*) FROM studies").fetchone()[0],
        "study_pair_rows": connection.execute("SELECT COUNT(*) FROM study_pairs").fetchone()[0],
        "unique_drug_disease_edges": connection.execute("SELECT COUNT(*) FROM pair_edges").fetchone()[0],
        "unique_drug_nodes": connection.execute("SELECT COUNT(DISTINCT drug) FROM pair_edges").fetchone()[0],
        "unique_disease_nodes": connection.execute("SELECT COUNT(DISTINCT disease) FROM pair_edges").fetchone()[0],
    }
    target_trials = target_trial_rows(connection, drug, disease)
    target_first_posted_before_2020 = sum(
        bool(row["first_posted"] and row["first_posted"] < "2020-01-01")
        for row in target_trials
    )
    target_first_posted_after_2019 = sum(
        bool(row["first_posted"] and row["first_posted"] >= "2020-01-01")
        for row in target_trials
    )
    amended_legacy_rows = apparent_pre_2020_disease_rows(connection, disease)
    metrics = {
        "target": {
            "drug_supplied": drug_label,
            "disease_supplied": disease_label,
            "drug_normalised": drug,
            "disease_normalised": disease,
        },
        "matching_rule": "whole-field equality after Unicode casefold and punctuation/whitespace normalisation; no synonyms, substring, fuzzy match, or ontology expansion",
        "graph_definition": "bipartite study-co-occurrence graph from INTERVENTIONAL studies and interventions explicitly typed DRUG; every drug field is crossed with every condition field in that study",
        "counts": counts,
        "direct_edge_present_in_snapshot": direct is not None,
        "direct_edge_exact_trial_count": 0 if direct is None else int(direct[0]),
        "direct_edge_first_posted": None if direct is None else direct[1],
        "direct_edge_last_posted": None if direct is None else direct[2],
        "target_trial_ids": [row["nct_id"] for row in target_trials],
        "target_trials_first_posted_before_2020": target_first_posted_before_2020,
        "target_trials_first_posted_after_2019": target_first_posted_after_2019,
        "leakage_control": "all exact HCQ--COVID study rows withheld before structural calculation",
        "withheld_graph_drug_degree": current_snapshot["drug_degree"],
        "withheld_graph_disease_degree": current_snapshot["disease_degree"],
        "withheld_graph_alternative_path_count_length_3": current_snapshot[
            "alternative_path_count_length_3"
        ],
        "withheld_graph_shortest_path": current_snapshot["shortest_path"],
        "withheld_graph_inverse_shortest_path": current_snapshot["inverse_shortest_path"],
        "withheld_graph_three_step_random_walk_probability": current_snapshot[
            "three_step_random_walk_probability"
        ],
        "withheld_graph_preferential_attachment": current_snapshot[
            "preferential_attachment"
        ],
        "temporal_withheld_graph_snapshots": temporal_snapshots,
        "historical_snapshot_warning": (
            "The temporal snapshots filter current record content by the study's "
            "original first-posted date. They are not valid historical registry "
            "snapshots because later amendments can add COVID-19 fields to older studies."
        ),
        "current_covid_rows_on_pre_2020_studies": len(amended_legacy_rows),
        "current_covid_pre_2020_study_ids": sorted(
            {row["nct_id"] for row in amended_legacy_rows}
        ),
        "interpretation": "structural relatedness/co-registration only; not biological validation or clinical efficacy",
    }
    return metrics, path_rows


def iter_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())


def first_text(element: ET.Element, paths: Iterable[str]) -> str:
    for path in paths:
        value = element.findtext(path)
        if value:
            return " ".join(value.split())
    return ""


def parse_pubmed_article(article: ET.Element, drug: str, disease: str) -> dict[str, Any]:
    medline = article.find("MedlineCitation")
    citation = None if medline is None else medline.find("Article")
    journal_issue = None if citation is None else citation.find("Journal/JournalIssue")
    pubdate = None if journal_issue is None else journal_issue.find("PubDate")
    title = "" if citation is None else iter_text(citation.find("ArticleTitle"))
    abstract_parts = [] if citation is None else citation.findall("Abstract/AbstractText")
    abstract = " ".join(iter_text(part) for part in abstract_parts).strip()
    combined = f"{title} {abstract}"
    pmid = "" if medline is None else str(medline.findtext("PMID") or "")

    year = "" if pubdate is None else str(pubdate.findtext("Year") or "")
    month = "" if pubdate is None else str(pubdate.findtext("Month") or "")
    day = "" if pubdate is None else str(pubdate.findtext("Day") or "")
    medline_date = "" if pubdate is None else str(pubdate.findtext("MedlineDate") or "")
    publication_date = "-".join(value for value in (year, month, day) if value) or medline_date

    article_ids = {
        str(item.attrib.get("IdType") or ""): str(item.text or "")
        for item in article.findall("PubmedData/ArticleIdList/ArticleId")
    }
    publication_types = [] if citation is None else [
        iter_text(item) for item in citation.findall("PublicationTypeList/PublicationType")
    ]
    mesh_terms = [] if medline is None else [
        iter_text(item.find("DescriptorName")) for item in medline.findall("MeshHeadingList/MeshHeading")
    ]
    return {
        "record_type": "PubmedArticle",
        "pmid": pmid,
        "doi": article_ids.get("doi", ""),
        "title": title,
        "abstract": abstract,
        "journal": "" if citation is None else first_text(citation, ("Journal/Title",)),
        "publication_date": publication_date,
        "publication_types": publication_types,
        "mesh_terms": mesh_terms,
        "exact_drug_phrase_verified": contains_exact_phrase(combined, drug),
        "exact_disease_phrase_verified": contains_exact_phrase(combined, disease),
    }


def parse_pubmed_book_article(
    article: ET.Element, drug: str, disease: str
) -> dict[str, Any]:
    title = iter_text(article.find(".//ArticleTitle")) or iter_text(
        article.find(".//BookTitle")
    )
    abstract = " ".join(
        iter_text(part) for part in article.findall(".//AbstractText")
    ).strip()
    combined = f"{title} {abstract}"
    article_ids = {
        str(item.attrib.get("IdType") or ""): str(item.text or "")
        for item in article.findall(".//ArticleId")
    }
    return {
        "record_type": "PubmedBookArticle",
        "pmid": str(article.findtext(".//PMID") or ""),
        "doi": article_ids.get("doi", ""),
        "title": title,
        "abstract": abstract,
        "journal": iter_text(article.find(".//BookTitle")),
        "publication_date": str(article.findtext(".//PubDate/Year") or ""),
        "publication_types": [],
        "mesh_terms": [],
        "exact_drug_phrase_verified": contains_exact_phrase(combined, drug),
        "exact_disease_phrase_verified": contains_exact_phrase(combined, disease),
    }


def pubmed_query(drug: str, disease: str) -> str:
    # Phrase verification below guards against any search-engine term expansion.
    return f'"{drug}"[Title/Abstract] AND "{disease}"[Title/Abstract]'


def fetch_pubmed(
    session: requests.Session,
    output_dir: Path,
    *,
    drug: str,
    disease: str,
    batch_size: int,
    pause: float,
    refresh: bool,
) -> dict[str, Any]:
    records_path = output_dir / "pubmed_records.jsonl.gz"
    exact_records_path = output_dir / "pubmed_exact_pair_records.jsonl.gz"
    audit_path = output_dir / "pubmed_audit.json"
    if records_path.exists() and audit_path.exists() and not refresh:
        print("PubMed records already present; reusing them.", flush=True)
        return json.loads(audit_path.read_text(encoding="utf-8"))

    query = pubmed_query(drug, disease)
    search_response = request_with_retries(
        session,
        NCBI_ESEARCH_URL,
        params={
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": 10000,
            "usehistory": "y",
            "sort": "pub date",
        },
    ).json()["esearchresult"]
    count = int(search_response["count"])
    if count > 10000:
        raise RuntimeError(
            f"PubMed query returned {count:,} records; pagination beyond 10,000 is not implemented."
        )
    ids = list(search_response.get("idlist") or [])
    if len(ids) != count:
        raise RuntimeError(f"PubMed ESearch returned {len(ids)} IDs for a reported count of {count}.")

    exact_both = 0
    exact_drug = 0
    exact_disease = 0
    parsed = 0
    parsed_ids: set[str] = set()
    record_types: Counter = Counter()
    publication_years: Counter = Counter()
    exact_publication_years: Counter = Counter()
    with gzip.open(records_path, "wt", encoding="utf-8", newline="") as handle, gzip.open(
        exact_records_path, "wt", encoding="utf-8", newline=""
    ) as exact_handle:
        for start in range(0, len(ids), batch_size):
            batch = ids[start : start + batch_size]
            xml_bytes = request_with_retries(
                session,
                NCBI_EFETCH_URL,
                params={"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
            ).content
            root = ET.fromstring(xml_bytes)
            containers = [
                (article, parse_pubmed_article)
                for article in root.findall("PubmedArticle")
            ] + [
                (article, parse_pubmed_book_article)
                for article in root.findall("PubmedBookArticle")
            ]
            for article, parser in containers:
                record = parser(article, drug, disease)
                parsed += 1
                if record["pmid"]:
                    parsed_ids.add(record["pmid"])
                record_types[record["record_type"]] += 1
                year_match = re.search(r"(?:19|20)\d{2}", record["publication_date"])
                publication_year = year_match.group(0) if year_match else "unknown"
                publication_years[publication_year] += 1
                exact_drug += int(record["exact_drug_phrase_verified"])
                exact_disease += int(record["exact_disease_phrase_verified"])
                exact_both += int(
                    record["exact_drug_phrase_verified"]
                    and record["exact_disease_phrase_verified"]
                )
                handle.write(json.dumps(record, ensure_ascii=False) + "\n")
                if (
                    record["exact_drug_phrase_verified"]
                    and record["exact_disease_phrase_verified"]
                ):
                    exact_publication_years[publication_year] += 1
                    exact_handle.write(json.dumps(record, ensure_ascii=False) + "\n")
            print(f"PubMed: parsed {parsed:,}/{count:,} records", flush=True)
            if pause:
                time.sleep(max(pause, 0.34))

    unresolved_ids = sorted(set(ids) - parsed_ids)
    write_json(output_dir / "pubmed_search_ids.json", ids)
    write_json(output_dir / "pubmed_unresolved_ids.json", unresolved_ids)
    audit = {
        "retrieved_at": utc_now(),
        "database": "PubMed",
        "query": query,
        "query_translation": search_response.get("querytranslation"),
        "reported_count": count,
        "ids_returned": len(ids),
        "records_parsed": parsed,
        "record_types": dict(record_types),
        "unresolved_id_count": len(unresolved_ids),
        "unresolved_ids": unresolved_ids,
        "exact_drug_phrase_verified": exact_drug,
        "exact_disease_phrase_verified": exact_disease,
        "exact_both_phrases_verified": exact_both,
        "publication_year_counts_all_query_results": dict(
            sorted(publication_years.items())
        ),
        "publication_year_counts_exact_pair": dict(
            sorted(exact_publication_years.items())
        ),
        "date_restriction": None,
        "result_cap": None,
        "evidence_direction_classification": "not performed",
        "classification_reason": "publication retrieval and pair verification do not establish therapeutic direction, efficacy, or safety",
        "records_file": records_path.name,
        "exact_pair_records_file": exact_records_path.name,
    }
    write_json(audit_path, audit)
    return audit


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_summary(
    *,
    registry_audit: dict[str, Any],
    network: dict[str, Any],
    pubmed: dict[str, Any] | None,
    database_path: Path,
) -> str:
    direct_count = network["direct_edge_exact_trial_count"]
    first_posted = network["direct_edge_first_posted"] or "not present"
    lines = [
        "# Exact-pair full-registry audit: hydroxychloroquine × COVID-19",
        "",
        f"Run completed: {utc_now()}",
        "",
        "## ClinicalTrials.gov and network",
        "",
        f"- Registry studies stored: {network['counts']['registry_studies']:,}",
        f"- Interventional drug–condition study rows: {network['counts']['study_pair_rows']:,}",
        f"- Unique exact drug–condition edges: {network['counts']['unique_drug_disease_edges']:,}",
        f"- Exact HCQ–COVID trials: {direct_count:,}",
        f"- Earliest posting among exact target trials: {first_posted}",
        f"- Exact target trials first posted before 2020: {network['target_trials_first_posted_before_2020']:,}",
        f"- Alternative length-3 paths after withholding the target edge: {network['withheld_graph_alternative_path_count_length_3']:,}",
        f"- Shortest withheld-graph path: {network['withheld_graph_shortest_path']}",
        f"- Three-step random-walk probability: {network['withheld_graph_three_step_random_walk_probability']:.12g}",
        "",
        "The direct HCQ–COVID edge and every supporting NCT record were excluded from structural scoring. "
        "The resulting metrics describe registry co-occurrence, not efficacy.",
        "",
        "Temporal warning: current registry fields include later amendments. "
        f"There are {network['current_covid_rows_on_pre_2020_studies']:,} current COVID-labelled drug rows "
        "on studies first posted before 2020, so filtering current records by first-posted date is not a valid historical snapshot.",
        "",
        "## PubMed",
        "",
    ]
    if pubmed is None:
        lines.append("PubMed retrieval was skipped.")
    else:
        lines.extend(
            [
                f"- Declared query: `{pubmed['query']}`",
                f"- Records returned: {pubmed['reported_count']:,}",
                f"- Records with both phrases verified in title/abstract: {pubmed['exact_both_phrases_verified']:,}",
                "- Therapeutic/adverse/irrelevant counts: not computed; retrieval volume alone cannot determine evidence direction.",
            ]
        )
    lines.extend(
        [
            "",
            "## Reproducibility files",
            "",
            f"- SQLite registry graph: `{database_path.name}`",
            "- Exact target trial table: `target_trials.csv`",
            "- Alternative paths: `target_alternative_paths.csv`",
            "- Machine-readable metrics: `network_audit.json`",
            "- Retrieval manifest: `run_manifest.json`",
            "",
            f"Registry page audit: {json.dumps(registry_audit, sort_keys=True)}",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drug", default="hydroxychloroquine")
    parser.add_argument("--disease", default="COVID-19")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "exact_pair_hcq_covid_current",
    )
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--pubmed-batch-size", type=int, default=200)
    parser.add_argument("--pause", type=float, default=0.05)
    parser.add_argument("--max-pages", type=int)
    parser.add_argument("--skip-pubmed", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--refresh-pubmed",
        action="store_true",
        help="refresh PubMed outputs while reusing a completed registry database",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    database_path = output_dir / "clinicaltrials_exact_graph.sqlite"
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    connection = initialise_database(database_path, refresh=args.refresh)
    try:
        registry_audit = download_registry(
            connection,
            session,
            page_size=args.page_size,
            max_pages=args.max_pages,
            pause=args.pause,
        )
        registry_complete = unpack_metadata(metadata_get(connection, "registry_complete"), False)
        if not registry_complete:
            print("Partial registry saved; rerun without --max-pages to resume.", flush=True)
            return 2

        print("Aggregating exact study pairs into graph edges...", flush=True)
        build_pair_edges(connection)
        network, path_rows = compute_network_audit(
            connection, drug_label=args.drug, disease_label=args.disease
        )
        trial_rows = target_trial_rows(
            connection, normalise_exact(args.drug), normalise_exact(args.disease)
        )
        write_json(output_dir / "network_audit.json", network)
        write_csv(output_dir / "target_trials.csv", trial_rows)
        write_csv(output_dir / "target_alternative_paths.csv", path_rows)
        write_csv(
            output_dir / "temporal_leakage_records.csv",
            apparent_pre_2020_disease_rows(
                connection, normalise_exact(args.disease)
            ),
        )

        pubmed = None
        if not args.skip_pubmed:
            pubmed = fetch_pubmed(
                session,
                output_dir,
                drug=args.drug,
                disease=args.disease,
                batch_size=args.pubmed_batch_size,
                pause=args.pause,
                refresh=args.refresh or args.refresh_pubmed,
            )

        manifest = {
            "run_completed_at": utc_now(),
            "script": str(Path(__file__).resolve()),
            "clinicaltrials_api": CTGOV_URL,
            "clinicaltrials_fields": CTGOV_FIELDS.split(","),
            "clinicaltrials_snapshot_started_at": unpack_metadata(
                metadata_get(connection, "registry_started_at")
            ),
            "clinicaltrials_snapshot_completed_at": unpack_metadata(
                metadata_get(connection, "registry_completed_at")
            ),
            "registry_complete": registry_complete,
            "registry_audit": registry_audit,
            "target": {"drug": args.drug, "disease": args.disease},
            "matching_rule": network["matching_rule"],
            "pubmed": pubmed,
        }
        write_json(output_dir / "run_manifest.json", manifest)
        (output_dir / "SUMMARY.md").write_text(
            make_summary(
                registry_audit=registry_audit,
                network=network,
                pubmed=pubmed,
                database_path=database_path,
            ),
            encoding="utf-8",
        )
        print(f"Complete. Results: {output_dir}", flush=True)
        return 0
    finally:
        connection.close()


if __name__ == "__main__":
    sys.exit(main())
