import requests

def pubmed_count_requests(drug, disease):
    term = f'"{drug}"[Title/Abstract] AND "{disease}"[Title/Abstract]'
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "email": "oseifrancis633@gmail.com" # Replace with your email
    }

    r = requests.get(url, params=params, timeout=10)
    r.raise_for_status()
    return int(r.json()["esearchresult"]["count"])
