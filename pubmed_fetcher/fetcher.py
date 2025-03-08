import requests
import pandas as pd
import re

SEARCH_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
FETCH_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

def fetch_paper_ids(query: str, max_results: int = 10):
    """Fetch PubMed paper IDs based on a search query."""
    params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": max_results}
    response = requests.get(SEARCH_API, params=params)
    response.raise_for_status()
    return response.json().get("esearchresult", {}).get("idlist", [])

def fetch_paper_details(paper_ids):
    """Fetch paper details including title, publication date, and authors."""
    if not paper_ids:
        return []

    params = {"db": "pubmed", "id": ",".join(paper_ids), "retmode": "json"}
    response = requests.get(FETCH_API, params=params)
    response.raise_for_status()
    data = response.json().get("result", {})

    papers = []
    for paper_id in paper_ids:
        paper = data.get(paper_id, {})
        authors = paper.get("authors", [])
        company_authors, company_affiliations = extract_company_authors(authors)

        papers.append({
            "PubmedID": paper_id,
            "Title": paper.get("title", "N/A"),
            "Publication Date": paper.get("pubdate", "N/A"),
            "Non-academic Author(s)": ", ".join(company_authors) if company_authors else "N/A",
            "Company Affiliation(s)": ", ".join(company_affiliations) if company_affiliations else "N/A",
            "Corresponding Author Email": extract_corresponding_author_email(paper.get("elocationid", ""))
        })

    return papers

def extract_company_authors(authors):
    """Identify authors affiliated with non-academic institutions."""
    company_authors = []
    company_affiliations = []

    for author in authors:
        name = author.get("name", "Unknown")
        affiliation = author.get("affiliation", "")

        if affiliation and not re.search(r"university|college|institute|school", affiliation, re.I):
            company_authors.append(name)
            company_affiliations.append(affiliation)

    return company_authors, company_affiliations

def extract_corresponding_author_email(elocationid):
    """Extract email if available (heuristic method)."""
    email_match = re.search(r"([\w\.-]+@[\w\.-]+\.\w+)", elocationid)
    return email_match.group(1) if email_match else "N/A"

def save_to_csv(data, filename="results.csv"):
    """Save fetched data to a CSV file."""
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
