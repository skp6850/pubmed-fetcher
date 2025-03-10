import logging
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import re

# Set your email (mandatory for using Entrez)
Entrez.email = "use_your_id@gmail.com"
API_KEY = "API_FROM_NCBI"  # Optional but recommended

def search_pubmed(query: str, start_year: int, end_year: int = None, max_results: int = 100) -> list:
    """Search PubMed for articles published in a specific year or date range and return a list of PMIDs."""
    try:
        if not end_year:
            end_year = start_year  # If only one year is provided, search for that year only.

        date_filter = f'("{start_year}/01/01"[Date - Publication] : "{end_year}/12/31"[Date - Publication])'
        full_query = f"{query} AND {date_filter}"  # Adds date filter to query

        params = {
            "db": "pubmed",
            "term": full_query,
            "retmax": max_results,  # Fetch up to 100 results
            "retmode": "xml",
            "api_key": API_KEY
        }
        handle = Entrez.esearch(**params)
        records = Entrez.read(handle)

        return records["IdList"]

    except Exception as e:
        logging.error(f"Error while fetching data: {e}")
        return []



def fetch_article_details(pubmed_ids: list) -> list:
    """Fetch detailed information for a list of PubMed IDs."""
    if not pubmed_ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml", api_key=API_KEY)
    articles = ET.parse(handle).getroot()

    results = []
    for article in articles.findall(".//PubmedArticle"):
        pmid = article.find(".//PMID").text
        title = article.find(".//ArticleTitle").text if article.find(".//ArticleTitle") is not None else "N/A"
        pub_date = article.find(".//PubDate/Year").text if article.find(".//PubDate/Year") is not None else "N/A"
        authors, affiliations = extract_authors_and_affiliations(article)

        results.append({
            "PubmedID": pmid,
            "Title": title,
            "Publication Date": pub_date,
            "Non-academic Author(s)": ", ".join(authors) if authors else "N/A",
            "Company Affiliation(s)": ", ".join(affiliations) if affiliations else "N/A",
            "Corresponding Author Email": extract_corresponding_author_email(article) or "N/A"
        })
    return results

def extract_authors_and_affiliations(article) -> tuple:
    """Extract non-academic authors and their company affiliations."""
    authors = []
    affiliations = set()

    # Define academic keywords to exclude
    ACADEMIC_KEYWORDS = ["university", "college", "hospital", "institute", "school", "research center"]

    for author in article.findall(".//Author"):
        last_name = author.find("LastName")
        fore_name = author.find("ForeName")
        full_name = f"{fore_name.text} {last_name.text}" if fore_name is not None and last_name is not None else "Unknown"

        aff_element = author.find(".//Affiliation")
        affiliation = aff_element.text.lower() if aff_element is not None else ""

        # Exclude authors if their affiliation contains academic keywords
        if affiliation and not any(keyword in affiliation for keyword in ACADEMIC_KEYWORDS):
            authors.append(full_name)
            affiliations.add(affiliation)

    return authors, affiliations


def extract_corresponding_author_email(article) -> str:
    """Extract corresponding author's email if available."""
    for affiliation in article.findall(".//Affiliation"):
        email_match = re.search(r"([\w\.-]+@[\w\.-]+\.\w+)", affiliation.text if affiliation.text else "")
        if email_match:
            return email_match.group(1)
    return ""

def save_to_csv(data: list, filename: str = "results.csv"):
    """Save results to a CSV file."""
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
