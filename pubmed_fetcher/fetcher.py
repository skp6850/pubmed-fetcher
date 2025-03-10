import logging
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import re

# Set your email (mandatory for using Entrez)
Entrez.email = "skp68500@gmail.com"
API_KEY = "a30a3cdded67bc4c65592abfad40aa57b509"  # Optional but recommended

def search_pubmed(query: str, start_year: int, end_year: int = None, max_results: int = 100) -> list:
    """Search PubMed for English-language articles published in a specific year or date range and return a unique list of PMIDs."""
    try:
        if not end_year:
            end_year = start_year  # If only one year is provided, search for that year only.

        date_filter = f'("{start_year}/01/01"[Date - Publication] : "{end_year}/12/31"[Date - Publication])'
        # Append language filter for English
        full_query = f"{query} AND {date_filter} AND eng[Language]"
        
        params = {
            "db": "pubmed",
            "term": full_query,
            "retmax": max_results,
            "retmode": "xml",
            "api_key": API_KEY
        }
        handle = Entrez.esearch(**params)
        record = Entrez.read(handle)
        handle.close()
        # Ensure unique PMIDs while preserving order
        unique_pmids = list(dict.fromkeys(record.get("IdList", [])))
        return unique_pmids

    except Exception as e:
        logging.error(f"Error while fetching data: {e}")
        return []

def fetch_article_details(pubmed_ids: list) -> list:
    """Fetch detailed information for a list of PubMed IDs and return unique results."""
    if not pubmed_ids:
        return []

    try:
        handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml", api_key=API_KEY)
        data = handle.read()
        handle.close()
        results = parse_articles_xml(data)
        
        # Filter duplicates based on PubmedID
        seen = set()
        unique_results = []
        for article in results:
            if article["PubmedID"] not in seen:
                seen.add(article["PubmedID"])
                unique_results.append(article)
        return unique_results
    except Exception as e:
        logging.error(f"Fetching articles failed: {e}")
        raise

def parse_articles_xml(xml_data: str) -> list:
    """Parse PubMed XML data into a list of ArticleData dictionaries."""
    root = ET.fromstring(xml_data)
    articles = []

    for article_elem in root.findall('.//PubmedArticle'):
        # Check language; skip non-English articles if Language tag exists.
        lang_elem = article_elem.find('.//Language')
        if lang_elem is not None and lang_elem.text.lower() != "eng":
            continue

        pmid_elem = article_elem.find('.//PMID')
        pmid = pmid_elem.text if pmid_elem is not None else ""

        title_elem = article_elem.find('.//ArticleTitle')
        title = title_elem.text.strip() if title_elem is not None and title_elem.text else "N/A"

        pub_date = parse_pub_date(article_elem)
        authors, affiliations = extract_authors_and_affiliations(article_elem)
        email = extract_corresponding_author_email(article_elem) or "N/A"

        articles.append({
            "PubmedID": pmid,
            "Title": title,
            "Publication Date": pub_date,
            "Non-academic Author(s)": ", ".join(authors) if authors else "N/A",
            "Company Affiliation(s)": ", ".join(affiliations) if affiliations else "N/A",
            "Corresponding Author Email": email
        })
    return articles

def parse_pub_date(article_elem: ET.Element) -> str:
    """Extract publication date from XML element and format as YYYY-MM-DD."""
    pub_date_elem = article_elem.find('.//PubMedPubDate[@PubStatus="pubmed"]')
    if pub_date_elem is not None:
        year = pub_date_elem.findtext('Year', '').strip()
        month = pub_date_elem.findtext('Month', '').strip().zfill(2)
        day = pub_date_elem.findtext('Day', '').strip().zfill(2)
        return f"{year}-{month}-{day}" if year else "N/A"
    return "N/A"

def extract_authors_and_affiliations(article: ET.Element) -> tuple:
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
        affiliation = aff_element.text.lower() if aff_element is not None and aff_element.text else ""
        
        # Exclude authors if their affiliation contains academic keywords
        if affiliation and not any(keyword in affiliation for keyword in ACADEMIC_KEYWORDS):
            authors.append(full_name)
            affiliations.add(affiliation)
            
    return authors, affiliations

def extract_corresponding_author_email(article: ET.Element) -> str:
    """Extract corresponding author's email if available."""
    for affiliation in article.findall(".//Affiliation"):
        if affiliation.text:
            email_match = re.search(r"([\w\.-]+@[\w\.-]+\.\w+)", affiliation.text)
            if email_match:
                return email_match.group(1)
    return ""

def save_to_csv(data: list, filename: str = "results.csv"):
    """Save results to a CSV file using UTF-8 encoding."""
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False, encoding="utf-8")
