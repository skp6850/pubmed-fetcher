# pubmed_fetcher/client.py
import logging
import re
from typing import List, Optional
from xml.etree import ElementTree as ET
from Bio import Entrez
from .models import ArticleData, Author
from .exceptions import PubMedFetchError

Entrez.email = "skp68500@gmail.com"  # Set your email for PubMed API

logger = logging.getLogger(__name__)

ACADEMIC_KEYWORDS = {
    'university', 'college', 'institute', 'school',
    'academy', 'hospital', 'lab', 'laboratory', 'foundation', 'clinic'
}
COMPANY_KEYWORDS = {
    'pharma', 'biotech', 'inc', 'ltd', 'corporation',
    'company', 'llc', 'co.', 'labs', 'research and development',
    'r&d', 'biopharma', 'biopharmaceutical', 'vaccine', 'therapeutics'
}


def search_pubmed(query: str, max_results: int = 100) -> List[str]:
    """Search PubMed and return a list of PMIDs."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        logger.error(f"PubMed search failed: {e}")
        raise PubMedFetchError(f"Search failed: {e}") from e


def fetch_articles(pmids: List[str]) -> List[ArticleData]:
    """Fetch and parse article details for given PMIDs."""
    if not pmids:
        return []

    try:
        handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")
        data = handle.read()
        handle.close()
        return parse_articles_xml(data)
    except Exception as e:
        logger.error(f"Fetching articles failed: {e}")
        raise PubMedFetchError(f"Fetch failed: {e}") from e


def parse_articles_xml(xml_data: str) -> List[ArticleData]:
    """Parse PubMed XML data into ArticleData objects."""
    root = ET.fromstring(xml_data)
    articles = []

    for article_elem in root.findall('.//PubmedArticle'):
        pmid_elem = article_elem.find('.//PMID')
        pmid = pmid_elem.text if pmid_elem is not None else ""

        title_elem = article_elem.find('.//ArticleTitle')
        title = ''.join(title_elem.itertext()).strip() if title_elem is not None else ""

        pub_date = parse_pub_date(article_elem)
        authors = parse_authors(article_elem)
        email = parse_email(article_elem)

        articles.append(ArticleData(pmid, title, pub_date, authors, email))

    return articles


def parse_pub_date(article_elem: ET.Element) -> str:
    """Extract publication date from XML element and format as yyyy-mm-dd."""
    pub_date_elem = article_elem.find('.//PubMedPubDate[@PubStatus="pubmed"]')
    if pub_date_elem is not None:
        year = pub_date_elem.findtext('Year', '').strip()
        month = pub_date_elem.findtext('Month', '').strip().zfill(2)  # Pad with zero
        day = pub_date_elem.findtext('Day', '').strip().zfill(2)  # Pad with zero
        return f"{year}-{month}-{day}" if year else ""
    return ""


def parse_authors(article_elem: ET.Element) -> List[Author]:
    """Parse authors and their affiliations from XML."""
    authors = []
    author_list = article_elem.find('.//AuthorList')
    if author_list is None:
        return authors

    for author_elem in author_list.findall('Author'):
        last_name = author_elem.findtext('LastName', '').strip()
        fore_name = author_elem.findtext('ForeName', '').strip()
        if not last_name and not fore_name:
            continue

        name = f"{fore_name} {last_name}".strip()
        affiliations = []
        for affil_elem in author_elem.findall('.//AffiliationInfo/Affiliation'):
            affil = affil_elem.text.strip() if affil_elem.text else ""
            if affil:
                affiliations.append(affil)

        authors.append(Author(name, affiliations))

    return authors


def parse_email(article_elem: ET.Element) -> Optional[str]:
    """Extract the email from the CorrespondingAuthor's Email tag."""
    for author_elem in article_elem.findall('.//Author'):
        if author_elem.get("ValidYN", "N") == "Y":  # Corresponding author flag
            email_elem = author_elem.find('Email')
            if email_elem is not None and email_elem.text:
                return email_elem.text.strip()
    return None