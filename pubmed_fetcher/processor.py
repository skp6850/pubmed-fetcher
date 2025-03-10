# pubmed_fetcher/processor.py
from typing import Dict
from .models import ArticleData

COMPANY_KEYWORDS = {
    "pharma", "biotech", "inc", "ltd", "corporation",
    "company", "llc", "co.", "labs", "research and development",
    "r&d", "biopharma", "biopharmaceutical", "vaccine", "therapeutics"
}
ACADEMIC_KEYWORDS = {
    "university", "college", "institute", "school",
    "academy", "hospital", "lab", "foundation", "department", "research center"
}

def process_article(article):
    """Strictly exclude academic affiliations and extract non-academic authors."""
    non_academic_authors = []
    company_affiliations = set()

    for author in article.authors:
        affiliations = [aff.lower() for aff in author.affiliations]
        is_academic = any(keyword in " ".join(affiliations) for keyword in ACADEMIC_KEYWORDS)

        if not is_academic:
            is_company = any(keyword in " ".join(affiliations) for keyword in COMPANY_KEYWORDS)
            if is_company:
                non_academic_authors.append(author.name)
                company_affiliations.update(author.affiliations)

    return {
        "PubmedID": article.pmid,
        "Title": article.title,
        "Publication Date": article.pub_date,
        "Non-academic Author(s)": "; ".join(non_academic_authors),
        "Company Affiliation(s)": "; ".join(company_affiliations),
        "Corresponding Author Email": article.email or ""
    }