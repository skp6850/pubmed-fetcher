# pubmed_fetcher/processor.py
from typing import Dict
from .models import ArticleData

ACADEMIC_KEYWORDS = {
    'university', 'college', 'institute', 'school',
    'academy', 'hospital', 'lab', 'laboratory', 'foundation', 'clinic'
}
COMPANY_KEYWORDS = {
    'pharma', 'biotech', 'inc', 'ltd', 'corporation',
    'company', 'llc', 'co.', 'labs', 'research and development',
    'r&d', 'biopharma', 'biopharmaceutical'
}


def process_article(article: ArticleData) -> Dict:
    """Process an article to extract non-academic authors and company affiliations."""
    non_academic_authors = []
    company_affiliations = set()

    for author in article.authors:
        author_company_affiliations = set()
        is_non_academic = False

        for affiliation in author.affiliations:
            affiliation_lower = affiliation.lower()
            is_academic = any(keyword in affiliation_lower for keyword in ACADEMIC_KEYWORDS)
            if not is_academic:
                is_company = any(keyword in affiliation_lower for keyword in COMPANY_KEYWORDS)
                if is_company:
                    author_company_affiliations.add(affiliation)
                    is_non_academic = True

        if is_non_academic:
            non_academic_authors.append(author.name)
            company_affiliations.update(author_company_affiliations)

    return {
        'PubmedID': article.pmid,
        'Title': article.title,
        'Publication Date': article.pub_date,
        'Non-academic Author(s)': '; '.join(non_academic_authors),
        'Company Affiliation(s)': '; '.join(company_affiliations),
        'Corresponding Author Email': article.email or ''
    }