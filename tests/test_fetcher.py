import pytest
from pubmed_fetcher.fetcher import extract_company_authors, extract_corresponding_author_email

def test_extract_company_authors():
    authors = [
        {"name": "Dr. John", "affiliation": "Pfizer Inc."},
        {"name": "Dr. Jane", "affiliation": "Harvard University"}
    ]
    company_authors, company_affiliations = extract_company_authors(authors)
    assert "Dr. John" in company_authors
    assert "Pfizer Inc." in company_affiliations
    assert "Dr. Jane" not in company_authors

def test_extract_corresponding_author_email():
    email = extract_corresponding_author_email("Contact: author@company.com")
    assert email == "author@company.com"
