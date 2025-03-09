import pytest
from xml.etree.ElementTree import Element, SubElement
from pubmed_fetcher.fetcher import extract_authors_and_affiliations, extract_corresponding_author_email


@pytest.fixture
def sample_article():
    """Mock XML article for testing"""
    article = Element("PubmedArticle")
    medline = SubElement(article, "MedlineCitation")
    article_info = SubElement(medline, "Article")

    # Add title
    title = SubElement(article_info, "ArticleTitle")
    title.text = "Test Paper on AI in Healthcare"

    # Add PubDate
    pub_date = SubElement(article_info, "PubDate")
    year = SubElement(pub_date, "Year")
    year.text = "2023"

    # Add AuthorList
    author_list = SubElement(article_info, "AuthorList")

    author1 = SubElement(author_list, "Author")
    fore_name1 = SubElement(author1, "ForeName")
    fore_name1.text = "John"
    last_name1 = SubElement(author1, "LastName")
    last_name1.text = "Doe"
    aff1 = SubElement(author1, "Affiliation")
    aff1.text = "Pfizer Inc."

    author2 = SubElement(author_list, "Author")
    fore_name2 = SubElement(author2, "ForeName")
    fore_name2.text = "Jane"
    last_name2 = SubElement(author2, "LastName")
    aff2 = SubElement(author2, "Affiliation")
    aff2.text = "Harvard University"

    return article


def test_excludes_academic_authors(sample_article):
    """Ensure academic authors are correctly excluded from results."""
    authors, affiliations = extract_authors_and_affiliations(sample_article)

    # Check that industry author is included
    assert "John Doe" in authors
    assert "Pfizer Inc." in affiliations

    # Ensure academic author (Harvard University) is excluded
    assert "Jane Doe" not in authors
    assert "Harvard University" not in affiliations


def test_extract_corresponding_author_email():
    """Ensure email is extracted correctly from an XML element."""
    mock_article = Element("PubmedArticle")
    aff1 = SubElement(mock_article, "Affiliation")
    aff1.text = "Pfizer Inc., New York, USA. contact@pfizer.com"

    email = extract_corresponding_author_email(mock_article)
    assert email == "contact@pfizer.com"

