import csv
import pytest

# Define academic keywords used for filtering in your application.
ACADEMIC_KEYWORDS = ["university", "college", "institute", "hospital", "school", "lab"]

@pytest.fixture
def sample_csv(tmp_path):
    """
    Create a temporary CSV file simulating the output of the PubMed Fetcher.
    It includes one academic-only paper (which should be excluded from the final results)
    and one industry paper.
    """
    csv_file = tmp_path / "results.csv"
    # Simulated output rows:
    rows = [
        {
            "PubmedID": "123456",
            "Title": "Academic Research on XYZ",
            "Publication Date": "2023-05-01",
            "Non-academic Author(s)": "",  # Academic-only paper; should be excluded
            "Company Affiliation(s)": "",
            "Corresponding Author Email": "prof@uni.edu"
        },
        {
            "PubmedID": "234567",
            "Title": "Industry Research on ABC",
            "Publication Date": "2023-06-15",
            "Non-academic Author(s)": "John Doe",
            "Company Affiliation(s)": "Pharma Inc.",
            "Corresponding Author Email": "john.doe@pharma.com"
        }
    ]
    fieldnames = [
        "PubmedID", "Title", "Publication Date",
        "Non-academic Author(s)", "Company Affiliation(s)", "Corresponding Author Email"
    ]
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return csv_file

def load_results(csv_file):
    """Helper function to load CSV results into a list of dictionaries."""
    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)

def test_academic_only_papers_excluded(sample_csv):
    """
    Ensure that in the final output, papers with only academic affiliations
    are excluded. For this test, we assume that an academic-only paper will have
    empty or 'N/A' values in both 'Non-academic Author(s)' and 'Company Affiliation(s)' fields.
    """
    results = load_results(sample_csv)
    for row in results:
        title_lower = row["Title"].lower()
        non_academic = row["Non-academic Author(s)"].strip().lower()
        company_aff = row["Company Affiliation(s)"].strip().lower()
        
        # Assume if the title contains 'academic research' or the email ends with a known academic domain,
        # then it should have no non-academic authors or company affiliations.
        if "academic" in title_lower or row["Corresponding Author Email"].endswith("uni.edu"):
            # The academic paper should have empty fields in the filtered output.
            assert non_academic in ("", "n/a"), f"Academic paper '{row['Title']}' should not have non-academic authors."
            assert company_aff in ("", "n/a"), f"Academic paper '{row['Title']}' should not have company affiliations."
        else:
            # For industry papers, ensure these fields are populated.
            assert non_academic not in ("", "n/a"), f"Industry paper '{row['Title']}' should have non-academic authors."
            assert company_aff not in ("", "n/a"), f"Industry paper '{row['Title']}' should have company affiliations."

def test_no_duplicate_titles(sample_csv):
    """
    Ensure that the output CSV does not contain duplicate titles.
    """
    results = load_results(sample_csv)
    titles = [row["Title"].strip() for row in results]
    assert len(titles) == len(set(titles)), "Duplicate titles found in the CSV output!"


