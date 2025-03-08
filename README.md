# PubMed Fetcher

## 📌 Overview
This tool fetches PubMed papers, filters out company-affiliated authors, and saves results in a CSV file.

## 📌 Installation
### 1️⃣ Install Poetry:
```sh
pip install poetry

### 2️⃣ Clone the repository:

git clone https://github.com/yourname/pubmed-fetcher.git
cd pubmed-fetcher
3️⃣ Install dependencies:

poetry install
📌 Usage
Run the following command:
poetry run get-papers-list "covid 19" -d -f results.csv

📌 How Non-Academic Authors Are Identified
Excludes affiliations with: "university", "college", "institute", "school"
Filters out company-affiliated authors.
Extracts corresponding author emails heuristically.

📌 External Dependencies
requests for API calls
pandas for data processing
pytest for testing

📌 License
MIT License.

🚀 **Your PubMed Fetcher is now fully documented, automated, and follows best practices!** 🚀