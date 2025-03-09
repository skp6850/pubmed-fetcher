# **PubMed Papers Fetcher**

## **📌 Overview**  
This **Python-based command-line tool** retrieves **research papers from PubMed** based on a user-specified query. It **filters out academic-only authors** and includes only papers where at least one author is affiliated with a **pharmaceutical or biotech company**. The results are exported as a **CSV file** or displayed in the console.  

---

## **🚀 Features**  
✔ **Fetches research papers using the PubMed API** (E-utilities: `esearch` and `efetch`).  
✔ **Filters non-academic authors** based on affiliation heuristics.  
✔ **Extracts key information** including PubMed ID, title, authors, affiliations, and emails.  
✔ **Outputs structured CSV data** for easy analysis.  
✔ **Command-line options for flexibility and control.**  

---

## **📌 How Filtering Works**  
🔹 **Non-Academic Authors** → Authors **without** keywords like:  
   - `University`, `College`, `Institute`, `Hospital`, `School`, `Lab`.  

🔹 **Company Affiliations** → Authors **with** keywords like:  
   - `pharma`, `pharmaceutical`, `biotech`, `inc`, `corporation`, `company`, `llc`, `labs`, `r&d`.  

🔹 **If at least one author has a company affiliation, the paper is included.**  

---

## **📌 Installation & Setup**  

### **1️⃣ Clone the Repository**  
```sh
git clone https://github.com/yourusername/pubmed_papers.git
cd pubmed_papers
```

### **2️⃣ Install Dependencies (Poetry)**  
Ensure you have **[Poetry](https://python-poetry.org/)** installed. Then run:  
```sh
poetry install
```
💪 This will install all required dependencies.  

---

## **📌 Running the Program**  
Run the command-line tool using:  
```sh
poetry run get-papers-list "COVID-19" -y 2023 -m 10 -f results.csv
```
OR, if installed via **pip**:
```sh
get-papers-list "COVID-19" -y 2023 -m 10 -f results.csv
```
💪 This fetches **10 research papers** on COVID-19 from **2023** and saves them to `results.csv`.  

---

## **📌 Command-Line Options**  
| **Option** | **Description** | **Example** |
|------------|----------------|-------------|
| `query` | **Search term** for PubMed | `"cancer treatment"` |
| `-y`, `--year` | **Filter by publication year or range** | `-y 2020-2023` |
| `-m`, `--max-results` | **Number of results to fetch** (default: `100`) | `-m 50` |
| `-f`, `--file` | **Filename for CSV output** (default: print to console) | `-f output.csv` |
| `-d`, `--debug` | **Enable debugging logs** | `-d` |

---

## **📌 Output Format (CSV)**  
Each row in the **CSV file** contains:  

| **Column** | **Description** |
|------------|----------------|
| `PubmedID` | Unique identifier for the paper |
| `Title` | Title of the research paper |
| `Publication Date` | Date of publication (YYYY-MM-DD) |
| `Non-academic Author(s)` | Names of industry-affiliated authors |
| `Company Affiliation(s)` | Names of pharmaceutical/biotech companies |
| `Corresponding Author Email` | Email of the corresponding author (if available) |

---

## **📌 Code Organization**  
```
pubmed_papers/
│— pubmed_fetcher/         # ✅ Main package
│   ├── __init__.py         # ✅ Marks this as a package
│   ├── cli.py              # ✅ Handles command-line interactions
│   ├── fetcher.py          # ✅ Fetches articles from PubMed
│   ├── client.py           # ✅ Handles API requests and XML parsing
│   ├── processor.py        # ✅ Filters and processes extracted data
│   ├── models.py           # ✅ Defines data models for articles & authors
│   ├── exceptions.py       # ✅ Custom error handling
│— tests/                  # ✅ Unit tests
│   ├── test_fetcher.py     # ✅ Tests for fetcher functionality
│— pyproject.toml          # ✅ Poetry config file
│— poetry.lock             # ✅ Poetry lock file
│— README.md               # ✅ Documentation
```

---

## **📌 Development & Testing**  
### **1️⃣ Run Unit Tests**  
To ensure correctness, run:  
```sh
pytest tests/
```

### **2️⃣ Debugging Issues**  
If you encounter problems, run with `--debug`:  
```sh
get-papers-list "diabetes" -y 2023 -m 5 -f diabetes.csv -d
```

---

## **📌 Publishing the Package**  
### **1️⃣ Build the Package**  
```sh
poetry build
```
### **2️⃣ Publish to TestPyPI (For Testing)**  
```sh
poetry publish --repository testpypi
```
### **3️⃣ Publish to PyPI (For Production)**  
```sh
poetry publish
```

---

## **📌 License**  
This project is licensed under the **MIT License**.  

---
Copyright (c) 2018 The Python Packaging Authority

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
