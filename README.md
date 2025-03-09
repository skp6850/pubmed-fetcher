# PubMed Papers Fetcher

This Python program fetches research papers from PubMed based on a user-specified query. It filters for papers that have at least one author affiliated with a pharmaceutical or biotech company and outputs the results as a CSV file (or prints to the console).

## Features

- **Source of Papers:** Uses the PubMed API (E-utilities: esearch and efetch).
- **Filtering:** Identifies non-academic authors using simple heuristics:
  - Non-academic affiliations are determined by the absence of keywords such as `University`, `College`, `Institute`, `Hospital`, or `School`.
  - Company affiliations are identified by keywords like `pharma`, `pharmaceutical`, `biotech`, `inc`, etc.
- **Output:** CSV with columns:
  - PubmedID
  - Title
  - Publication Date
  - Non-academic Author(s)
  - Company Affiliation(s)
  - Corresponding Author Email
- **Command-line Options:**
  - `-h` / `--help`: Display usage instructions.
  - `-d` / `--debug`: Print debug information.
  - `-f` / `--file`: Specify the filename to save the CSV results.
- **Code Organization:**
  - The logic is split into a module (`pubmed_fetcher.py`) and a CLI entry point (`main.py`).
- **Environment & Dependencies:**
  - Managed with [Poetry](https://python-poetry.org/).
  - To install dependencies, run: `poetry install`
- **Executable Command:**
  - The project is configured (via `pyproject.toml`) to expose an executable command `get-papers-list` after installation.
- **Bonus:**
  - The module is separated from the CLI.
  - (Optionally, you can publish the module on TestPyPI.)

## Setup and Execution

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/pubmed_papers.git
   cd pubmed_papers









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
