import argparse
import logging
import sys
from pubmed_fetcher.fetcher import search_pubmed, fetch_article_details, save_to_csv

def setup_logging(debug: bool):
    """Configures logging."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers with industry authors.")
    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-f", "--file", help="Output CSV file name", default="results.csv")
    parser.add_argument("-m", "--max-results", type=int, default=100, help="Maximum number of results to fetch")
    parser.add_argument("-y", "--year", type=str, help="Year or range of publication (e.g., 2023 or 2020-2023)", required=True)

    args = parser.parse_args()
    setup_logging(args.debug)

    try:
        logging.info(f"Searching PubMed for query: {args.query}")

        # 🔹 Add Year Parsing Here (NEW)
        if "-" in args.year:
            start_year, end_year = map(int, args.year.split("-"))
        else:
            start_year = int(args.year)
            end_year = None  # Single year search

        # 🔹 Call search_pubmed with the parsed year values
        pmids = search_pubmed(args.query, start_year, end_year, args.max_results)

        if not pmids:
            logging.warning("No articles found matching the criteria.")
            sys.exit(0)

        articles = fetch_article_details(pmids)

        if args.debug:
            for article in articles:
                logging.debug(article)

        save_to_csv(articles, args.file)
        logging.info(f"Results saved to {args.file}")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)
