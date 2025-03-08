import argparse
from .fetcher import fetch_paper_ids, fetch_paper_details, save_to_csv


def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers and extract company-affiliated authors.")
    parser.add_argument("query", type=str, help="Search query for PubMed")
    parser.add_argument("-f", "--file", type=str, help="Output CSV filename", default="results.csv")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")

    args = parser.parse_args()

    paper_ids = fetch_paper_ids(args.query)

    if args.debug:
        print(f"Found {len(paper_ids)} paper IDs: {paper_ids}")

    papers = fetch_paper_details(paper_ids)

    if args.debug:
        for paper in papers:
            print(paper)

    save_to_csv(papers, args.file)
    print(f"Results saved to {args.file}")


if __name__ == "__main__":
    main()
