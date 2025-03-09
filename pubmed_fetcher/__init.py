# pubmed_fetcher/__init__.py
from .client import search_pubmed, fetch_articles
from .processor import process_article
from .models import ArticleData, Author

__version__ = "0.1.0"

__all__ = [
    "search_pubmed",
    "fetch_articles",
    "process_article",
    "ArticleData",
    "Author",
]