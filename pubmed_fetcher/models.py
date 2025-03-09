# pubmed_fetcher/models.py
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class Author:
    name: str
    affiliations: List[str]

@dataclass
class ArticleData:
    pmid: str
    title: str
    pub_date: str
    authors: List[Author]
    email: Optional[str]