"""
Academic Literature Monitor for Virus Surveillance

Daily monitoring of academic publications (PubMed, J-STAGE) for
target virus research in Japanese pig populations.

Data Sources:
- PubMed (NCBI): https://pubmed.ncbi.nlm.nih.gov/
- J-STAGE (Japanese): https://www.jstage.jst.go.jp/
"""

import os
import json
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
from pathlib import Path
import requests
from xml.etree import ElementTree as ET
import boto3
from botocore.exceptions import ClientError


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AcademicMonitor:
    """
    Monitor academic publications for target virus research

    Searches PubMed and J-STAGE for recent publications related to
    hantavirus, polyomavirus, spumavirus, and EEEV in pigs.
    """

    # PubMed E-utilities API
    PUBMED_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    PUBMED_SEARCH_URL = f"{PUBMED_BASE_URL}/esearch.fcgi"
    PUBMED_FETCH_URL = f"{PUBMED_BASE_URL}/efetch.fcgi"
    PUBMED_SUMMARY_URL = f"{PUBMED_BASE_URL}/esummary.fcgi"

    # J-STAGE API (if available - may need to use web scraping)
    JSTAGE_SEARCH_URL = "https://www.jstage.jst.go.jp/AF06S010SryTopHyj"

    # Search queries for each virus
    SEARCH_QUERIES = {
        "hantavirus": {
            "pubmed": '("hantavirus"[All Fields] OR "hantaan"[All Fields] OR "seoul virus"[All Fields]) AND ("swine"[All Fields] OR "pig"[All Fields] OR "porcine"[All Fields]) AND "Japan"[All Fields]',
            "jstage_jp": "ハンタウイルス AND (豚 OR ブタ)",
            "jstage_en": "hantavirus AND (swine OR pig)"
        },
        "polyomavirus": {
            "pubmed": '("polyomavirus"[All Fields] OR "polyoma virus"[All Fields]) AND ("swine"[All Fields] OR "pig"[All Fields] OR "porcine"[All Fields] OR "Sus scrofa"[All Fields])',
            "jstage_jp": "ポリオーマウイルス AND (豚 OR ブタ)",
            "jstage_en": "polyomavirus AND (swine OR pig)"
        },
        "spumavirus": {
            "pubmed": '("spumavirus"[All Fields] OR "foamy virus"[All Fields] OR "porcine foamy virus"[All Fields]) AND ("swine"[All Fields] OR "pig"[All Fields] OR "porcine"[All Fields])',
            "jstage_jp": "スピューマウイルス OR フォーミーウイルス AND (豚 OR ブタ)",
            "jstage_en": "spumavirus OR foamy virus AND (swine OR pig)"
        },
        "eeev": {
            "pubmed": '("eastern equine encephalitis"[All Fields] OR "EEEV"[All Fields] OR "EEE virus"[All Fields]) AND ("swine"[All Fields] OR "pig"[All Fields] OR "porcine"[All Fields])',
            "jstage_jp": "東部ウマ脳炎 OR EEEV",
            "jstage_en": "eastern equine encephalitis OR EEEV"
        }
    }

    def __init__(self, email: str = "surveillance@example.com", region: str = "ap-northeast-1"):
        """
        Initialize academic monitor

        Args:
            email: Email for PubMed API (required by NCBI)
            region: AWS region for S3/DynamoDB storage
        """
        self.email = email
        self.region = region

        # Initialize AWS clients
        self.s3_client = boto3.client('s3', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)

        # Session for API requests
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'MinION-Surveillance-System/1.0'
        })

        logger.info("Initialized Academic Monitor")

    def search_pubmed(
        self,
        query: str,
        days_back: int = 1,
        max_results: int = 100
    ) -> List[str]:
        """
        Search PubMed for recent publications

        Args:
            query: PubMed search query
            days_back: Number of days to search back
            max_results: Maximum number of results

        Returns:
            List of PubMed IDs (PMIDs)
        """
        try:
            # Calculate date range
            end_date = datetime.now()
            start_date = end_date - timedelta(days=days_back)

            date_filter = f"{start_date.strftime('%Y/%m/%d')}:{end_date.strftime('%Y/%m/%d')}[pdat]"
            full_query = f"{query} AND {date_filter}"

            params = {
                'db': 'pubmed',
                'term': full_query,
                'retmax': max_results,
                'retmode': 'json',
                'email': self.email,
                'tool': 'MinION-Surveillance'
            }

            response = self.session.get(self.PUBMED_SEARCH_URL, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()
            pmids = data.get('esearchresult', {}).get('idlist', [])

            logger.info(f"Found {len(pmids)} PubMed results for query: {query[:50]}...")
            return pmids

        except requests.exceptions.RequestException as e:
            logger.error(f"PubMed search failed: {e}")
            return []

    def fetch_pubmed_details(self, pmids: List[str]) -> List[Dict[str, Any]]:
        """
        Fetch detailed information for PubMed articles

        Args:
            pmids: List of PubMed IDs

        Returns:
            List of article details
        """
        if not pmids:
            return []

        try:
            params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'xml',
                'email': self.email,
                'tool': 'MinION-Surveillance'
            }

            response = self.session.get(self.PUBMED_FETCH_URL, params=params, timeout=60)
            response.raise_for_status()

            # Parse XML response
            root = ET.fromstring(response.content)
            articles = []

            for article_elem in root.findall('.//PubmedArticle'):
                article = self._parse_pubmed_article(article_elem)
                if article:
                    articles.append(article)

            logger.info(f"Fetched details for {len(articles)} articles")
            return articles

        except Exception as e:
            logger.error(f"Failed to fetch PubMed details: {e}")
            return []

    def _parse_pubmed_article(self, article_elem: ET.Element) -> Optional[Dict[str, Any]]:
        """
        Parse PubMed article XML element

        Args:
            article_elem: XML element for article

        Returns:
            Article metadata dictionary
        """
        try:
            pmid_elem = article_elem.find('.//PMID')
            title_elem = article_elem.find('.//ArticleTitle')
            abstract_elem = article_elem.find('.//AbstractText')
            journal_elem = article_elem.find('.//Journal/Title')
            pub_date_elem = article_elem.find('.//PubDate')

            # Extract publication date
            pub_date = None
            if pub_date_elem is not None:
                year = pub_date_elem.findtext('Year')
                month = pub_date_elem.findtext('Month', '01')
                day = pub_date_elem.findtext('Day', '01')
                if year:
                    pub_date = f"{year}-{month}-{day}"

            # Extract authors
            authors = []
            for author_elem in article_elem.findall('.//Author'):
                last_name = author_elem.findtext('LastName', '')
                fore_name = author_elem.findtext('ForeName', '')
                if last_name:
                    authors.append(f"{last_name} {fore_name}".strip())

            return {
                'pmid': pmid_elem.text if pmid_elem is not None else None,
                'title': title_elem.text if title_elem is not None else 'No title',
                'abstract': abstract_elem.text if abstract_elem is not None else '',
                'journal': journal_elem.text if journal_elem is not None else '',
                'pub_date': pub_date,
                'authors': authors,
                'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid_elem.text}/" if pmid_elem is not None else None,
                'source': 'pubmed',
                'discovered_at': datetime.now().isoformat()
            }

        except Exception as e:
            logger.error(f"Failed to parse PubMed article: {e}")
            return None

    def search_jstage(self, query_jp: str, query_en: str, days_back: int = 1) -> List[Dict[str, Any]]:
        """
        Search J-STAGE for recent publications via web scraping

        Args:
            query_jp: Japanese query
            query_en: English query
            days_back: Number of days to search back

        Returns:
            List of article metadata
        """
        logger.info(f"Searching J-STAGE for: {query_jp}")

        articles = []

        # Try Japanese query first
        try:
            jp_articles = self._scrape_jstage_search(query_jp, 'ja')
            articles.extend(jp_articles)
        except Exception as e:
            logger.error(f"J-STAGE Japanese search failed: {e}")

        # Then English query
        try:
            en_articles = self._scrape_jstage_search(query_en, 'en')
            articles.extend(en_articles)
        except Exception as e:
            logger.error(f"J-STAGE English search failed: {e}")

        # Filter by date if needed (J-STAGE doesn't have good date filtering)
        threshold_date = datetime.now() - timedelta(days=days_back)

        # Deduplicate by title
        seen_titles = set()
        unique_articles = []
        for article in articles:
            title = article.get('title', '').lower()
            if title and title not in seen_titles:
                seen_titles.add(title)
                unique_articles.append(article)

        logger.info(f"Found {len(unique_articles)} unique J-STAGE articles")
        return unique_articles

    def _scrape_jstage_search(self, query: str, lang: str = 'ja') -> List[Dict[str, Any]]:
        """
        Scrape J-STAGE search results

        Args:
            query: Search query
            lang: Language (ja or en)

        Returns:
            List of article metadata
        """
        from bs4 import BeautifulSoup
        from urllib.parse import quote

        articles = []

        # J-STAGE search URL
        # Using global search endpoint
        encoded_query = quote(query)
        search_url = f"https://www.jstage.jst.go.jp/result/global/-char/{lang}?globalSearchKey={encoded_query}&item=1"

        try:
            # Add delay to be respectful
            import time
            time.sleep(1)

            response = self.session.get(search_url, timeout=30)
            response.raise_for_status()
            response.encoding = 'utf-8'

            soup = BeautifulSoup(response.text, 'html.parser')

            # Find article entries
            # J-STAGE uses different structures, try multiple selectors
            article_items = soup.find_all('li', class_='searchlist-item')

            if not article_items:
                # Alternative: Try finding by article tag
                article_items = soup.find_all('div', class_='result-list-item')

            if not article_items:
                # Another alternative: look for title links
                title_links = soup.find_all('a', href=lambda x: x and '/article/' in x if x else False)
                if title_links:
                    logger.info(f"Found {len(title_links)} article links via alternative method")

            for item in article_items:
                try:
                    article_data = self._parse_jstage_article(item)
                    if article_data:
                        articles.append(article_data)
                except Exception as e:
                    logger.debug(f"Failed to parse J-STAGE article item: {e}")
                    continue

            # If no articles found with structured parsing, try simple link extraction
            if not articles:
                articles = self._fallback_jstage_parsing(soup)

        except Exception as e:
            logger.error(f"J-STAGE scraping failed for query '{query}': {e}")

        return articles

    def _parse_jstage_article(self, item) -> Optional[Dict[str, Any]]:
        """
        Parse J-STAGE article item from HTML

        Args:
            item: BeautifulSoup element

        Returns:
            Article metadata dict
        """
        from bs4 import BeautifulSoup

        article = {}

        # Extract title
        title_elem = item.find('h3') or item.find('h4') or item.find('a', class_='article-title')
        if title_elem:
            article['title'] = title_elem.get_text(strip=True)
        else:
            return None

        # Extract URL
        link_elem = item.find('a', href=lambda x: x and '/article/' in x if x else False)
        if link_elem:
            href = link_elem.get('href', '')
            if href.startswith('/'):
                href = f"https://www.jstage.jst.go.jp{href}"
            article['url'] = href

            # Try to extract DOI from URL
            if '/article/' in href:
                parts = href.split('/article/')
                if len(parts) > 1:
                    article['jstage_id'] = parts[1].split('/')[0]

        # Extract authors
        author_elem = item.find('p', class_='author') or item.find('div', class_='author')
        if author_elem:
            article['authors'] = author_elem.get_text(strip=True)

        # Extract journal name
        journal_elem = item.find('p', class_='journal') or item.find('span', class_='journal-title')
        if journal_elem:
            article['journal'] = journal_elem.get_text(strip=True)

        # Extract publication date
        date_elem = item.find('p', class_='date') or item.find('span', class_='pub-date')
        if date_elem:
            article['pub_date'] = date_elem.get_text(strip=True)

        # Extract abstract preview if available
        abstract_elem = item.find('p', class_='abstract')
        if abstract_elem:
            article['abstract_preview'] = abstract_elem.get_text(strip=True)

        # Add metadata
        article['source'] = 'jstage'
        article['discovered_at'] = datetime.now().isoformat()

        return article

    def _fallback_jstage_parsing(self, soup) -> List[Dict[str, Any]]:
        """
        Fallback parsing method for J-STAGE when structured parsing fails

        Args:
            soup: BeautifulSoup object

        Returns:
            List of article metadata
        """
        articles = []

        # Find all article links
        article_links = soup.find_all('a', href=lambda x: x and '/article/' in x if x else False)

        for link in article_links[:20]:  # Limit to first 20
            try:
                href = link.get('href', '')
                if href.startswith('/'):
                    href = f"https://www.jstage.jst.go.jp{href}"

                title = link.get_text(strip=True)
                if not title or len(title) < 10:
                    continue

                article = {
                    'title': title,
                    'url': href,
                    'source': 'jstage',
                    'discovered_at': datetime.now().isoformat(),
                    'parsing_method': 'fallback'
                }

                articles.append(article)

            except Exception as e:
                logger.debug(f"Fallback parsing error: {e}")
                continue

        logger.info(f"Fallback parsing found {len(articles)} articles")
        return articles

    def daily_pubmed_search(self, days_back: int = 1) -> Dict[str, Any]:
        """
        Daily PubMed search for all target viruses

        Args:
            days_back: Number of days to search back (default: 1 for daily)

        Returns:
            Summary of search results
        """
        logger.info(f"Starting daily PubMed search (last {days_back} days)")

        results = {
            'date': datetime.now().strftime("%Y-%m-%d"),
            'source': 'pubmed',
            'days_searched': days_back,
            'viruses_checked': 0,
            'total_articles': 0,
            'articles_by_virus': {},
            'new_publications': []
        }

        for virus_name, queries in self.SEARCH_QUERIES.items():
            try:
                logger.info(f"Searching PubMed for {virus_name}")

                # Search PubMed
                pmids = self.search_pubmed(
                    query=queries['pubmed'],
                    days_back=days_back,
                    max_results=100
                )

                results['viruses_checked'] += 1
                results['articles_by_virus'][virus_name] = len(pmids)
                results['total_articles'] += len(pmids)

                if pmids:
                    # Fetch article details
                    articles = self.fetch_pubmed_details(pmids)

                    for article in articles:
                        article['virus_type'] = virus_name
                        results['new_publications'].append(article)

                    logger.info(f"Found {len(articles)} new {virus_name} articles")

            except Exception as e:
                logger.error(f"Error searching for {virus_name}: {e}")

        return results

    def save_to_s3(self, bucket_name: str, results: Dict[str, Any]) -> List[str]:
        """
        Save search results to S3

        Args:
            bucket_name: S3 bucket name
            results: Search results

        Returns:
            List of S3 paths
        """
        s3_paths = []
        today = datetime.now()

        try:
            # Save summary
            summary_key = f"external/academic/{today.year}/{today.month:02d}/{today.day:02d}/summary_{results['date']}.json"
            self.s3_client.put_object(
                Bucket=bucket_name,
                Key=summary_key,
                Body=json.dumps(results, ensure_ascii=False, indent=2),
                ContentType="application/json"
            )
            s3_paths.append(f"s3://{bucket_name}/{summary_key}")

            # Save individual articles
            for article in results.get('new_publications', []):
                article_key = f"external/academic/{today.year}/{today.month:02d}/{today.day:02d}/articles/{article['virus_type']}_{article['pmid']}.json"
                self.s3_client.put_object(
                    Bucket=bucket_name,
                    Key=article_key,
                    Body=json.dumps(article, ensure_ascii=False, indent=2),
                    ContentType="application/json"
                )

            logger.info(f"Saved {len(results.get('new_publications', []))} articles to S3")

        except ClientError as e:
            logger.error(f"Failed to save to S3: {e}")

        return s3_paths

    def save_to_dynamodb(self, results: Dict[str, Any]) -> None:
        """
        Save search results summary to DynamoDB

        Args:
            results: Search results
        """
        try:
            table = self.dynamodb.Table('surveillance-external-updates')

            table.put_item(Item={
                'source#date': f"academic#{results['date']}",
                'update_id': datetime.now().isoformat(),
                'source': results['source'],
                'update_date': results['date'],
                'days_searched': results['days_searched'],
                'viruses_checked': results['viruses_checked'],
                'total_articles': results['total_articles'],
                'articles_by_virus': results.get('articles_by_virus', {}),
                'new_items_count': results['total_articles'],
                'jstage_articles': results.get('jstage_articles', 0),
                'pubmed_articles': results.get('pubmed_articles', 0),
                'timestamp': int(datetime.now().timestamp())
            })

            logger.info(f"Saved academic update to DynamoDB: academic#{results['date']}")

        except ClientError as e:
            logger.error(f"Failed to save to DynamoDB: {e}")

    def daily_jstage_search(self, days_back: int = 1) -> Dict[str, Any]:
        """
        Daily J-STAGE search for all target viruses

        Args:
            days_back: Number of days to search back (default: 1 for daily)

        Returns:
            Summary of search results
        """
        logger.info(f"Starting daily J-STAGE search (last {days_back} days)")

        results = {
            'date': datetime.now().strftime("%Y-%m-%d"),
            'source': 'jstage',
            'days_searched': days_back,
            'viruses_checked': 0,
            'total_articles': 0,
            'articles_by_virus': {},
            'new_publications': []
        }

        for virus_name, queries in self.SEARCH_QUERIES.items():
            try:
                logger.info(f"Searching J-STAGE for {virus_name}")

                # Search J-STAGE
                articles = self.search_jstage(
                    query_jp=queries['jstage_jp'],
                    query_en=queries['jstage_en'],
                    days_back=days_back
                )

                results['viruses_checked'] += 1
                results['articles_by_virus'][virus_name] = len(articles)
                results['total_articles'] += len(articles)

                # Add virus type to articles
                for article in articles:
                    article['virus_type'] = virus_name
                    results['new_publications'].append(article)

                logger.info(f"Found {len(articles)} J-STAGE articles for {virus_name}")

            except Exception as e:
                logger.error(f"Error searching J-STAGE for {virus_name}: {e}")

        return results

    def daily_monitor(self, bucket_name: str, days_back: int = 1) -> Dict[str, Any]:
        """
        Daily execution: Monitor academic publications from PubMed and J-STAGE

        Args:
            bucket_name: S3 bucket for storing results
            days_back: Number of days to search back

        Returns:
            Summary of daily monitoring
        """
        logger.info("Starting daily academic monitoring (PubMed + J-STAGE)")

        # Search PubMed
        pubmed_results = self.daily_pubmed_search(days_back=days_back)

        # Search J-STAGE
        jstage_results = self.daily_jstage_search(days_back=days_back)

        # Combine results
        combined_results = {
            'date': pubmed_results['date'],
            'source': 'academic',
            'days_searched': days_back,
            'viruses_checked': pubmed_results['viruses_checked'],
            'total_articles': pubmed_results['total_articles'] + jstage_results['total_articles'],
            'pubmed_articles': pubmed_results['total_articles'],
            'jstage_articles': jstage_results['total_articles'],
            'articles_by_virus': {},
            'new_publications': []
        }

        # Merge articles by virus
        for virus in self.SEARCH_QUERIES.keys():
            pubmed_count = pubmed_results['articles_by_virus'].get(virus, 0)
            jstage_count = jstage_results['articles_by_virus'].get(virus, 0)
            combined_results['articles_by_virus'][virus] = pubmed_count + jstage_count

        # Combine publications
        combined_results['new_publications'].extend(pubmed_results['new_publications'])
        combined_results['new_publications'].extend(jstage_results['new_publications'])

        # Save combined results
        s3_paths = self.save_to_s3(bucket_name, combined_results)
        combined_results['s3_paths'] = s3_paths

        self.save_to_dynamodb(combined_results)

        logger.info(f"Daily academic monitoring completed: {combined_results['total_articles']} articles "
                   f"(PubMed: {combined_results['pubmed_articles']}, J-STAGE: {combined_results['jstage_articles']})")

        return combined_results


if __name__ == "__main__":
    # Example usage
    monitor = AcademicMonitor(email="surveillance@example.com")

    # Test: Daily PubMed search (last 7 days for testing)
    print("=== Testing PubMed Search ===")
    pubmed_results = monitor.daily_pubmed_search(days_back=7)
    print(f"PubMed: {pubmed_results['total_articles']} articles found")
    print(json.dumps(pubmed_results, indent=2, ensure_ascii=False))

    # Test: J-STAGE search
    print("\n=== Testing J-STAGE Search ===")
    jstage_results = monitor.daily_jstage_search(days_back=7)
    print(f"J-STAGE: {jstage_results['total_articles']} articles found")
    print(json.dumps(jstage_results, indent=2, ensure_ascii=False))

    # Test: Full daily monitor (requires S3 bucket)
    # print("\n=== Testing Full Daily Monitor ===")
    # results = monitor.daily_monitor(bucket_name="surveillance-data", days_back=1)
    # print(json.dumps(results, indent=2, ensure_ascii=False))
