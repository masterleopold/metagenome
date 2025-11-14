"""
MAFF (Ministry of Agriculture, Forestry and Fisheries) Web Scraper

Daily monitoring of livestock disease surveillance reports from MAFF website.
Detects new reports and scans for target virus keywords.

Target URL: https://www.maff.go.jp/j/syouan/douei/kansi_densen/
"""

import os
import re
import json
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Set, Any
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from bs4 import BeautifulSoup
import boto3
from botocore.exceptions import ClientError
import hashlib


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MAFFScraper:
    """
    Scraper for MAFF livestock disease surveillance reports

    Monitors daily updates, downloads new reports (PDF/Excel),
    and scans for target virus keywords.
    """

    # MAFF URLs
    BASE_URL = "https://www.maff.go.jp"
    TARGET_URL = "https://www.maff.go.jp/j/syouan/douei/kansi_densen/"

    # Alternative URLs for comprehensive monitoring
    ADDITIONAL_URLS = [
        "https://www.maff.go.jp/j/syouan/douei/kansi_densen/kansi_densen.html",  # Main surveillance page
        "https://www.maff.go.jp/j/syouan/douei/kansi_densen/attach/pdf/index.html",  # PDF archive
    ]

    # Target virus keywords (Japanese)
    VIRUS_KEYWORDS = {
        "hantavirus": ["ハンタウイルス", "ハンタ", "hantavirus", "HANTAV"],
        "polyomavirus": ["ポリオーマウイルス", "ポリオーマウィルス", "ポリオーマ", "polyomavirus", "polyoma"],
        "spumavirus": ["スピューマウイルス", "スピューマウィルス", "スピューマ", "フォーミーウイルス", "spumavirus", "foamy virus"],
        "eeev": ["東部ウマ脳炎", "EEEV", "Eastern equine encephalitis", "ウマ脳炎ウイルス"]
    }

    # Additional disease keywords for context
    DISEASE_KEYWORDS = [
        "豚", "ブタ", "pig", "swine", "porcine",
        "感染", "発生", "outbreak", "infection",
        "疾病", "disease"
    ]

    def __init__(self, region: str = "ap-northeast-1"):
        """
        Initialize MAFF scraper

        Args:
            region: AWS region for S3/DynamoDB storage
        """
        self.region = region

        # Initialize AWS clients
        self.s3_client = boto3.client('s3', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)

        # Session for HTTP requests
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        })

        logger.info("Initialized MAFF scraper")

    def fetch_page(self, url: str) -> Optional[BeautifulSoup]:
        """
        Fetch and parse HTML page

        Args:
            url: Target URL

        Returns:
            BeautifulSoup object or None if failed
        """
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            response.encoding = response.apparent_encoding  # Handle Japanese encoding

            soup = BeautifulSoup(response.text, 'html.parser')
            logger.info(f"Successfully fetched: {url}")
            return soup

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch {url}: {e}")
            return None

    def extract_report_links(self, soup: BeautifulSoup, base_url: str) -> List[Dict[str, str]]:
        """
        Extract PDF and Excel report links from page

        Args:
            soup: BeautifulSoup object
            base_url: Base URL for resolving relative links

        Returns:
            List of report metadata (url, title, file_type, date)
        """
        reports = []

        # Find all links to PDF and Excel files
        for link in soup.find_all('a', href=True):
            href = link['href']

            # Check if link is PDF or Excel
            if not (href.endswith('.pdf') or href.endswith('.xls') or href.endswith('.xlsx')):
                continue

            # Resolve relative URLs
            full_url = urljoin(base_url, href)

            # Extract title (link text or parent text)
            title = link.get_text(strip=True)
            if not title:
                parent = link.find_parent()
                title = parent.get_text(strip=True) if parent else "Untitled"

            # Determine file type
            file_type = 'pdf' if href.endswith('.pdf') else 'excel'

            # Try to extract date from title or URL
            date_match = re.search(r'(\d{4})年?(\d{1,2})月?(\d{1,2})?日?|(\d{4})-(\d{2})-(\d{2})|(\d{4})(\d{2})(\d{2})', title + href)
            extracted_date = None
            if date_match:
                groups = [g for g in date_match.groups() if g]
                if len(groups) >= 2:
                    extracted_date = '-'.join(groups[:3]) if len(groups) >= 3 else '-'.join(groups[:2])

            reports.append({
                'url': full_url,
                'title': title,
                'file_type': file_type,
                'date': extracted_date,
                'discovered_at': datetime.now().isoformat()
            })

        logger.info(f"Found {len(reports)} report links")
        return reports

    def download_report(self, url: str, save_dir: Path) -> Optional[Path]:
        """
        Download report file

        Args:
            url: Report URL
            save_dir: Directory to save file

        Returns:
            Path to downloaded file or None if failed
        """
        try:
            # Create filename from URL
            filename = os.path.basename(urlparse(url).path)
            if not filename:
                filename = hashlib.md5(url.encode()).hexdigest() + '.pdf'

            filepath = save_dir / filename

            # Download file
            response = self.session.get(url, timeout=60, stream=True)
            response.raise_for_status()

            # Save to disk
            filepath.parent.mkdir(parents=True, exist_ok=True)
            with open(filepath, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(f"Downloaded: {filename} ({filepath.stat().st_size} bytes)")
            return filepath

        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            return None

    def scan_text_for_keywords(self, text: str) -> Dict[str, List[str]]:
        """
        Scan text for virus and disease keywords

        Args:
            text: Text to scan

        Returns:
            Dictionary of found keywords by category
        """
        found = {
            'virus_keywords': [],
            'disease_keywords': []
        }

        # Check virus keywords
        for virus_type, keywords in self.VIRUS_KEYWORDS.items():
            for keyword in keywords:
                if keyword in text:
                    found['virus_keywords'].append(f"{virus_type}:{keyword}")
                    logger.warning(f"Found target virus keyword: {keyword} ({virus_type})")

        # Check disease keywords for context
        for keyword in self.DISEASE_KEYWORDS:
            if keyword in text:
                found['disease_keywords'].append(keyword)

        return found

    def check_if_new_report(self, report_url: str, table_name: str = 'surveillance-external-updates') -> bool:
        """
        Check if report has been processed before

        Args:
            report_url: Report URL
            table_name: DynamoDB table name

        Returns:
            True if report is new, False if already processed
        """
        try:
            table = self.dynamodb.Table(table_name)

            # Create hash of URL as identifier
            url_hash = hashlib.md5(report_url.encode()).hexdigest()

            response = table.get_item(
                Key={
                    'source#date': f"maff#{url_hash}",
                    'update_id': 'report_check'
                }
            )

            return 'Item' not in response  # New if not found

        except ClientError as e:
            logger.warning(f"DynamoDB check failed, assuming new report: {e}")
            return True

    def mark_report_as_processed(self, report: Dict[str, str], table_name: str = 'surveillance-external-updates') -> None:
        """
        Mark report as processed in DynamoDB

        Args:
            report: Report metadata
            table_name: DynamoDB table name
        """
        try:
            table = self.dynamodb.Table(table_name)

            url_hash = hashlib.md5(report['url'].encode()).hexdigest()

            table.put_item(Item={
                'source#date': f"maff#{url_hash}",
                'update_id': 'report_check',
                'source': 'maff',
                'report_url': report['url'],
                'report_title': report['title'],
                'processed_at': datetime.now().isoformat(),
                'timestamp': int(datetime.now().timestamp())
            })

        except ClientError as e:
            logger.error(f"Failed to mark report as processed: {e}")

    def daily_check(self, bucket_name: str) -> Dict[str, Any]:
        """
        Daily execution: Check MAFF website for new reports

        Args:
            bucket_name: S3 bucket for storing reports

        Returns:
            Summary of daily check results
        """
        logger.info("Starting daily MAFF surveillance check")

        today = datetime.now()
        date_str = today.strftime("%Y-%m-%d")

        results = {
            'date': date_str,
            'source': 'maff',
            'pages_checked': 0,
            'reports_found': 0,
            'new_reports': 0,
            'virus_keywords_found': [],
            's3_paths': [],
            'errors': []
        }

        all_reports = []

        # Check main URL and additional URLs
        urls_to_check = [self.TARGET_URL] + self.ADDITIONAL_URLS

        for url in urls_to_check:
            try:
                logger.info(f"Checking URL: {url}")
                soup = self.fetch_page(url)

                if not soup:
                    results['errors'].append(f"Failed to fetch {url}")
                    continue

                results['pages_checked'] += 1

                # Extract report links
                reports = self.extract_report_links(soup, url)
                all_reports.extend(reports)
                results['reports_found'] += len(reports)

                # Scan page text for keywords
                page_text = soup.get_text()
                keywords_found = self.scan_text_for_keywords(page_text)

                if keywords_found['virus_keywords']:
                    results['virus_keywords_found'].extend(keywords_found['virus_keywords'])
                    logger.warning(f"Found virus keywords on page {url}: {keywords_found['virus_keywords']}")

            except Exception as e:
                error_msg = f"Error checking {url}: {str(e)}"
                logger.error(error_msg)
                results['errors'].append(error_msg)

        # Process new reports
        temp_dir = Path(f"/tmp/maff_reports_{date_str}")
        temp_dir.mkdir(parents=True, exist_ok=True)

        for report in all_reports:
            try:
                # Check if already processed
                if not self.check_if_new_report(report['url']):
                    logger.debug(f"Skipping already processed report: {report['title']}")
                    continue

                results['new_reports'] += 1
                logger.info(f"Processing new report: {report['title']}")

                # Download report
                filepath = self.download_report(report['url'], temp_dir)

                if filepath:
                    # Upload to S3
                    s3_key = f"external/maff/{today.year}/{today.month:02d}/{today.day:02d}/{filepath.name}"
                    self._upload_to_s3(bucket_name, s3_key, filepath)
                    results['s3_paths'].append(f"s3://{bucket_name}/{s3_key}")

                    # Mark as processed
                    self.mark_report_as_processed(report)

                    # TODO: Extract text from PDF/Excel and scan for keywords
                    # This requires additional libraries (PyPDF2, openpyxl, pdfplumber)

            except Exception as e:
                error_msg = f"Error processing report {report['title']}: {str(e)}"
                logger.error(error_msg)
                results['errors'].append(error_msg)

        # Save summary to DynamoDB
        self._save_summary_to_dynamodb(results)

        logger.info(f"Daily MAFF check completed: {results['pages_checked']} pages, "
                   f"{results['new_reports']} new reports, "
                   f"{len(results['virus_keywords_found'])} virus keywords found")

        return results

    def _upload_to_s3(self, bucket_name: str, key: str, filepath: Path) -> None:
        """
        Upload file to S3

        Args:
            bucket_name: S3 bucket name
            key: S3 object key
            filepath: Local file path
        """
        try:
            with open(filepath, 'rb') as f:
                self.s3_client.put_object(
                    Bucket=bucket_name,
                    Key=key,
                    Body=f,
                    ContentType=self._get_content_type(filepath)
                )
            logger.info(f"Uploaded to S3: s3://{bucket_name}/{key}")

        except ClientError as e:
            logger.error(f"Failed to upload to S3: {e}")
            raise

    def _get_content_type(self, filepath: Path) -> str:
        """Get content type based on file extension"""
        ext = filepath.suffix.lower()
        content_types = {
            '.pdf': 'application/pdf',
            '.xls': 'application/vnd.ms-excel',
            '.xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        }
        return content_types.get(ext, 'application/octet-stream')

    def _save_summary_to_dynamodb(self, results: Dict[str, Any]) -> None:
        """
        Save daily check summary to DynamoDB

        Args:
            results: Daily check summary
        """
        try:
            table = self.dynamodb.Table('surveillance-external-updates')

            table.put_item(Item={
                'source#date': f"maff#{results['date']}",
                'update_id': datetime.now().isoformat(),
                'source': results['source'],
                'update_date': results['date'],
                'pages_checked': results['pages_checked'],
                'reports_found': results['reports_found'],
                'new_reports': results['new_reports'],
                'virus_keywords_found': results.get('virus_keywords_found', []),
                's3_paths': results.get('s3_paths', []),
                'errors': results.get('errors', []),
                'timestamp': int(datetime.now().timestamp())
            })

            logger.info(f"Saved MAFF update to DynamoDB: maff#{results['date']}")

        except ClientError as e:
            logger.error(f"Failed to save to DynamoDB: {e}")


if __name__ == "__main__":
    # Example usage
    scraper = MAFFScraper()

    # Test: Daily check (requires S3 bucket)
    # results = scraper.daily_check(bucket_name="surveillance-data")
    # print(json.dumps(results, indent=2, ensure_ascii=False))

    # Test: Fetch page
    soup = scraper.fetch_page(scraper.TARGET_URL)
    if soup:
        reports = scraper.extract_report_links(soup, scraper.TARGET_URL)
        print(f"Found {len(reports)} reports")
        for report in reports[:5]:  # Show first 5
            print(f"  - {report['title'][:50]}... ({report['file_type']})")
