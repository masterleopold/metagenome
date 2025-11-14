"""
E-Stat API Client for Daily Livestock Statistics Monitoring

Fetches livestock disease statistics from Japanese government e-Stat portal.
Focuses on pig health indicators and regional disease trends.

API Documentation: https://www.e-stat.go.jp/api/
"""

import os
import json
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
from pathlib import Path
import requests
import boto3
from botocore.exceptions import ClientError


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class EStatClient:
    """
    Client for e-Stat API (Japanese Government Statistics Portal)

    Provides access to livestock statistics including pig health data,
    disease occurrence, and regional trends.
    """

    # E-Stat API Configuration
    APP_ID = "bae1f981a6d093a9676b03c8eea37324b8de421b"
    BASE_URL = "https://api.e-stat.go.jp/rest/3.0/app"

    # Statistical Table IDs (examples - need to be confirmed)
    STAT_IDS = {
        "livestock_disease": "0003103532",  # Livestock disease statistics
        "pig_population": "0003103533",     # Pig population by region
        "veterinary_inspection": "0003103534"  # Veterinary inspection results
    }

    # Target virus keywords for filtering
    VIRUS_KEYWORDS = [
        "ハンタウイルス", "ハンタ",
        "ポリオーマウイルス", "ポリオーマ", "ポリオーマウィルス",
        "スピューマウイルス", "スピューマ", "スピューマウィルス", "フォーミーウイルス",
        "東部ウマ脳炎", "EEEV", "ウマ脳炎"
    ]

    def __init__(self, app_id: Optional[str] = None, region: str = "ap-northeast-1"):
        """
        Initialize E-Stat API client

        Args:
            app_id: E-Stat application ID (defaults to class constant)
            region: AWS region for S3/DynamoDB storage
        """
        self.app_id = app_id or self.APP_ID
        self.region = region

        # Initialize AWS clients
        self.s3_client = boto3.client('s3', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)

        # Session for API requests
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'MinION-Surveillance-System/1.0'
        })

        logger.info(f"Initialized E-Stat client with APP_ID: {self.app_id[:8]}...")

    def get_stats_list(self, search_keyword: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Get list of available statistical tables

        Args:
            search_keyword: Filter tables by keyword (e.g., "豚", "家畜")

        Returns:
            List of statistical table metadata
        """
        endpoint = f"{self.BASE_URL}/getStatsList"
        params = {
            "appId": self.app_id,
            "lang": "J",  # Japanese
            "searchWord": search_keyword or "豚 疾病"  # Pig disease
        }

        try:
            response = self.session.get(endpoint, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Check for API errors
            if data.get("GET_STATS_LIST", {}).get("RESULT", {}).get("STATUS") != 0:
                error_msg = data.get("GET_STATS_LIST", {}).get("RESULT", {}).get("ERROR_MSG")
                logger.error(f"E-Stat API error: {error_msg}")
                return []

            stats_list = data.get("GET_STATS_LIST", {}).get("DATALIST_INF", {}).get("TABLE_INF", [])
            logger.info(f"Found {len(stats_list)} statistical tables")

            return stats_list

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch stats list: {e}")
            return []

    def get_stats_data(
        self,
        stats_id: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        cdArea: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get statistical data for a specific table

        Args:
            stats_id: Statistical table ID
            start_date: Start date (YYYY or YYYYMM)
            end_date: End date (YYYY or YYYYMM)
            cdArea: Area code (e.g., "00000" for national, "13000" for Tokyo)

        Returns:
            Statistical data in JSON format
        """
        endpoint = f"{self.BASE_URL}/getStatsData"
        params = {
            "appId": self.app_id,
            "lang": "J",
            "statsDataId": stats_id,
            "metaGetFlg": "Y",  # Include metadata
            "cntGetFlg": "N",   # Don't include count
            "sectionHeaderFlg": "1"  # Section headers
        }

        # Add optional filters
        if start_date:
            params["cdTimeFrom"] = start_date
        if end_date:
            params["cdTimeTo"] = end_date
        if cdArea:
            params["cdArea"] = cdArea

        try:
            response = self.session.get(endpoint, params=params, timeout=60)
            response.raise_for_status()

            data = response.json()

            # Check for API errors
            if data.get("GET_STATS_DATA", {}).get("RESULT", {}).get("STATUS") != 0:
                error_msg = data.get("GET_STATS_DATA", {}).get("RESULT", {}).get("ERROR_MSG")
                logger.error(f"E-Stat API error: {error_msg}")
                return {}

            logger.info(f"Successfully fetched data for stats_id: {stats_id}")
            return data

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch stats data: {e}")
            return {}

    def search_virus_keywords(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Search for virus-related keywords in statistical data

        Args:
            data: E-Stat API response data

        Returns:
            List of matching records with virus keywords
        """
        matches = []
        data_values = data.get("GET_STATS_DATA", {}).get("STATISTICAL_DATA", {}).get("DATA_INF", {}).get("VALUE", [])

        if not data_values:
            return matches

        for record in data_values:
            # Convert record to string for keyword search
            record_str = json.dumps(record, ensure_ascii=False)

            for keyword in self.VIRUS_KEYWORDS:
                if keyword in record_str:
                    matches.append({
                        "keyword": keyword,
                        "record": record,
                        "timestamp": datetime.now().isoformat()
                    })
                    logger.warning(f"Found virus keyword '{keyword}' in E-Stat data")
                    break

        return matches

    def daily_fetch_livestock_stats(self, bucket_name: str) -> Dict[str, Any]:
        """
        Daily execution: Fetch livestock statistics and check for virus keywords

        Args:
            bucket_name: S3 bucket for storing results

        Returns:
            Summary of daily fetch results
        """
        logger.info("Starting daily E-Stat livestock statistics fetch")

        today = datetime.now()
        date_str = today.strftime("%Y-%m-%d")

        results = {
            "date": date_str,
            "source": "estat",
            "stats_checked": 0,
            "virus_keywords_found": [],
            "s3_paths": [],
            "errors": []
        }

        # Fetch data for each relevant statistical table
        for stat_name, stat_id in self.STAT_IDS.items():
            try:
                logger.info(f"Fetching {stat_name} (ID: {stat_id})")

                # Get last 30 days of data
                end_date = today.strftime("%Y%m")
                start_date = (today - timedelta(days=30)).strftime("%Y%m")

                data = self.get_stats_data(
                    stats_id=stat_id,
                    start_date=start_date,
                    end_date=end_date
                )

                if not data:
                    results["errors"].append(f"Failed to fetch {stat_name}")
                    continue

                results["stats_checked"] += 1

                # Search for virus keywords
                keyword_matches = self.search_virus_keywords(data)
                if keyword_matches:
                    results["virus_keywords_found"].extend([m["keyword"] for m in keyword_matches])

                # Save to S3
                s3_key = f"external/estat/{today.year}/{today.month:02d}/{today.day:02d}/{stat_name}_{date_str}.json"
                self._save_to_s3(bucket_name, s3_key, data)
                results["s3_paths"].append(f"s3://{bucket_name}/{s3_key}")

                logger.info(f"Saved {stat_name} to S3: {s3_key}")

            except Exception as e:
                error_msg = f"Error processing {stat_name}: {str(e)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)

        # Save summary to DynamoDB
        self._save_to_dynamodb(results)

        logger.info(f"Daily E-Stat fetch completed: {results['stats_checked']} tables checked, "
                   f"{len(results['virus_keywords_found'])} virus keywords found")

        return results

    def _save_to_s3(self, bucket_name: str, key: str, data: Dict[str, Any]) -> None:
        """
        Save data to S3

        Args:
            bucket_name: S3 bucket name
            key: S3 object key
            data: Data to save
        """
        try:
            self.s3_client.put_object(
                Bucket=bucket_name,
                Key=key,
                Body=json.dumps(data, ensure_ascii=False, indent=2),
                ContentType="application/json"
            )
        except ClientError as e:
            logger.error(f"Failed to save to S3: {e}")
            raise

    def _save_to_dynamodb(self, results: Dict[str, Any]) -> None:
        """
        Save daily fetch results to DynamoDB

        Args:
            results: Daily fetch summary
        """
        try:
            table = self.dynamodb.Table('surveillance-external-updates')

            item = {
                'source#date': f"estat#{results['date']}",
                'update_id': datetime.now().isoformat(),
                'source': results['source'],
                'update_date': results['date'],
                'stats_checked': results['stats_checked'],
                'virus_keywords_found': results.get('virus_keywords_found', []),
                's3_paths': results.get('s3_paths', []),
                'errors': results.get('errors', []),
                'timestamp': int(datetime.now().timestamp())
            }

            table.put_item(Item=item)
            logger.info(f"Saved E-Stat update to DynamoDB: estat#{results['date']}")

        except ClientError as e:
            logger.error(f"Failed to save to DynamoDB: {e}")
            # Don't raise - this is non-critical


if __name__ == "__main__":
    # Example usage
    client = EStatClient()

    # Test: Get stats list
    stats_list = client.get_stats_list(search_keyword="豚")
    print(f"Found {len(stats_list)} statistical tables")

    # Test: Daily fetch (requires S3 bucket)
    # results = client.daily_fetch_livestock_stats(bucket_name="surveillance-data")
    # print(json.dumps(results, indent=2, ensure_ascii=False))
