"""
Lambda Function: Daily External Information Collector

Orchestrates daily collection of surveillance data from external sources:
- MAFF (Ministry of Agriculture)
- E-Stat (Government Statistics Portal)
- Academic publications (PubMed/J-STAGE)

Triggered by EventBridge daily schedule (cron: 0 2 * * ? * = 11:00 JST)
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

# Add surveillance modules to path
sys.path.insert(0, '/opt/python')  # Lambda Layer path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from surveillance.external.maff_scraper import MAFFScraper
from surveillance.external.estat_client import EStatClient
from surveillance.external.academic_monitor import AcademicMonitor
from surveillance.alerting.severity_engine import SeverityEngine
from surveillance.alerting.notification_router import NotificationRouter


# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)


# Environment variables
S3_BUCKET = os.environ.get('SURVEILLANCE_BUCKET', 'surveillance-data')
E_STAT_APP_ID = os.environ.get('E_STAT_APP_ID', 'bae1f981a6d093a9676b03c8eea37324b8de421b')
PUBMED_EMAIL = os.environ.get('PUBMED_EMAIL', 'surveillance@example.com')
AWS_REGION = os.environ.get('AWS_REGION', 'ap-northeast-1')


def lambda_handler(event, context):
    """
    Lambda handler for daily external information collection

    Args:
        event: EventBridge schedule event
        context: Lambda context

    Returns:
        Summary of collection results
    """
    logger.info(f"Starting daily external collection at {datetime.now().isoformat()}")
    logger.info(f"Event: {json.dumps(event)}")

    results = {
        'execution_time': datetime.now().isoformat(),
        'bucket': S3_BUCKET,
        'region': AWS_REGION,
        'sources': {},
        'alerts_generated': 0,
        'errors': []
    }

    try:
        # 1. Collect from MAFF
        logger.info("=== Starting MAFF collection ===")
        try:
            maff_scraper = MAFFScraper(region=AWS_REGION)
            maff_results = maff_scraper.daily_check(bucket_name=S3_BUCKET)

            results['sources']['maff'] = {
                'status': 'success',
                'pages_checked': maff_results.get('pages_checked', 0),
                'new_reports': maff_results.get('new_reports', 0),
                'virus_keywords': maff_results.get('virus_keywords_found', []),
                's3_paths': maff_results.get('s3_paths', [])
            }

            # Check for alerts
            if maff_results.get('virus_keywords_found'):
                results['alerts_generated'] += len(maff_results['virus_keywords_found'])
                logger.warning(f"MAFF keywords found: {maff_results['virus_keywords_found']}")

        except Exception as e:
            error_msg = f"MAFF collection failed: {str(e)}"
            logger.error(error_msg)
            results['sources']['maff'] = {'status': 'error', 'error': error_msg}
            results['errors'].append(error_msg)

        # 2. Collect from E-Stat
        logger.info("=== Starting E-Stat collection ===")
        try:
            estat_client = EStatClient(app_id=E_STAT_APP_ID, region=AWS_REGION)
            estat_results = estat_client.daily_fetch_livestock_stats(bucket_name=S3_BUCKET)

            results['sources']['estat'] = {
                'status': 'success',
                'stats_checked': estat_results.get('stats_checked', 0),
                'virus_keywords': estat_results.get('virus_keywords_found', []),
                's3_paths': estat_results.get('s3_paths', [])
            }

            # Check for alerts
            if estat_results.get('virus_keywords_found'):
                results['alerts_generated'] += len(estat_results['virus_keywords_found'])
                logger.warning(f"E-Stat keywords found: {estat_results['virus_keywords_found']}")

        except Exception as e:
            error_msg = f"E-Stat collection failed: {str(e)}"
            logger.error(error_msg)
            results['sources']['estat'] = {'status': 'error', 'error': error_msg}
            results['errors'].append(error_msg)

        # 3. Collect from Academic sources
        logger.info("=== Starting Academic collection ===")
        try:
            academic_monitor = AcademicMonitor(email=PUBMED_EMAIL, region=AWS_REGION)
            academic_results = academic_monitor.daily_monitor(
                bucket_name=S3_BUCKET,
                days_back=1  # Daily search
            )

            results['sources']['academic'] = {
                'status': 'success',
                'viruses_checked': academic_results.get('viruses_checked', 0),
                'total_articles': academic_results.get('total_articles', 0),
                'articles_by_virus': academic_results.get('articles_by_virus', {}),
                's3_paths': academic_results.get('s3_paths', [])
            }

            # New publications are informational, not immediate alerts
            if academic_results.get('total_articles', 0) > 0:
                logger.info(f"Found {academic_results['total_articles']} new academic publications")

        except Exception as e:
            error_msg = f"Academic collection failed: {str(e)}"
            logger.error(error_msg)
            results['sources']['academic'] = {'status': 'error', 'error': error_msg}
            results['errors'].append(error_msg)

        # 4. Process detections through severity engine if keywords found
        if results['alerts_generated'] > 0:
            logger.info(f"Processing {results['alerts_generated']} keyword alerts")

            try:
                severity_engine = SeverityEngine()
                notification_router = NotificationRouter(region=AWS_REGION)

                # Create detections for keyword matches
                detections = []

                # MAFF keywords
                for keyword in results['sources'].get('maff', {}).get('virus_keywords', []):
                    detection = create_detection_from_keyword(keyword, 'maff')
                    classified = severity_engine.classify_detection(detection)
                    detections.append(classified)

                # E-Stat keywords
                for keyword in results['sources'].get('estat', {}).get('virus_keywords', []):
                    detection = create_detection_from_keyword(keyword, 'estat')
                    classified = severity_engine.classify_detection(detection)
                    detections.append(classified)

                # Route alerts
                for detection in detections:
                    if detection['severity'] in ['critical', 'high']:
                        notification_router.route_alert(detection)
                        logger.info(f"Alert routed for {detection['virus_type']} ({detection['severity']})")

                results['detections_processed'] = len(detections)

            except Exception as e:
                error_msg = f"Alert processing failed: {str(e)}"
                logger.error(error_msg)
                results['errors'].append(error_msg)

        # 5. Log summary
        logger.info("=== Collection Summary ===")
        logger.info(f"MAFF: {results['sources'].get('maff', {}).get('status', 'not_run')}")
        logger.info(f"E-Stat: {results['sources'].get('estat', {}).get('status', 'not_run')}")
        logger.info(f"Academic: {results['sources'].get('academic', {}).get('status', 'not_run')}")
        logger.info(f"Alerts generated: {results['alerts_generated']}")
        logger.info(f"Errors: {len(results['errors'])}")

        # Return success
        return {
            'statusCode': 200,
            'body': json.dumps(results, default=str)
        }

    except Exception as e:
        logger.error(f"Fatal error in Lambda execution: {str(e)}", exc_info=True)

        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e),
                'results': results
            }, default=str)
        }


def create_detection_from_keyword(keyword: str, source: str) -> Dict[str, Any]:
    """
    Create detection record from keyword match

    Args:
        keyword: Matched keyword (e.g., "hantavirus:ハンタウイルス")
        source: Source name (maff, estat, academic)

    Returns:
        Detection dictionary
    """
    import uuid

    # Parse virus type from keyword
    virus_type = 'unknown'
    if 'hanta' in keyword.lower() or 'ハンタ' in keyword:
        virus_type = 'hantavirus'
    elif 'polyoma' in keyword.lower() or 'ポリオーマ' in keyword:
        virus_type = 'polyomavirus'
    elif 'spuma' in keyword.lower() or 'スピューマ' in keyword or 'foamy' in keyword.lower():
        virus_type = 'spumavirus'
    elif 'eeev' in keyword.lower() or 'equine' in keyword.lower() or 'ウマ脳炎' in keyword:
        virus_type = 'eeev'

    return {
        'detection_id': str(uuid.uuid4()),
        'timestamp': datetime.now().isoformat(),
        'virus_type': virus_type,
        'source': f'external_{source}',
        'detected': True,
        'keyword_match': keyword,
        'metadata': {
            'keyword_match': True,
            'original_keyword': keyword
        }
    }


if __name__ == "__main__":
    # Local testing
    test_event = {
        'source': 'aws.events',
        'detail-type': 'Scheduled Event',
        'time': datetime.now().isoformat()
    }

    class TestContext:
        def __init__(self):
            self.function_name = 'external_collector'
            self.memory_limit_in_mb = 512
            self.invoked_function_arn = 'arn:aws:lambda:ap-northeast-1:123456789012:function:external_collector'
            self.aws_request_id = 'test-request-id'

    result = lambda_handler(test_event, TestContext())
    print(json.dumps(result, indent=2))
