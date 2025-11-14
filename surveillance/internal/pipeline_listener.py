"""
Internal Pipeline Listener for Real-time Pathogen Detection

Monitors MinION pipeline Phase 4 results for target virus detections.
Triggered by S3 events and DynamoDB Streams.

Integration Points:
- S3: s3://minion-data/results/phase4/*/
- DynamoDB: pipeline-metadata table
- Existing: scripts/phase4_pathogen/perv_typing.py pattern
"""

import os
import json
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any
from pathlib import Path
import boto3
from botocore.exceptions import ClientError


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PipelineListener:
    """
    Listener for MinION pipeline Phase 4 pathogen detection results

    Monitors S3 for new Phase 4 result files and processes them
    to detect target viruses.
    """

    # Target viruses configuration
    TARGET_VIRUSES = {
        "hantavirus": {
            "taxa_names": [
                "Hantaan virus", "Hantavirus", "Seoul virus", "Puumala virus",
                "Dobrava-Belgrade virus", "Sin Nombre virus"
            ],
            "ncbi_taxids": [1980519, 11594, 11596, 1980416],  # Hantavirus family
            "keywords": ["hanta", "hantaan", "seoul virus"]
        },
        "polyomavirus": {
            "taxa_names": [
                "Sus scrofa polyomavirus", "Polyomavirus", "Porcine polyomavirus"
            ],
            "ncbi_taxids": [1891763],  # Sus scrofa polyomavirus 2
            "keywords": ["polyoma", "polyomavirus"]
        },
        "spumavirus": {
            "taxa_names": [
                "Porcine foamy virus", "Spumavirus", "Simian foamy virus",
                "Bovine foamy virus", "Feline foamy virus"
            ],
            "ncbi_taxids": [11791, 35268],  # Spumaretrovirinae
            "keywords": ["spuma", "foamy virus", "retrovirus"]
        },
        "eeev": {
            "taxa_names": [
                "Eastern equine encephalitis virus", "EEEV",
                "Eastern equine encephalomyelitis virus"
            ],
            "ncbi_taxids": [11021],  # EEEV
            "keywords": ["eastern equine", "EEEV", "alphavirus"]
        }
    }

    # Phase 4 result file patterns
    RESULT_FILES = {
        "kraken2": "kraken2_report.txt",
        "blast": "blast_results.xml",
        "perv_typing": "perv_typing_results.json",
        "pmda_4virus": "detect_pmda_4viruses_results.json"  # If exists
    }

    def __init__(self, region: str = "ap-northeast-1"):
        """
        Initialize pipeline listener

        Args:
            region: AWS region
        """
        self.region = region

        # Initialize AWS clients
        self.s3_client = boto3.client('s3', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)

        logger.info("Initialized Pipeline Listener")

    def process_s3_event(self, event: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process S3 event notification (Lambda trigger)

        Args:
            event: S3 event from Lambda

        Returns:
            Processing summary
        """
        logger.info(f"Processing S3 event: {json.dumps(event)}")

        results = {
            'timestamp': datetime.now().isoformat(),
            'files_processed': 0,
            'detections': [],
            'errors': []
        }

        # Parse S3 event records
        for record in event.get('Records', []):
            try:
                bucket = record['s3']['bucket']['name']
                key = record['s3']['object']['key']

                logger.info(f"Processing S3 object: s3://{bucket}/{key}")

                # Check if this is a Phase 4 result file
                if 'phase4' not in key:
                    logger.debug(f"Skipping non-Phase4 file: {key}")
                    continue

                # Determine file type and process accordingly
                if 'kraken2_report' in key:
                    detections = self.process_kraken2_report(bucket, key)
                elif 'blast_results' in key:
                    detections = self.process_blast_results(bucket, key)
                elif 'pmda_4viruses' in key:
                    detections = self.process_pmda_4virus_results(bucket, key)
                else:
                    logger.debug(f"Unknown file type: {key}")
                    continue

                results['files_processed'] += 1
                if detections:
                    results['detections'].extend(detections)

            except Exception as e:
                error_msg = f"Error processing record: {str(e)}"
                logger.error(error_msg)
                results['errors'].append(error_msg)

        # Save detections to database
        if results['detections']:
            self._save_detections(results['detections'])

        return results

    def process_kraken2_report(self, bucket: str, key: str) -> List[Dict[str, Any]]:
        """
        Process Kraken2 classification report

        Args:
            bucket: S3 bucket name
            key: S3 object key

        Returns:
            List of virus detections
        """
        try:
            # Download Kraken2 report
            response = self.s3_client.get_object(Bucket=bucket, Key=key)
            report_content = response['Body'].read().decode('utf-8')

            detections = []

            # Parse Kraken2 report format:
            # percentage  reads_clade  reads_taxon  rank  taxid  name
            for line in report_content.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) < 6:
                    continue

                percentage = float(parts[0].strip())
                reads_clade = int(parts[1].strip())
                reads_taxon = int(parts[2].strip())
                rank = parts[3].strip()
                taxid = int(parts[4].strip())
                name = parts[5].strip()

                # Check against target viruses
                for virus_type, config in self.TARGET_VIRUSES.items():
                    # Check by taxid
                    if taxid in config['ncbi_taxids']:
                        detection = self._create_detection(
                            virus_type=virus_type,
                            source='kraken2',
                            taxid=taxid,
                            taxa_name=name,
                            reads=reads_taxon,
                            percentage=percentage,
                            s3_path=f"s3://{bucket}/{key}",
                            metadata={'rank': rank, 'reads_clade': reads_clade}
                        )
                        detections.append(detection)
                        logger.warning(f"Detected {virus_type} in Kraken2 report: {name} ({reads_taxon} reads)")

                    # Check by name keywords
                    elif any(keyword.lower() in name.lower() for keyword in config['keywords']):
                        detection = self._create_detection(
                            virus_type=virus_type,
                            source='kraken2',
                            taxid=taxid,
                            taxa_name=name,
                            reads=reads_taxon,
                            percentage=percentage,
                            s3_path=f"s3://{bucket}/{key}",
                            metadata={'rank': rank, 'reads_clade': reads_clade, 'keyword_match': True}
                        )
                        detections.append(detection)
                        logger.warning(f"Detected {virus_type} (keyword match) in Kraken2: {name}")

            return detections

        except Exception as e:
            logger.error(f"Failed to process Kraken2 report {key}: {e}")
            return []

    def process_blast_results(self, bucket: str, key: str) -> List[Dict[str, Any]]:
        """
        Process BLAST alignment results (XML format)

        Args:
            bucket: S3 bucket name
            key: S3 object key

        Returns:
            List of virus detections
        """
        try:
            from xml.etree import ElementTree as ET

            # Download BLAST results
            response = self.s3_client.get_object(Bucket=bucket, Key=key)
            xml_content = response['Body'].read()

            root = ET.fromstring(xml_content)
            detections = []

            # Parse BLAST XML
            for iteration in root.findall('.//Iteration'):
                for hit in iteration.findall('.//Hit'):
                    hit_def = hit.findtext('Hit_def', '')
                    hit_accession = hit.findtext('Hit_accession', '')
                    hit_len = int(hit.findtext('Hit_len', 0))

                    # Get best HSP (High-scoring Segment Pair)
                    hsp = hit.find('.//Hsp')
                    if hsp is None:
                        continue

                    identity = int(hsp.findtext('Hsp_identity', 0))
                    align_len = int(hsp.findtext('Hsp_align-len', 0))
                    evalue = float(hsp.findtext('Hsp_evalue', 1.0))

                    identity_pct = (identity / align_len * 100) if align_len > 0 else 0

                    # Check against target viruses
                    for virus_type, config in self.TARGET_VIRUSES.items():
                        if any(name.lower() in hit_def.lower() for name in config['taxa_names']):
                            detection = self._create_detection(
                                virus_type=virus_type,
                                source='blast',
                                taxid=None,
                                taxa_name=hit_def,
                                reads=None,
                                percentage=identity_pct,
                                s3_path=f"s3://{bucket}/{key}",
                                metadata={
                                    'accession': hit_accession,
                                    'evalue': evalue,
                                    'identity': identity,
                                    'align_length': align_len,
                                    'hit_length': hit_len
                                }
                            )
                            detections.append(detection)
                            logger.warning(f"Detected {virus_type} in BLAST: {hit_def} ({identity_pct:.2f}% identity)")

            return detections

        except Exception as e:
            logger.error(f"Failed to process BLAST results {key}: {e}")
            return []

    def process_pmda_4virus_results(self, bucket: str, key: str) -> List[Dict[str, Any]]:
        """
        Process PMDA 4-virus detection results (if exists from Protocol 11)

        Args:
            bucket: S3 bucket name
            key: S3 object key

        Returns:
            List of virus detections
        """
        try:
            # Download results JSON
            response = self.s3_client.get_object(Bucket=bucket, Key=key)
            results = json.loads(response['Body'].read().decode('utf-8'))

            detections = []

            for virus_type in ['polyomavirus', 'hantavirus', 'eeev', 'spumavirus']:
                if virus_type in results and results[virus_type].get('detected'):
                    detection = self._create_detection(
                        virus_type=virus_type,
                        source='pmda_4virus_protocol',
                        taxid=None,
                        taxa_name=results[virus_type].get('name', virus_type),
                        reads=results[virus_type].get('reads'),
                        percentage=None,
                        s3_path=f"s3://{bucket}/{key}",
                        metadata={
                            'copies_per_ml': results[virus_type].get('copies_per_ml'),
                            'lod_threshold': results[virus_type].get('lod_threshold'),
                            'detection_method': results[virus_type].get('method')
                        }
                    )
                    detections.append(detection)
                    logger.warning(f"Detected {virus_type} in PMDA 4-virus protocol: {results[virus_type].get('copies_per_ml')} copies/mL")

            return detections

        except Exception as e:
            logger.error(f"Failed to process PMDA 4-virus results {key}: {e}")
            return []

    def _create_detection(
        self,
        virus_type: str,
        source: str,
        taxid: Optional[int],
        taxa_name: str,
        reads: Optional[int],
        percentage: Optional[float],
        s3_path: str,
        metadata: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Create standardized detection record

        Args:
            virus_type: Type of virus detected
            source: Detection source (kraken2, blast, etc.)
            taxid: NCBI taxonomy ID
            taxa_name: Taxonomic name
            reads: Number of reads
            percentage: Percentage/identity
            s3_path: S3 path to source file
            metadata: Additional metadata

        Returns:
            Detection record
        """
        import uuid

        detection = {
            'detection_id': str(uuid.uuid4()),
            'timestamp': datetime.now().isoformat(),
            'virus_type': virus_type,
            'source': f"internal_pipeline_{source}",
            'taxid': taxid,
            'taxa_name': taxa_name,
            'reads': reads,
            'percentage': percentage,
            's3_path': s3_path,
            'metadata': metadata,
            'severity': None,  # Will be calculated by severity engine
            'alert_sent': False
        }

        return detection

    def _save_detections(self, detections: List[Dict[str, Any]]) -> None:
        """
        Save detections to DynamoDB

        Args:
            detections: List of detection records
        """
        try:
            table = self.dynamodb.Table('surveillance-detections')

            for detection in detections:
                table.put_item(Item={
                    **detection,
                    'timestamp_sort': int(datetime.fromisoformat(detection['timestamp']).timestamp())
                })

                logger.info(f"Saved detection to DynamoDB: {detection['virus_type']} ({detection['detection_id']})")

        except ClientError as e:
            logger.error(f"Failed to save detections to DynamoDB: {e}")


def lambda_handler(event, context):
    """
    AWS Lambda handler for S3 event processing

    Args:
        event: S3 event notification
        context: Lambda context

    Returns:
        Processing results
    """
    listener = PipelineListener()
    results = listener.process_s3_event(event)

    logger.info(f"Lambda processing complete: {results['files_processed']} files, "
                f"{len(results['detections'])} detections")

    return {
        'statusCode': 200,
        'body': json.dumps(results)
    }


if __name__ == "__main__":
    # Example: Test S3 event processing
    test_event = {
        'Records': [{
            's3': {
                'bucket': {'name': 'minion-data'},
                'object': {'key': 'results/phase4/RUN-001/kraken2_report.txt'}
            }
        }]
    }

    listener = PipelineListener()
    results = listener.process_s3_event(test_event)
    print(json.dumps(results, indent=2))
