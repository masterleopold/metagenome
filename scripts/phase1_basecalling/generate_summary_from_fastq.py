#!/usr/bin/env python3
"""
MinION Metagenomics Pipeline - Generate Summary from FASTQ
Extracts quality metrics and statistics from basecalled FASTQ files
"""

import argparse
import gzip
import json
import logging
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio import SeqIO
import boto3
import psycopg2
from psycopg2.extras import RealDictCursor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class FASTQAnalyzer:
    """Analyzes FASTQ files to generate sequencing summaries."""

    def __init__(self, run_id: str, db_config: Dict = None):
        self.run_id = run_id
        self.db_config = db_config
        self.metrics = defaultdict(list)
        self.summary = {}

    def process_fastq(self, fastq_path: Path) -> Dict:
        """
        Process FASTQ file and extract metrics.

        Args:
            fastq_path: Path to FASTQ file (can be gzipped)

        Returns:
            Dictionary containing extracted metrics
        """
        logger.info(f"Processing {fastq_path}")

        # Determine if file is gzipped
        if str(fastq_path).endswith('.gz'):
            handle = gzip.open(fastq_path, 'rt')
        else:
            handle = open(fastq_path, 'r')

        read_count = 0
        total_bases = 0
        read_lengths = []
        quality_scores = []

        try:
            for record in SeqIO.parse(handle, "fastq"):
                read_count += 1
                seq_len = len(record.seq)
                total_bases += seq_len
                read_lengths.append(seq_len)

                # Extract quality scores
                quals = record.letter_annotations["phred_quality"]
                quality_scores.extend(quals)

                # Extract additional metadata from header if available
                if '|' in record.description:
                    self._parse_header_metadata(record.description)

                # Log progress
                if read_count % 10000 == 0:
                    logger.info(f"Processed {read_count:,} reads...")

        finally:
            handle.close()

        logger.info(f"Total reads processed: {read_count:,}")

        # Calculate statistics
        self.summary = self._calculate_statistics(
            read_count, total_bases, read_lengths, quality_scores
        )

        return self.summary

    def _parse_header_metadata(self, header: str):
        """Parse ONT-specific metadata from FASTQ header."""
        # Example: @read_id runid=abc123 sampleid=sample1 ...
        parts = header.split()
        for part in parts[1:]:
            if '=' in part:
                key, value = part.split('=', 1)
                self.metrics[key].append(value)

    def _calculate_statistics(
        self,
        read_count: int,
        total_bases: int,
        read_lengths: List[int],
        quality_scores: List[int]
    ) -> Dict:
        """Calculate comprehensive statistics from collected metrics."""

        if read_count == 0:
            return {"error": "No reads found in file"}

        read_lengths = np.array(read_lengths)
        quality_scores = np.array(quality_scores)

        # Calculate N50
        sorted_lengths = np.sort(read_lengths)[::-1]
        cumsum = np.cumsum(sorted_lengths)
        n50_idx = np.where(cumsum >= total_bases / 2)[0][0]
        n50 = sorted_lengths[n50_idx]

        # Calculate quality score distributions
        q_scores = {
            f"q{threshold}": np.sum(quality_scores >= threshold) / len(quality_scores) * 100
            for threshold in [5, 7, 10, 15, 20, 30]
        }

        summary = {
            "run_id": self.run_id,
            "analysis_timestamp": datetime.now().isoformat(),
            "read_metrics": {
                "total_reads": int(read_count),
                "total_bases": int(total_bases),
                "total_gigabases": round(total_bases / 1e9, 2),
            },
            "length_metrics": {
                "mean_length": round(float(np.mean(read_lengths)), 1),
                "median_length": round(float(np.median(read_lengths)), 1),
                "min_length": int(np.min(read_lengths)),
                "max_length": int(np.max(read_lengths)),
                "n50": int(n50),
                "length_percentiles": {
                    "p10": int(np.percentile(read_lengths, 10)),
                    "p25": int(np.percentile(read_lengths, 25)),
                    "p50": int(np.percentile(read_lengths, 50)),
                    "p75": int(np.percentile(read_lengths, 75)),
                    "p90": int(np.percentile(read_lengths, 90)),
                }
            },
            "quality_metrics": {
                "mean_qscore": round(float(np.mean(quality_scores)), 2),
                "median_qscore": round(float(np.median(quality_scores)), 2),
                "min_qscore": round(float(np.min(quality_scores)), 2),
                "max_qscore": round(float(np.max(quality_scores)), 2),
                "qscore_percentiles": {
                    "p10": round(float(np.percentile(quality_scores, 10)), 2),
                    "p25": round(float(np.percentile(quality_scores, 25)), 2),
                    "p50": round(float(np.percentile(quality_scores, 50)), 2),
                    "p75": round(float(np.percentile(quality_scores, 75)), 2),
                    "p90": round(float(np.percentile(quality_scores, 90)), 2),
                },
                "qscore_thresholds": q_scores
            }
        }

        # Add duplex metrics if available
        if 'duplex' in self.metrics:
            duplex_count = sum(1 for x in self.metrics['duplex'] if x == 'true')
            summary['duplex_metrics'] = {
                "duplex_reads": duplex_count,
                "duplex_rate": round(duplex_count / read_count * 100, 2)
            }

        return summary

    def update_database(self):
        """Update database with calculated metrics."""
        if not self.db_config or not self.summary:
            return

        try:
            conn = psycopg2.connect(**self.db_config)
            cur = conn.cursor()

            # Update workflow execution
            cur.execute("""
                UPDATE workflow_executions
                SET total_reads = %(total_reads)s,
                    total_bases_gb = %(total_gb)s,
                    mean_qscore = %(mean_q)s,
                    median_qscore = %(median_q)s,
                    n50 = %(n50)s,
                    median_read_length = %(median_len)s,
                    updated_at = CURRENT_TIMESTAMP
                WHERE run_id = %(run_id)s
            """, {
                'run_id': self.run_id,
                'total_reads': self.summary['read_metrics']['total_reads'],
                'total_gb': self.summary['read_metrics']['total_gigabases'],
                'mean_q': self.summary['quality_metrics']['mean_qscore'],
                'median_q': self.summary['quality_metrics']['median_qscore'],
                'n50': self.summary['length_metrics']['n50'],
                'median_len': self.summary['length_metrics']['median_length']
            })

            # Insert QC metrics
            for metric_name, metric_value in [
                ('total_reads', self.summary['read_metrics']['total_reads']),
                ('mean_qscore', self.summary['quality_metrics']['mean_qscore']),
                ('n50', self.summary['length_metrics']['n50']),
                ('q30_percent', self.summary['quality_metrics']['qscore_thresholds'].get('q30', 0))
            ]:
                cur.execute("""
                    INSERT INTO qc_metrics (run_id, phase, metric_name, metric_value, passed)
                    VALUES (%s, %s, %s, %s, %s)
                    ON CONFLICT (run_id, phase, metric_name) DO UPDATE
                    SET metric_value = EXCLUDED.metric_value,
                        passed = EXCLUDED.passed
                """, (
                    self.run_id,
                    'basecalling',
                    metric_name,
                    metric_value,
                    self._check_qc_pass(metric_name, metric_value)
                ))

            conn.commit()
            logger.info("Database updated successfully")

        except Exception as e:
            logger.error(f"Database update failed: {e}")
        finally:
            if 'conn' in locals():
                conn.close()

    def _check_qc_pass(self, metric_name: str, value: float) -> bool:
        """Check if QC metric passes threshold."""
        thresholds = {
            'total_reads': 100000,  # Minimum 100k reads
            'mean_qscore': 9,       # Minimum Q9
            'n50': 200,             # Minimum N50 of 200bp
            'q30_percent': 5        # At least 5% Q30 reads
        }
        return value >= thresholds.get(metric_name, 0)

    def save_summary(self, output_path: Path):
        """Save summary to JSON file."""
        with open(output_path, 'w') as f:
            json.dump(self.summary, f, indent=2)
        logger.info(f"Summary saved to {output_path}")

    def upload_to_s3(self, output_path: Path, s3_bucket: str):
        """Upload summary to S3."""
        if not s3_bucket:
            return

        try:
            s3 = boto3.client('s3')
            s3_key = f"runs/{self.run_id}/basecalling/sequencing_summary.json"

            s3.upload_file(
                str(output_path),
                s3_bucket,
                s3_key,
                ExtraArgs={
                    'ContentType': 'application/json',
                    'Metadata': {
                        'run_id': self.run_id,
                        'phase': 'basecalling'
                    }
                }
            )
            logger.info(f"Summary uploaded to s3://{s3_bucket}/{s3_key}")

        except Exception as e:
            logger.error(f"S3 upload failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate summary statistics from basecalled FASTQ files"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=Path,
        help='Input FASTQ file (can be gzipped)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        type=Path,
        help='Output JSON file for summary'
    )
    parser.add_argument(
        '-r', '--run-id',
        required=True,
        help='Unique run identifier'
    )
    parser.add_argument(
        '--db-host',
        help='Database host'
    )
    parser.add_argument(
        '--db-name',
        default='minion_metadata',
        help='Database name'
    )
    parser.add_argument(
        '--db-user',
        default='minion_user',
        help='Database user'
    )
    parser.add_argument(
        '--db-password',
        help='Database password'
    )
    parser.add_argument(
        '--s3-bucket',
        help='S3 bucket for uploading summary'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Prepare database config
    db_config = None
    if args.db_host:
        db_config = {
            'host': args.db_host,
            'database': args.db_name,
            'user': args.db_user,
            'password': args.db_password
        }

    # Process FASTQ
    analyzer = FASTQAnalyzer(args.run_id, db_config)
    summary = analyzer.process_fastq(args.input)

    # Save and upload results
    analyzer.save_summary(args.output)
    analyzer.update_database()
    analyzer.upload_to_s3(args.output, args.s3_bucket)

    # Print summary
    print("\n" + "="*60)
    print("Sequencing Summary")
    print("="*60)
    print(f"Total reads: {summary['read_metrics']['total_reads']:,}")
    print(f"Total bases: {summary['read_metrics']['total_gigabases']:.2f} Gb")
    print(f"Mean length: {summary['length_metrics']['mean_length']:.1f} bp")
    print(f"N50: {summary['length_metrics']['n50']:,} bp")
    print(f"Mean Q-score: {summary['quality_metrics']['mean_qscore']:.2f}")
    if 'duplex_metrics' in summary:
        print(f"Duplex rate: {summary['duplex_metrics']['duplex_rate']:.2f}%")
    print("="*60)


if __name__ == '__main__':
    main()