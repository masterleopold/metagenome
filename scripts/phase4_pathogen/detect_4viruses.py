#!/usr/bin/env python3
"""
4-Virus Detection Module for Surveillance System

Detects Hantavirus, Polyomavirus, Spumavirus, and EEEV from MinION sequencing data.
Integrates with Phase 4 pathogen detection pipeline.

Pattern follows: scripts/phase4_pathogen/perv_typing.py
"""

import argparse
import json
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Any
import pysam
import numpy as np
import boto3
from datetime import datetime


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Target virus markers and detection criteria
VIRUS_MARKERS = {
    'hantavirus': {
        'taxa_names': ['Hantaan', 'Seoul', 'Puumala', 'Hantavirus'],
        'ncbi_taxids': [1980519, 11594, 11596, 1980416],
        'genome_size': 12000,  # ~12kb tri-segmented
        'min_reads': 10,
        'min_coverage': 1.0,  # %
        'specific_motifs': [
            'GCCATGGA',  # L segment conserved region
            'TAGTAGTAGACT'  # Nucleocapsid protein
        ]
    },
    'polyomavirus': {
        'taxa_names': ['Polyomavirus', 'Sus scrofa polyomavirus'],
        'ncbi_taxids': [1891763],  # Sus scrofa polyomavirus 2
        'genome_size': 5000,  # ~5kb circular
        'min_reads': 5,
        'min_coverage': 0.5,
        'specific_motifs': [
            'GAGGCGCCAT',  # VP1 capsid protein
            'CTGGAGGCGG'   # Large T antigen
        ]
    },
    'spumavirus': {
        'taxa_names': ['Spumavirus', 'foamy virus', 'Porcine foamy'],
        'ncbi_taxids': [11791, 35268],  # Spumaretrovirinae
        'genome_size': 12000,  # ~12kb retrovirus
        'min_reads': 10,
        'min_coverage': 1.0,
        'specific_motifs': [
            'TGGAAGACCT',  # pol gene conserved
            'ATGGCAGCCA'   # env region
        ]
    },
    'eeev': {
        'taxa_names': ['Eastern equine encephalitis', 'EEEV', 'Alphavirus'],
        'ncbi_taxids': [11021],
        'genome_size': 11700,  # ~11.7kb
        'min_reads': 10,
        'min_coverage': 1.0,
        'specific_motifs': [
            'ATGGAGGACG',  # nsP1 conserved
            'CACAGACAAG'   # E1 glycoprotein
        ]
    }
}


# Detection thresholds (aligned with Protocol 12 LOD: <100-500 copies/mL)
DETECTION_THRESHOLDS = {
    'critical': 500,   # CRITICAL alert (PERV-level)
    'high': 100,       # HIGH alert
    'medium': 50,      # MEDIUM alert
    'low': 10          # LOW alert (detection only)
}


def detect_4viruses_from_bam(
    bam_file: Path,
    run_id: str,
    sample_id: str,
    output_dir: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Detect target viruses from BAM file

    Args:
        bam_file: Path to BAM file from Phase 4
        run_id: MinION run identifier
        sample_id: Sample identifier
        output_dir: Output directory for results

    Returns:
        Detection results dictionary
    """
    # Validate BAM file
    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if bam_file.stat().st_size == 0:
        raise ValueError(f"Empty BAM file: {bam_file}")

    logger.info(f"Processing BAM file: {bam_file}")

    try:
        bam = pysam.AlignmentFile(str(bam_file), "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file {bam_file}: {e}")

    # Initialize results
    results = {
        'run_id': run_id,
        'sample_id': sample_id,
        'bam_file': str(bam_file),
        'timestamp': datetime.now().isoformat(),
        'viruses': {}
    }

    # Process each target virus
    for virus_name, markers in VIRUS_MARKERS.items():
        logger.info(f"Scanning for {virus_name}...")

        virus_reads = defaultdict(list)
        virus_coverage = defaultdict(int)
        motif_matches = []

        # Scan BAM for virus-aligned reads
        for read in bam.fetch():
            if read.is_unmapped:
                continue

            reference = read.reference_name
            seq = read.query_sequence

            # Check if reference matches target virus
            virus_detected = False
            for taxa_name in markers['taxa_names']:
                if taxa_name.lower() in reference.lower():
                    virus_detected = True
                    break

            if not virus_detected:
                continue

            # Store read information
            virus_reads[read.query_name].append({
                'reference': reference,
                'start': read.reference_start,
                'end': read.reference_end,
                'mapq': read.mapping_quality,
                'length': read.query_length
            })

            # Track coverage
            for pos in range(read.reference_start, read.reference_end):
                virus_coverage[pos] += 1

            # Check for specific motifs
            if seq:
                for motif in markers['specific_motifs']:
                    if motif in seq:
                        motif_matches.append({
                            'motif': motif,
                            'read': read.query_name,
                            'position': seq.find(motif)
                        })

        # Calculate detection statistics
        unique_reads = len(virus_reads)
        detected = unique_reads >= markers['min_reads']

        # Calculate coverage metrics
        coverage_depth = 0
        coverage_breadth = 0
        if virus_coverage:
            coverage_array = list(virus_coverage.values())
            coverage_depth = np.mean(coverage_array)
            coverage_breadth = (len([x for x in coverage_array if x > 0]) /
                               markers['genome_size'] * 100)

        # Estimate copies/mL (simplified - would need spike-in calibration)
        # Rough estimate based on read count and sequencing depth
        estimated_copies_per_ml = unique_reads * 100  # Placeholder formula

        # Determine severity
        severity = determine_severity(estimated_copies_per_ml)

        # Store virus results
        results['viruses'][virus_name] = {
            'detected': detected,
            'unique_reads': unique_reads,
            'total_alignments': sum(len(v) for v in virus_reads.values()),
            'coverage_depth': float(coverage_depth),
            'coverage_breadth': float(coverage_breadth),
            'motif_matches': len(motif_matches),
            'estimated_copies_per_ml': estimated_copies_per_ml,
            'severity': severity,
            'detection_threshold_met': detected,
            'metadata': {
                'genome_size': markers['genome_size'],
                'min_reads_threshold': markers['min_reads'],
                'min_coverage_threshold': markers['min_coverage']
            }
        }

        if detected:
            logger.warning(f"DETECTED: {virus_name} - {unique_reads} reads, "
                         f"{estimated_copies_per_ml} copies/mL (estimated), "
                         f"severity: {severity}")

    bam.close()

    # Determine if any critical detections require immediate alert
    critical_detections = [
        v for v, data in results['viruses'].items()
        if data['detected'] and data['severity'] == 'critical'
    ]

    results['requires_immediate_alert'] = len(critical_detections) > 0
    results['critical_viruses'] = critical_detections

    # Save results
    if output_dir:
        save_results(results, output_dir)

    # Send SNS alert if critical detection
    if results['requires_immediate_alert']:
        send_sns_alert(results)

    return results


def determine_severity(copies_per_ml: float) -> str:
    """
    Determine severity level based on viral load

    Args:
        copies_per_ml: Estimated viral load

    Returns:
        Severity level (critical/high/medium/low)
    """
    if copies_per_ml >= DETECTION_THRESHOLDS['critical']:
        return 'critical'
    elif copies_per_ml >= DETECTION_THRESHOLDS['high']:
        return 'high'
    elif copies_per_ml >= DETECTION_THRESHOLDS['medium']:
        return 'medium'
    else:
        return 'low'


def save_results(results: Dict[str, Any], output_dir: Path) -> None:
    """
    Save detection results to JSON file

    Args:
        results: Detection results
        output_dir: Output directory
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"detect_4viruses_{results['run_id']}.json"

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to: {output_file}")


def send_sns_alert(results: Dict[str, Any]) -> None:
    """
    Send SNS alert for critical detections
    Pattern follows: perv_typing.py:34 SNS alert

    Args:
        results: Detection results with critical viruses
    """
    try:
        sns_client = boto3.client('sns', region_name='ap-northeast-1')

        # SNS Topic ARN (to be configured)
        topic_arn = "arn:aws:sns:ap-northeast-1:ACCOUNT:4virus-critical-alerts"

        message = {
            'alert_type': 'CRITICAL_VIRUS_DETECTION',
            'timestamp': results['timestamp'],
            'run_id': results['run_id'],
            'sample_id': results['sample_id'],
            'critical_viruses': results['critical_viruses'],
            'details': {
                virus: results['viruses'][virus]
                for virus in results['critical_viruses']
            }
        }

        subject = f"CRITICAL: {', '.join(results['critical_viruses'])} detected - {results['run_id']}"

        sns_client.publish(
            TopicArn=topic_arn,
            Subject=subject,
            Message=json.dumps(message, indent=2)
        )

        logger.warning(f"SNS alert sent for critical detection: {results['critical_viruses']}")

    except Exception as e:
        logger.error(f"Failed to send SNS alert: {e}")


def detect_4viruses_from_kraken2(
    kraken_report: Path,
    run_id: str,
    sample_id: str
) -> Dict[str, Any]:
    """
    Detect target viruses from Kraken2 classification report

    Args:
        kraken_report: Path to Kraken2 report
        run_id: Run identifier
        sample_id: Sample identifier

    Returns:
        Detection results
    """
    if not kraken_report.exists():
        raise FileNotFoundError(f"Kraken2 report not found: {kraken_report}")

    logger.info(f"Processing Kraken2 report: {kraken_report}")

    results = {
        'run_id': run_id,
        'sample_id': sample_id,
        'kraken_report': str(kraken_report),
        'timestamp': datetime.now().isoformat(),
        'viruses': {}
    }

    with open(kraken_report, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            percentage = float(parts[0].strip())
            reads_clade = int(parts[1].strip())
            reads_taxon = int(parts[2].strip())
            rank = parts[3].strip()
            taxid = int(parts[4].strip())
            name = parts[5].strip()

            # Check against target viruses
            for virus_name, markers in VIRUS_MARKERS.items():
                # Check by taxonomy ID
                if taxid in markers['ncbi_taxids']:
                    results['viruses'][virus_name] = {
                        'detected': True,
                        'source': 'kraken2',
                        'taxid': taxid,
                        'taxa_name': name,
                        'reads': reads_taxon,
                        'percentage': percentage,
                        'rank': rank
                    }
                    logger.warning(f"Detected {virus_name}: {name} ({reads_taxon} reads)")

    return results


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description="Detect 4 surveillance viruses from MinION Phase 4 data"
    )
    parser.add_argument(
        '--bam', type=Path, required=True,
        help="BAM file from Phase 4 pathogen detection"
    )
    parser.add_argument(
        '--run-id', type=str, required=True,
        help="MinION run identifier"
    )
    parser.add_argument(
        '--sample-id', type=str, required=True,
        help="Sample identifier"
    )
    parser.add_argument(
        '--output-dir', type=Path, default=Path('.'),
        help="Output directory for results"
    )
    parser.add_argument(
        '--kraken-report', type=Path,
        help="Optional: Kraken2 report for additional validation"
    )

    args = parser.parse_args()

    # Run BAM-based detection
    results = detect_4viruses_from_bam(
        bam_file=args.bam,
        run_id=args.run_id,
        sample_id=args.sample_id,
        output_dir=args.output_dir
    )

    # Optional: Cross-check with Kraken2
    if args.kraken_report:
        kraken_results = detect_4viruses_from_kraken2(
            kraken_report=args.kraken_report,
            run_id=args.run_id,
            sample_id=args.sample_id
        )
        results['kraken2_validation'] = kraken_results

    # Print summary
    print(json.dumps(results, indent=2))

    # Exit with error code if critical detection
    if results.get('requires_immediate_alert'):
        logger.error("CRITICAL DETECTION - immediate action required")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
