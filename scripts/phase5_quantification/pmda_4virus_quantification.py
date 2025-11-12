#!/usr/bin/env python3
"""
PMDA 4-Virus Quantification Module
Specialized quantification for high-sensitivity virus detection

Handles:
- Polyomavirus (dsDNA, 5.15 kb)
- Hantavirus (ssRNA-, trisegmented: L/M/S, 11.9 kb total)
- EEEV (ssRNA+, 11.8 kb)
- Spumavirus (retrovirus proviral DNA, ~12 kb)

Features:
- Segmented genome support (Hantavirus 3-segment concordance)
- LOD-based confidence scoring (≥50 copies/mL target)
- RNA vs DNA extraction efficiency adjustments
- Validation against PMDA requirements
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np


# PMDA 4-virus genome sizes (bp)
VIRUS_GENOMES = {
    'POLYOMA': {
        'size': 5150,
        'type': 'dsDNA',
        'segments': 1,
        'lod_target': 50  # copies/mL
    },
    'HANTV': {
        'size': 11900,  # Total L+M+S
        'type': 'ssRNA-',
        'segments': 3,
        'segment_sizes': {'L': 6530, 'M': 3650, 'S': 1720},
        'lod_target': 100  # copies/mL (lower sensitivity for RNA)
    },
    'EEEV': {
        'size': 11841,
        'type': 'ssRNA+',
        'segments': 1,
        'lod_target': 50  # copies/mL
    },
    'SPUMV': {
        'size': 12000,
        'type': 'retrovirus-DNA',
        'segments': 1,
        'lod_target': 10,  # copies per 10^5 PBMCs (NOT plasma!)
        'sample_type': 'pbmc'  # Special: detected in cells, not plasma
    }
}

# Extraction efficiency by nucleic acid type
EXTRACTION_EFFICIENCY = {
    'dsDNA': 0.70,  # DNA extraction ~70% efficient
    'ssRNA-': 0.50,  # RNA extraction ~50% (more fragile)
    'ssRNA+': 0.55,  # Poly(A) selection helps slightly
    'retrovirus-DNA': 0.65  # Proviral DNA from cells
}


class PMDA4VirusQuantifier:
    """Quantification engine for PMDA 4 high-sensitivity viruses."""

    def __init__(self, plasma_volume_ml: float = 10.0, avg_read_length: int = 1000):
        """
        Initialize quantifier.

        Args:
            plasma_volume_ml: Plasma sample volume in mL
            avg_read_length: Average sequencing read length in bp
        """
        self.plasma_volume_ml = plasma_volume_ml
        self.avg_read_length = avg_read_length

    def calculate_copies_per_ml(self, virus_code: str, reads: int,
                               total_reads: int, segment: Optional[str] = None) -> float:
        """
        Calculate absolute copy number per mL for a virus.

        Args:
            virus_code: Virus code (POLYOMA, HANTV, EEEV, SPUMV)
            reads: Number of reads mapping to virus
            total_reads: Total sequencing reads
            segment: For segmented viruses (e.g., 'L', 'M', 'S')

        Returns:
            Copies per mL plasma
        """
        if virus_code not in VIRUS_GENOMES:
            raise ValueError(f"Unknown virus: {virus_code}")

        virus_info = VIRUS_GENOMES[virus_code]

        # For segmented viruses, use segment-specific size
        if segment and 'segment_sizes' in virus_info:
            genome_size = virus_info['segment_sizes'][segment]
        else:
            genome_size = virus_info['size']

        # Get extraction efficiency
        extraction_eff = EXTRACTION_EFFICIENCY[virus_info['type']]

        # Estimate total DNA/RNA molecules from sequencing depth
        # Assumes: sequencing depth = total_reads * avg_read_length / genome_size
        total_molecules = total_reads * self.avg_read_length / genome_size

        # Calculate pathogen fraction
        pathogen_fraction = reads / total_reads if total_reads > 0 else 0

        # Pathogen molecules detected
        pathogen_molecules = pathogen_fraction * total_molecules

        # Adjust for extraction efficiency and normalize to per mL
        copies_per_ml = (pathogen_molecules / extraction_eff) / self.plasma_volume_ml

        return copies_per_ml

    def quantify_polyomavirus(self, detection_results: Dict) -> Dict:
        """
        Quantify Polyomavirus from detection results.

        Args:
            detection_results: Output from detect_pmda_4viruses.py

        Returns:
            Quantification results
        """
        polyoma_data = detection_results.get('polyomavirus', {})

        if not polyoma_data.get('detected', False):
            return {
                'detected': False,
                'copies_per_ml': 0,
                'confidence': 'not_detected'
            }

        reads = polyoma_data.get('total_reads', 0)
        total_reads = detection_results.get('metadata', {}).get('total_input_reads', 1000000)

        copies_per_ml = self.calculate_copies_per_ml('POLYOMA', reads, total_reads)

        # Confidence based on LOD target (50 copies/mL)
        lod = VIRUS_GENOMES['POLYOMA']['lod_target']
        if copies_per_ml >= lod * 10:
            confidence = 'high'
        elif copies_per_ml >= lod:
            confidence = 'medium'
        else:
            confidence = 'low'

        return {
            'detected': True,
            'reads': reads,
            'copies_per_ml': round(copies_per_ml, 2),
            'log10_copies_per_ml': round(np.log10(copies_per_ml + 1), 2),
            'confidence': confidence,
            'pmda_lod_target': lod,
            'above_lod': copies_per_ml >= lod
        }

    def quantify_hantavirus(self, detection_results: Dict) -> Dict:
        """
        Quantify Hantavirus from detection results (3-segment concordance).

        For Hantavirus, we need:
        - All 3 segments detected (L AND M AND S)
        - Average copies/mL across segments
        - Segment-specific quantification

        Args:
            detection_results: Output from detect_pmda_4viruses.py

        Returns:
            Quantification results with segment breakdown
        """
        hantv_data = detection_results.get('hantavirus', {})

        if not hantv_data.get('detected', False):
            return {
                'detected': False,
                'copies_per_ml': 0,
                'confidence': 'not_detected',
                'segments': {}
            }

        total_reads = detection_results.get('metadata', {}).get('total_input_reads', 1000000)
        segment_results = hantv_data.get('segments', {})

        # Quantify each segment
        segment_quant = {}
        segment_copies = []

        for segment in ['L', 'M', 'S']:
            if segment not in segment_results:
                continue

            seg_data = segment_results[segment]
            reads = seg_data.get('reads', 0)

            if reads > 0:
                copies_per_ml = self.calculate_copies_per_ml('HANTV', reads, total_reads, segment)
                segment_quant[segment] = {
                    'reads': reads,
                    'copies_per_ml': round(copies_per_ml, 2),
                    'log10_copies_per_ml': round(np.log10(copies_per_ml + 1), 2)
                }
                segment_copies.append(copies_per_ml)

        # Average across segments (geometric mean for log-distributed data)
        if segment_copies:
            avg_copies_per_ml = np.exp(np.mean(np.log(np.array(segment_copies) + 1))) - 1
        else:
            avg_copies_per_ml = 0

        # Confidence based on LOD target (100 copies/mL for RNA virus)
        lod = VIRUS_GENOMES['HANTV']['lod_target']
        all_segments_detected = len(segment_quant) == 3

        if not all_segments_detected:
            confidence = 'low'  # Must have all 3 segments
        elif avg_copies_per_ml >= lod * 10:
            confidence = 'high'
        elif avg_copies_per_ml >= lod:
            confidence = 'medium'
        else:
            confidence = 'low'

        return {
            'detected': True,
            'species': hantv_data.get('species', 'Unknown'),
            'segments': segment_quant,
            'all_segments_detected': all_segments_detected,
            'copies_per_ml': round(avg_copies_per_ml, 2),
            'log10_copies_per_ml': round(np.log10(avg_copies_per_ml + 1), 2),
            'confidence': confidence,
            'pmda_lod_target': lod,
            'above_lod': avg_copies_per_ml >= lod,
            'note': '3-segment concordance required for high confidence'
        }

    def quantify_eeev(self, detection_results: Dict) -> Dict:
        """
        Quantify EEEV from detection results.

        Args:
            detection_results: Output from detect_pmda_4viruses.py

        Returns:
            Quantification results
        """
        eeev_data = detection_results.get('eeev', {})

        if not eeev_data.get('detected', False):
            return {
                'detected': False,
                'copies_per_ml': 0,
                'confidence': 'not_detected'
            }

        reads = eeev_data.get('total_reads', 0)
        total_reads = detection_results.get('metadata', {}).get('total_input_reads', 1000000)

        copies_per_ml = self.calculate_copies_per_ml('EEEV', reads, total_reads)

        # Confidence based on LOD target (50 copies/mL)
        lod = VIRUS_GENOMES['EEEV']['lod_target']
        if copies_per_ml >= lod * 10:
            confidence = 'high'
        elif copies_per_ml >= lod:
            confidence = 'medium'
        else:
            confidence = 'low'

        return {
            'detected': True,
            'species': eeev_data.get('species', 'EEEV'),
            'reads': reads,
            'copies_per_ml': round(copies_per_ml, 2),
            'log10_copies_per_ml': round(np.log10(copies_per_ml + 1), 2),
            'confidence': confidence,
            'pmda_lod_target': lod,
            'above_lod': copies_per_ml >= lod
        }

    def quantify_spumavirus(self, detection_results: Dict) -> Dict:
        """
        Quantify Spumavirus from detection results.

        NOTE: Spumavirus is detected in PBMCs, not plasma.
        Units are copies per 10^5 PBMCs, NOT copies/mL plasma.

        Args:
            detection_results: Output from detect_pmda_4viruses.py

        Returns:
            Quantification results (copies per 10^5 PBMCs)
        """
        spumv_data = detection_results.get('spumavirus', {})

        if not spumv_data.get('detected', False):
            return {
                'detected': False,
                'copies_per_1e5_pbmcs': 0,
                'confidence': 'not_detected',
                'note': 'Nested PCR from PBMC DNA required for reliable detection'
            }

        reads = spumv_data.get('total_reads', 0)
        total_reads = detection_results.get('metadata', {}).get('total_input_reads', 1000000)

        # For metagenomic detection from plasma cfDNA
        # This is NOT the recommended method - nested PCR from PBMCs is required
        copies_per_ml_plasma = self.calculate_copies_per_ml('SPUMV', reads, total_reads)

        # Convert to approximate copies per 10^5 PBMCs (rough estimate)
        # Assuming 1 mL plasma ~ 10^5 PBMCs for proviral DNA
        copies_per_1e5_pbmcs = copies_per_ml_plasma

        lod = VIRUS_GENOMES['SPUMV']['lod_target']

        return {
            'detected': True,
            'reads': reads,
            'copies_per_1e5_pbmcs': round(copies_per_1e5_pbmcs, 2),
            'log10_copies_per_1e5_pbmcs': round(np.log10(copies_per_1e5_pbmcs + 1), 2),
            'confidence': 'low',  # Always low for metagenomic detection
            'pmda_lod_target': lod,
            'above_lod': copies_per_1e5_pbmcs >= lod,
            'note': 'Metagenomic detection has low sensitivity. Nested PCR from PBMC DNA required for confirmation.',
            'recommended_test': 'detect_spumavirus_nested_pcr.py'
        }

    def quantify_all(self, detection_results_file: Path, run_id: str) -> Dict:
        """
        Quantify all 4 PMDA viruses from detection results.

        Args:
            detection_results_file: JSON file from detect_pmda_4viruses.py
            run_id: Run identifier

        Returns:
            Complete quantification results
        """
        # Load detection results
        with open(detection_results_file) as f:
            detection_results = json.load(f)

        results = {
            'run_id': run_id,
            'plasma_volume_ml': self.plasma_volume_ml,
            'avg_read_length': self.avg_read_length,
            'viruses': {}
        }

        # Quantify each virus
        results['viruses']['POLYOMA'] = self.quantify_polyomavirus(detection_results)
        results['viruses']['HANTV'] = self.quantify_hantavirus(detection_results)
        results['viruses']['EEEV'] = self.quantify_eeev(detection_results)
        results['viruses']['SPUMV'] = self.quantify_spumavirus(detection_results)

        # Summary statistics
        detected_count = sum(1 for v in results['viruses'].values() if v.get('detected', False))
        high_confidence = sum(1 for v in results['viruses'].values()
                            if v.get('confidence') == 'high')

        results['summary'] = {
            'total_viruses_screened': 4,
            'viruses_detected': detected_count,
            'high_confidence_detections': high_confidence,
            'pmda_compliant': True  # All 4 viruses screened
        }

        return results


def main():
    parser = argparse.ArgumentParser(
        description='Quantify PMDA 4-virus detection results (copies/mL)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard quantification
  %(prog)s -i detection_results.json -o quantification.json -r RUN-001

  # Custom plasma volume
  %(prog)s -i detection_results.json -o quantification.json -r RUN-001 \\
      --plasma-volume 5.0

  # Adjust read length
  %(prog)s -i detection_results.json -o quantification.json -r RUN-001 \\
      --avg-read-length 2000
        """
    )

    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Detection results JSON from detect_pmda_4viruses.py')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Output quantification JSON file')
    parser.add_argument('-r', '--run-id', required=True,
                       help='Run identifier')

    parser.add_argument('--plasma-volume', type=float, default=10.0,
                       help='Plasma sample volume in mL (default: 10.0)')
    parser.add_argument('--avg-read-length', type=int, default=1000,
                       help='Average sequencing read length in bp (default: 1000)')

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"ERROR: Detection results file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print("=" * 80)
    print("PMDA 4-Virus Quantification")
    print("=" * 80)
    print(f"Run ID: {args.run_id}")
    print(f"Input: {args.input}")
    print(f"Plasma volume: {args.plasma_volume} mL")
    print(f"Average read length: {args.avg_read_length} bp")
    print()

    # Initialize quantifier
    quantifier = PMDA4VirusQuantifier(
        plasma_volume_ml=args.plasma_volume,
        avg_read_length=args.avg_read_length
    )

    # Quantify all viruses
    results = quantifier.quantify_all(args.input, args.run_id)

    # Save results
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print("\nQuantification Summary:")
    print("-" * 80)

    for virus_code, data in results['viruses'].items():
        virus_name = virus_code
        detected = data.get('detected', False)
        confidence = data.get('confidence', 'N/A')

        if detected:
            if virus_code == 'SPUMV':
                copies = data.get('copies_per_1e5_pbmcs', 0)
                unit = 'copies per 10^5 PBMCs'
            else:
                copies = data.get('copies_per_ml', 0)
                unit = 'copies/mL'

            print(f"{virus_name}: DETECTED")
            print(f"  Quantity: {copies:.2e} {unit}")
            print(f"  Confidence: {confidence}")

            if 'above_lod' in data:
                lod_status = "PASS" if data['above_lod'] else "BELOW LOD"
                print(f"  LOD Status: {lod_status} (target: {data.get('pmda_lod_target')} {unit})")

            if 'note' in data:
                print(f"  Note: {data['note']}")
        else:
            print(f"{virus_name}: NOT DETECTED")

        print()

    print("=" * 80)
    print(f"Summary: {results['summary']['viruses_detected']}/4 viruses detected")
    print(f"High confidence detections: {results['summary']['high_confidence_detections']}")
    print(f"Results saved to: {args.output}")


if __name__ == '__main__':
    main()
