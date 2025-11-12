#!/usr/bin/env python3
"""
Phase 0: Sample Preparation Router
Routes samples to appropriate nucleic acid extraction and prep workflows based on target pathogens.

Handles:
- DNA viruses (Polyomavirus, etc.)
- RNA viruses with poly(A) (EEEV, etc.)
- RNA viruses without poly(A) (Hantavirus, etc.)
- Retrovirus proviral DNA (Spumavirus, PERV, etc.)
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List
from enum import Enum

class VirusType(Enum):
    """Virus genome types requiring different sample prep workflows."""
    DNA = "dsDNA"
    RNA_POLYA_PLUS = "ssRNA_polyA+"
    RNA_POLYA_MINUS = "ssRNA_polyA-"
    PROVIRUS = "proviral_DNA"

class SampleType(Enum):
    """Sample types for different virus detection."""
    PLASMA_CFDNA = "plasma_cfDNA"
    PLASMA_CFRNA = "plasma_cfRNA"
    PBMC_GENOMIC_DNA = "PBMC_genomic_DNA"

# PMDA 4-virus classification
PMDA_4VIRUS_CONFIG = {
    "polyomavirus": {
        "genome_type": VirusType.DNA,
        "sample_type": SampleType.PLASMA_CFDNA,
        "extraction_method": "dual_dna_rna",
        "host_depletion": "cpg_methylation",
        "library_prep": "dna_ligation",
        "pmda_line": "35/58"
    },
    "hantavirus": {
        "genome_type": VirusType.RNA_POLYA_MINUS,
        "sample_type": SampleType.PLASMA_CFRNA,
        "extraction_method": "dual_dna_rna",
        "host_depletion": "rrna_depletion",
        "library_prep": "amplicon_rt_pcr",  # or smart9n_cdna
        "pmda_line": "31/54"
    },
    "eeev": {
        "genome_type": VirusType.RNA_POLYA_PLUS,
        "sample_type": SampleType.PLASMA_CFRNA,
        "extraction_method": "dual_dna_rna",
        "host_depletion": "polya_selection",
        "library_prep": "direct_rna",  # or direct_cdna
        "pmda_line": "32/55"
    },
    "spumavirus": {
        "genome_type": VirusType.PROVIRUS,
        "sample_type": SampleType.PBMC_GENOMIC_DNA,
        "extraction_method": "pbmc_genomic_dna",
        "host_depletion": "none",  # Cannot deplete host for integrated provirus
        "library_prep": "nested_pcr",
        "pmda_line": "154"
    }
}

def determine_extraction_workflow(target_viruses: List[str]) -> Dict[str, any]:
    """
    Determine the optimal extraction and prep workflow based on target viruses.

    Args:
        target_viruses: List of target virus names (e.g., ['polyomavirus', 'hantavirus'])

    Returns:
        Dict containing workflow recommendations
    """
    virus_configs = [PMDA_4VIRUS_CONFIG[v] for v in target_viruses if v in PMDA_4VIRUS_CONFIG]

    if not virus_configs:
        return {
            "error": "No valid target viruses specified",
            "valid_targets": list(PMDA_4VIRUS_CONFIG.keys())
        }

    # Analyze required sample types
    sample_types_needed = set(config["sample_type"] for config in virus_configs)
    genome_types = set(config["genome_type"] for config in virus_configs)

    # Determine extraction strategy
    needs_plasma = any(st in [SampleType.PLASMA_CFDNA, SampleType.PLASMA_CFRNA]
                      for st in sample_types_needed)
    needs_pbmc = SampleType.PBMC_GENOMIC_DNA in sample_types_needed
    needs_dna = VirusType.DNA in genome_types or VirusType.PROVIRUS in genome_types
    needs_rna = any(vt in [VirusType.RNA_POLYA_PLUS, VirusType.RNA_POLYA_MINUS]
                   for vt in genome_types)

    workflow = {
        "target_viruses": target_viruses,
        "blood_collection": {
            "volume_ml": 10 if needs_rna else 5,  # RNA viruses need more volume
            "tube_type": "EDTA",
            "rnase_inhibitor": needs_rna,
            "rnase_inhibitor_concentration": "20 U/mL" if needs_rna else None
        },
        "sample_separation": {
            "plasma_required": needs_plasma,
            "pbmc_required": needs_pbmc,
            "plasma_volume_ml": 5 if needs_plasma else 0,
            "pbmc_cell_count": "1-5×10⁶" if needs_pbmc else "0"
        },
        "extraction": {
            "plasma_extraction": "dual_dna_rna" if (needs_dna and needs_rna) else (
                "cfDNA_only" if needs_dna else "cfRNA_only"
            ),
            "pbmc_extraction": "genomic_dna" if needs_pbmc else None
        },
        "processing_paths": []
    }

    # Define processing paths for each virus
    for virus in target_viruses:
        if virus in PMDA_4VIRUS_CONFIG:
            config = PMDA_4VIRUS_CONFIG[virus]
            path = {
                "virus": virus,
                "genome_type": config["genome_type"].value,
                "sample_type": config["sample_type"].value,
                "host_depletion_method": config["host_depletion"],
                "library_prep_method": config["library_prep"],
                "estimated_lod_copies_ml": get_lod_estimate(virus)
            }
            workflow["processing_paths"].append(path)

    return workflow

def get_lod_estimate(virus: str) -> str:
    """Get estimated LOD for each virus with current protocol."""
    lod_map = {
        "polyomavirus": "50-500 copies/mL (with probe capture)",
        "hantavirus": "100-500 copies/mL (amplicon approach)",
        "eeev": "50-500 copies/mL (poly(A) selection + capture)",
        "spumavirus": "1-10 copies/10⁵ PBMCs (nested PCR)"
    }
    return lod_map.get(virus, "Unknown")

def validate_input_requirements(workflow: Dict, sample_metadata: Dict) -> List[str]:
    """
    Validate that sample metadata meets workflow requirements.

    Returns:
        List of validation errors (empty if valid)
    """
    errors = []

    # Check blood volume
    required_volume = workflow["blood_collection"]["volume_ml"]
    actual_volume = sample_metadata.get("blood_volume_ml", 0)
    if actual_volume < required_volume:
        errors.append(f"Insufficient blood volume: {actual_volume} mL < {required_volume} mL required")

    # Check RNase inhibitor for RNA viruses
    if workflow["blood_collection"]["rnase_inhibitor"]:
        if not sample_metadata.get("rnase_inhibitor_added", False):
            errors.append("RNase inhibitor required but not added to blood collection")

    # Check plasma volume if needed
    if workflow["sample_separation"]["plasma_required"]:
        plasma_vol = sample_metadata.get("plasma_volume_ml", 0)
        required_plasma = workflow["sample_separation"]["plasma_volume_ml"]
        if plasma_vol < required_plasma:
            errors.append(f"Insufficient plasma: {plasma_vol} mL < {required_plasma} mL required")

    # Check PBMC count if needed
    if workflow["sample_separation"]["pbmc_required"]:
        pbmc_count = sample_metadata.get("pbmc_count", 0)
        if pbmc_count < 1e6:
            errors.append(f"Insufficient PBMCs: {pbmc_count} < 1×10⁶ required")

    return errors

def generate_extraction_protocol(workflow: Dict, output_dir: Path) -> Dict[str, Path]:
    """
    Generate detailed extraction protocol files for each processing path.

    Returns:
        Dict mapping virus names to protocol file paths
    """
    protocol_files = {}

    for path in workflow["processing_paths"]:
        virus = path["virus"]
        protocol_file = output_dir / f"{virus}_extraction_protocol.json"

        protocol = {
            "virus": virus,
            "genome_type": path["genome_type"],
            "sample_type": path["sample_type"],
            "steps": generate_protocol_steps(path)
        }

        with open(protocol_file, 'w') as f:
            json.dump(protocol, f, indent=2)

        protocol_files[virus] = protocol_file

    return protocol_files

def generate_protocol_steps(path: Dict) -> List[Dict]:
    """Generate detailed protocol steps for a processing path."""
    steps = []

    # Step 1: Extraction
    if path["sample_type"] == "plasma_cfDNA":
        steps.append({
            "step": 1,
            "name": "cfDNA Extraction",
            "kit": "Zymo Quick-cfDNA/RNA Serum & Plasma Kit",
            "input": "500 μL plasma",
            "output": "50 μL cfDNA elute",
            "expected_yield": "5-50 ng"
        })
    elif path["sample_type"] == "plasma_cfRNA":
        steps.append({
            "step": 1,
            "name": "cfRNA Extraction",
            "kit": "Zymo Quick-cfDNA/RNA Serum & Plasma Kit",
            "input": "500 μL plasma",
            "output": "50 μL cfRNA elute",
            "expected_yield": "1-10 ng",
            "critical": "RNase-free technique required"
        })
    elif path["sample_type"] == "PBMC_genomic_DNA":
        steps.append({
            "step": 1,
            "name": "PBMC Genomic DNA Extraction",
            "kit": "QIAGEN QIAamp DNA Blood Mini Kit",
            "input": "1-5×10⁶ PBMCs",
            "output": "50 μL genomic DNA",
            "expected_yield": "100-500 ng"
        })

    # Step 2: Host depletion
    if path["host_depletion_method"] == "cpg_methylation":
        steps.append({
            "step": 2,
            "name": "CpG Methylation-based Host DNA Depletion",
            "kit": "NEBNext Microbiome DNA Enrichment Kit",
            "input": "50 μL cfDNA",
            "output": "45 μL depleted DNA",
            "depletion_efficiency": "60-90% host DNA removed",
            "viral_recovery": "80-95%"
        })
    elif path["host_depletion_method"] == "rrna_depletion":
        steps.append({
            "step": 2,
            "name": "rRNA Depletion (RNase H method)",
            "kit": "NEBNext rRNA Depletion Kit (Human/Mouse/Rat)",
            "input": "50 μL cfRNA (after DNase treatment)",
            "output": "40 μL rRNA-depleted RNA",
            "depletion_efficiency": "80-95% rRNA removed",
            "viral_recovery": "70-90%",
            "note": "Partial cross-reactivity with pig rRNA"
        })
    elif path["host_depletion_method"] == "polya_selection":
        steps.append({
            "step": 2,
            "name": "Poly(A) Selection",
            "kit": "NEBNext Poly(A) mRNA Magnetic Isolation Module",
            "input": "50 μL cfRNA (after DNase treatment)",
            "output": "45 μL poly(A)+ RNA",
            "depletion_efficiency": "98% rRNA removed",
            "viral_recovery": "80-90% (poly(A)+ viruses only)",
            "advantage": "EEEV has poly(A) tail → high enrichment"
        })
    elif path["host_depletion_method"] == "none":
        steps.append({
            "step": 2,
            "name": "No Host Depletion",
            "reason": "Proviral DNA integrated into host genome",
            "alternative": "Nested PCR for specific detection"
        })

    # Step 3: Library preparation
    if path["library_prep_method"] == "dna_ligation":
        steps.append({
            "step": 3,
            "name": "DNA Ligation Library Prep",
            "kit": "Oxford Nanopore Ligation Sequencing Kit V14 (SQK-LSK114)",
            "input": "45 μL DNA",
            "output": "15 μL library (50-200 fmol)",
            "sequencing_time": "24-48 hours"
        })
    elif path["library_prep_method"] == "direct_rna":
        steps.append({
            "step": 3,
            "name": "Direct RNA Sequencing Library Prep",
            "kit": "Oxford Nanopore Direct RNA Sequencing Kit (SQK-RNA002)",
            "input": "45 μL poly(A)+ RNA (≥50 ng)",
            "output": "21 μL library",
            "sequencing_time": "24-48 hours",
            "advantage": "Native strand, RNA modifications detectable"
        })
    elif path["library_prep_method"] == "amplicon_rt_pcr":
        steps.append({
            "step": 3,
            "name": "Amplicon RT-PCR (Tiled, ARTIC-style)",
            "description": "Multiplex RT-PCR for Hantavirus L/M/S segments",
            "primers": "36 amplicons (400 bp, 100 bp overlap)",
            "input": "10 μL rRNA-depleted RNA",
            "output": "20 μL pooled amplicons",
            "sensitivity": "10-100 copies/mL",
            "kit_for_sequencing": "ONT Native Barcoding Kit (SQK-NBD114)"
        })
    elif path["library_prep_method"] == "nested_pcr":
        steps.append({
            "step": 3,
            "name": "Nested PCR (Spumavirus pol gene)",
            "description": "Two-round PCR with degenerate primers",
            "target": "pol gene (reverse transcriptase domain)",
            "input": "100 ng genomic DNA",
            "output": "~400 bp PCR product",
            "cycles": "1st: 35 cycles, 2nd: 25 cycles",
            "sensitivity": "1-10 copies/10⁵ PBMCs",
            "confirmation": "Sanger sequencing recommended"
        })

    return steps

def main():
    parser = argparse.ArgumentParser(
        description='Phase 0: Sample Preparation Router for PMDA 4-Virus Detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Route sample for polyomavirus and EEEV detection
  python sample_router.py --targets polyomavirus eeev --metadata sample_001.json --output workflows/

  # Generate workflow for all 4 PMDA viruses
  python sample_router.py --targets polyomavirus hantavirus eeev spumavirus \\
      --metadata sample_001.json --output workflows/
        """
    )

    parser.add_argument('--targets', nargs='+', required=True,
                       choices=['polyomavirus', 'hantavirus', 'eeev', 'spumavirus'],
                       help='Target virus(es) for detection')
    parser.add_argument('--metadata', type=Path, required=True,
                       help='Sample metadata JSON file')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output directory for workflow files')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate input, do not generate workflows')

    args = parser.parse_args()

    # Load sample metadata
    try:
        with open(args.metadata) as f:
            sample_metadata = json.load(f)
    except Exception as e:
        print(f"ERROR: Failed to load metadata file: {e}", file=sys.stderr)
        sys.exit(1)

    # Determine workflow
    print(f"Determining workflow for targets: {', '.join(args.targets)}...")
    workflow = determine_extraction_workflow(args.targets)

    if "error" in workflow:
        print(f"ERROR: {workflow['error']}", file=sys.stderr)
        print(f"Valid targets: {', '.join(workflow['valid_targets'])}", file=sys.stderr)
        sys.exit(1)

    # Validate requirements
    print("Validating sample requirements...")
    validation_errors = validate_input_requirements(workflow, sample_metadata)

    if validation_errors:
        print("VALIDATION ERRORS:", file=sys.stderr)
        for error in validation_errors:
            print(f"  - {error}", file=sys.stderr)
        if args.validate_only:
            sys.exit(1)
        else:
            print("\n[WARNING] Validation errors found but proceeding with workflow generation...\n")
    else:
        print("✓ Validation passed")

    if args.validate_only:
        print("Validation complete.")
        sys.exit(0)

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Write workflow summary
    workflow_file = args.output / "workflow_summary.json"
    with open(workflow_file, 'w') as f:
        json.dump(workflow, f, indent=2)
    print(f"✓ Workflow summary written to: {workflow_file}")

    # Generate detailed protocols
    print("Generating detailed extraction protocols...")
    protocol_files = generate_extraction_protocol(workflow, args.output)

    for virus, protocol_file in protocol_files.items():
        print(f"  ✓ {virus}: {protocol_file}")

    # Print workflow summary
    print("\n" + "="*70)
    print("WORKFLOW SUMMARY")
    print("="*70)
    print(f"Target viruses: {', '.join(args.targets)}")
    print(f"Blood volume required: {workflow['blood_collection']['volume_ml']} mL")
    print(f"RNase inhibitor: {'YES' if workflow['blood_collection']['rnase_inhibitor'] else 'NO'}")
    print(f"Plasma extraction: {workflow['extraction']['plasma_extraction']}")
    print(f"PBMC extraction: {workflow['extraction']['pbmc_extraction'] or 'None'}")

    print("\nProcessing Paths:")
    for i, path in enumerate(workflow['processing_paths'], 1):
        print(f"\n{i}. {path['virus'].upper()}")
        print(f"   Genome type: {path['genome_type']}")
        print(f"   Sample type: {path['sample_type']}")
        print(f"   Host depletion: {path['host_depletion_method']}")
        print(f"   Library prep: {path['library_prep_method']}")
        print(f"   Estimated LOD: {path['estimated_lod_copies_ml']}")

    print("\n" + "="*70)
    print(f"All workflow files written to: {args.output}")
    print("="*70)

if __name__ == '__main__':
    main()
