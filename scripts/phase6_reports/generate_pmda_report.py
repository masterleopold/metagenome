#!/usr/bin/env python3
"""
MinION Metagenomics Pipeline - PMDA Compliance Report Generator
Generates comprehensive reports for PMDA 91 pathogen screening
"""

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd
from jinja2 import Template
import boto3
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

logger = logging.getLogger(__name__)


class PMDAReportGenerator:
    """Generate PMDA-compliant pathogen screening reports."""

    # PMDA 91 pathogens list (abbreviated for example)
    PMDA_PATHOGENS = {
        'PERV-A': {'name_ja': 'ブタ内在性レトロウイルスA', 'risk': 'Critical'},
        'PERV-B': {'name_ja': 'ブタ内在性レトロウイルスB', 'risk': 'Critical'},
        'PERV-C': {'name_ja': 'ブタ内在性レトロウイルスC', 'risk': 'Critical'},
        'HEV': {'name_ja': 'E型肝炎ウイルス', 'risk': 'High'},
        'JEV': {'name_ja': '日本脳炎ウイルス', 'risk': 'High'},
        'SS': {'name_ja': 'ブタレンサ球菌', 'risk': 'High'},
        # ... (complete list of 91 pathogens)
    }

    def __init__(self, run_id: str, output_dir: Path):
        self.run_id = run_id
        self.output_dir = output_dir
        self.results = {}

    def load_results(self, results_dir: Path):
        """Load all analysis results."""
        # Load QC results
        qc_file = results_dir / 'qc' / 'qc_summary.json'
        if qc_file.exists():
            with open(qc_file) as f:
                self.results['qc'] = json.load(f)

        # Load Kraken2 results
        kraken_file = results_dir / 'kraken2' / 'pmda_pathogens.json'
        if kraken_file.exists():
            with open(kraken_file) as f:
                self.results['kraken2'] = json.load(f)

        # Load PERV results
        perv_file = results_dir / 'perv' / 'perv_summary.json'
        if perv_file.exists():
            with open(perv_file) as f:
                self.results['perv'] = json.load(f)

        # Load quantification results
        quant_file = results_dir / 'quantification' / 'absolute_quantification.json'
        if quant_file.exists():
            with open(quant_file) as f:
                self.results['quantification'] = json.load(f)

    def generate_pdf_report(self):
        """Generate PDF report for PMDA submission."""
        pdf_path = self.output_dir / f'{self.run_id}_pmda_report.pdf'
        doc = SimpleDocTemplate(str(pdf_path), pagesize=A4)
        story = []
        styles = getSampleStyleSheet()

        # Title
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#003366'),
            alignment=1  # Center
        )

        story.append(Paragraph('病原体スクリーニング報告書', title_style))
        story.append(Paragraph('Pathogen Screening Report', styles['Heading2']))
        story.append(Spacer(1, 0.5*inch))

        # Report metadata
        metadata = [
            ['検査ID / Run ID:', self.run_id],
            ['検査日 / Test Date:', datetime.now().strftime('%Y-%m-%d')],
            ['検査方法 / Method:', 'MinION Metagenomics (Duplex)'],
            ['準拠基準 / Compliance:', 'PMDA Guideline for Xenotransplantation']
        ]

        metadata_table = Table(metadata, colWidths=[2*inch, 4*inch])
        metadata_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
        ]))

        story.append(metadata_table)
        story.append(Spacer(1, 0.3*inch))

        # QC Summary
        story.append(Paragraph('品質管理結果 / Quality Control Results', styles['Heading2']))

        if 'qc' in self.results:
            qc_data = [
                ['指標 / Metric', '値 / Value', '基準 / Threshold', '判定 / Pass'],
                ['総リード数 / Total Reads', f"{self.results['qc']['metrics']['total_reads']:,}", '≥100,000', 'PASS' if self.results['qc']['metrics']['total_reads'] >= 100000 else 'FAIL'],
                ['平均品質スコア / Mean Q-score', f"{self.results['qc']['metrics']['mean_qscore']:.1f}", '≥9.0', 'PASS' if self.results['qc']['metrics']['mean_qscore'] >= 9 else 'FAIL'],
                ['N50', f"{self.results['qc']['metrics']['n50']:,}", '≥200', 'PASS' if self.results['qc']['metrics']['n50'] >= 200 else 'FAIL'],
            ]

            qc_table = Table(qc_data, colWidths=[2*inch, 1.5*inch, 1.5*inch, 0.5*inch])
            qc_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ]))

            story.append(qc_table)

        story.append(PageBreak())

        # Pathogen Detection Results
        story.append(Paragraph('病原体検出結果 / Pathogen Detection Results', styles['Heading2']))

        # Create detection summary table
        detection_data = [['病原体 / Pathogen', '検出 / Detected', 'リード数 / Reads', '定量値 / Copies/mL']]

        for pathogen_code, info in self.PMDA_PATHOGENS.items():
            detected = self._check_pathogen_detected(pathogen_code)
            reads = self._get_pathogen_reads(pathogen_code)
            copies = self._get_pathogen_copies(pathogen_code)

            row_data = [
                f"{info['name_ja']}\n{pathogen_code}",
                '検出' if detected else '不検出',
                str(reads) if reads > 0 else '-',
                f"{copies:.2e}" if copies > 0 else '-'
            ]

            detection_data.append(row_data)

        detection_table = Table(detection_data, colWidths=[2.5*inch, 1*inch, 1*inch, 1.5*inch])
        detection_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ]))

        story.append(detection_table)
        story.append(Spacer(1, 0.3*inch))

        # Critical alerts
        if self._has_critical_findings():
            story.append(Paragraph('[WARNING] 重要所見 / Critical Findings', styles['Heading3']))
            findings = self._get_critical_findings()
            for finding in findings:
                story.append(Paragraph(f"• {finding}", styles['BodyText']))

        # Conclusion
        story.append(Spacer(1, 0.5*inch))
        story.append(Paragraph('結論 / Conclusion', styles['Heading2']))

        conclusion = self._generate_conclusion()
        story.append(Paragraph(conclusion, styles['BodyText']))

        # Build PDF
        doc.build(story)
        logger.info(f"PDF report generated: {pdf_path}")

        return pdf_path

    def generate_json_report(self):
        """Generate machine-readable JSON report."""
        report = {
            'run_id': self.run_id,
            'report_timestamp': datetime.now().isoformat(),
            'compliance': 'PMDA',
            'qc_metrics': self.results.get('qc', {}),
            'pathogen_detections': [],
            'critical_findings': self._get_critical_findings(),
            'conclusion': self._generate_conclusion()
        }

        # Add pathogen detection results
        for pathogen_code, info in self.PMDA_PATHOGENS.items():
            detection = {
                'pathogen_code': pathogen_code,
                'pathogen_name_ja': info['name_ja'],
                'risk_level': info['risk'],
                'detected': self._check_pathogen_detected(pathogen_code),
                'read_count': self._get_pathogen_reads(pathogen_code),
                'copies_per_ml': self._get_pathogen_copies(pathogen_code)
            }
            report['pathogen_detections'].append(detection)

        json_path = self.output_dir / f'{self.run_id}_pmda_report.json'
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)

        logger.info(f"JSON report generated: {json_path}")
        return json_path

    def _check_pathogen_detected(self, pathogen_code: str) -> bool:
        """Check if pathogen was detected."""
        # Check in Kraken2 results
        if 'kraken2' in self.results:
            if pathogen_code in self.results['kraken2'].get('detected_pathogens', []):
                return True

        # Special check for PERV
        if pathogen_code.startswith('PERV') and 'perv' in self.results:
            return self.results['perv'].get('perv_detected', False)

        return False

    def _get_pathogen_reads(self, pathogen_code: str) -> int:
        """Get read count for pathogen."""
        if 'kraken2' in self.results:
            pathogens = self.results['kraken2'].get('pathogen_details', {})
            if pathogen_code in pathogens:
                return pathogens[pathogen_code].get('read_count', 0)
        return 0

    def _get_pathogen_copies(self, pathogen_code: str) -> float:
        """Get absolute quantification for pathogen."""
        if 'quantification' in self.results:
            quants = self.results['quantification'].get('pathogens', {})
            if pathogen_code in quants:
                return quants[pathogen_code].get('copies_per_ml', 0.0)
        return 0.0

    def _has_critical_findings(self) -> bool:
        """Check for critical findings requiring attention."""
        # Check for PERV
        if 'perv' in self.results and self.results['perv'].get('perv_detected'):
            return True

        # Check for high-risk pathogens
        for pathogen_code, info in self.PMDA_PATHOGENS.items():
            if info['risk'] in ['Critical', 'High'] and self._check_pathogen_detected(pathogen_code):
                return True

        return False

    def _get_critical_findings(self) -> List[str]:
        """Get list of critical findings."""
        findings = []

        if 'perv' in self.results and self.results['perv'].get('perv_detected'):
            findings.append('PERV (ブタ内在性レトロウイルス) が検出されました。異種移植の安全性評価が必要です。')

        for pathogen_code, info in self.PMDA_PATHOGENS.items():
            if info['risk'] in ['Critical', 'High'] and self._check_pathogen_detected(pathogen_code):
                findings.append(f"{info['name_ja']} ({pathogen_code}) が検出されました。")

        return findings

    def _generate_conclusion(self) -> str:
        """Generate report conclusion."""
        if not self._has_critical_findings():
            return "PMDA指定91種病原体のスクリーニングが完了しました。臨床的に重要な病原体は検出されませんでした。"
        else:
            return "PMDA指定91種病原体のスクリーニングが完了しました。重要な所見が認められたため、詳細な評価が必要です。"


def main():
    parser = argparse.ArgumentParser(description="Generate PMDA compliance report")
    parser.add_argument('-r', '--run-id', required=True, help='Run ID')
    parser.add_argument('-i', '--input-dir', required=True, type=Path, help='Results directory')
    parser.add_argument('-o', '--output-dir', required=True, type=Path, help='Output directory')
    parser.add_argument('--format', choices=['pdf', 'json', 'both'], default='both')

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate report
    generator = PMDAReportGenerator(args.run_id, args.output_dir)
    generator.load_results(args.input_dir)

    if args.format in ['pdf', 'both']:
        generator.generate_pdf_report()

    if args.format in ['json', 'both']:
        generator.generate_json_report()

    print(f"Reports generated in {args.output_dir}")


if __name__ == '__main__':
    main()