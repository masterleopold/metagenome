import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { CodeBlock } from "@/components/ui/CodeBlock";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import {
  ZapIcon,
  CheckCircleIcon,
  DatabaseIcon,
  SearchIcon,
  CalculatorIcon,
  FileTextIcon,
  ClockIcon,
  ServerIcon,
  InfoIcon,
} from "lucide-react";

const phases = [
  {
    id: 1,
    name: "Basecalling",
    icon: ZapIcon,
    description: "Convert raw signal (FAST5/POD5) to sequences (FASTQ)",
    tool: "Dorado Duplex",
    duration: "2-4 hours",
    instance: "g4dn.xlarge (GPU)",
    color: "text-yellow-600",
    bgColor: "bg-yellow-50 dark:bg-yellow-900/10",
    borderColor: "border-yellow-200 dark:border-yellow-800",
  },
  {
    id: 2,
    name: "Quality Control",
    icon: CheckCircleIcon,
    description: "Assess read quality metrics and filter low-quality reads",
    tool: "PycoQC, NanoPlot",
    duration: "10-15 minutes",
    instance: "t3.large",
    color: "text-green-600",
    bgColor: "bg-green-50 dark:bg-green-900/10",
    borderColor: "border-green-200 dark:border-green-800",
  },
  {
    id: 3,
    name: "Host Removal",
    icon: DatabaseIcon,
    description: "Align reads to Sus scrofa genome and remove host DNA",
    tool: "Minimap2, SAMtools",
    duration: "30-60 minutes",
    instance: "r5.xlarge",
    color: "text-blue-600",
    bgColor: "bg-blue-50 dark:bg-blue-900/10",
    borderColor: "border-blue-200 dark:border-blue-800",
  },
  {
    id: 4,
    name: "Pathogen Detection",
    icon: SearchIcon,
    description: "Multi-database screening for PMDA pathogens",
    tool: "Kraken2, BLAST, Diamond",
    duration: "1-2 hours",
    instance: "r5.4xlarge",
    color: "text-red-600",
    bgColor: "bg-red-50 dark:bg-red-900/10",
    borderColor: "border-red-200 dark:border-red-800",
  },
  {
    id: 5,
    name: "Quantification",
    icon: CalculatorIcon,
    description: "Absolute copy number calculation with spike-in normalization",
    tool: "Custom Python scripts",
    duration: "15-30 minutes",
    instance: "t3.large",
    color: "text-purple-600",
    bgColor: "bg-purple-50 dark:bg-purple-900/10",
    borderColor: "border-purple-200 dark:border-purple-800",
  },
  {
    id: 6,
    name: "Reporting",
    icon: FileTextIcon,
    description: "Generate PMDA-compliant reports in PDF, JSON, and HTML",
    tool: "ReportLab, WeasyPrint",
    duration: "10-15 minutes",
    instance: "t3.medium",
    color: "text-teal-600",
    bgColor: "bg-teal-50 dark:bg-teal-900/10",
    borderColor: "border-teal-200 dark:border-teal-800",
  },
];

export default function PipelinePhasesPage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Core Concepts</Badge>
              <h1 className="text-4xl font-bold mb-4">Pipeline Phases</h1>
              <p className="text-lg text-muted-foreground">
                Deep dive into each of the 6 analysis phases from basecalling to reporting.
              </p>
            </div>

            <Alert className="mb-8">
              <InfoIcon className="h-4 w-4" />
              <AlertTitle>Total Pipeline Duration</AlertTitle>
              <AlertDescription>
                Complete analysis typically takes 4-8 hours depending on sample size and complexity. Each phase runs sequentially with automated transitions.
              </AlertDescription>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase Overview</h2>

              <div className="space-y-6">
                {phases.map((phase) => {
                  const Icon = phase.icon;
                  return (
                    <Card key={phase.id} className={`${phase.borderColor} border-2`}>
                      <CardHeader>
                        <div className="flex items-start gap-4">
                          <div className={`${phase.bgColor} p-3 rounded-lg`}>
                            <Icon className={`h-8 w-8 ${phase.color}`} />
                          </div>
                          <div className="flex-1">
                            <div className="flex items-center gap-2 mb-2">
                              <Badge>{`Phase ${phase.id}`}</Badge>
                              <CardTitle>{phase.name}</CardTitle>
                            </div>
                            <CardDescription className="text-base">
                              {phase.description}
                            </CardDescription>
                          </div>
                        </div>
                      </CardHeader>
                      <CardContent>
                        <div className="grid grid-cols-3 gap-4 text-sm">
                          <div>
                            <div className="text-muted-foreground mb-1">Tool</div>
                            <div className="font-medium">{phase.tool}</div>
                          </div>
                          <div>
                            <div className="text-muted-foreground mb-1">Duration</div>
                            <div className="font-medium flex items-center gap-1">
                              <ClockIcon className="h-4 w-4" />
                              {phase.duration}
                            </div>
                          </div>
                          <div>
                            <div className="text-muted-foreground mb-1">Instance Type</div>
                            <div className="font-medium flex items-center gap-1">
                              <ServerIcon className="h-4 w-4" />
                              {phase.instance}
                            </div>
                          </div>
                        </div>
                      </CardContent>
                    </Card>
                  );
                })}
              </div>
            </section>

            {/* Phase 1: Basecalling */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 1: Basecalling</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Dorado Duplex Mode</CardTitle>
                  <CardDescription>Q30+ accuracy basecalling</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Converts raw electrical signals from FAST5/POD5 files into nucleotide sequences (FASTQ) using Oxford Nanopore's Dorado basecaller in duplex mode for maximum accuracy.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Key Features</h4>
                    <ul className="text-sm space-y-1 text-muted-foreground">
                      <li>• <strong>Duplex Mode:</strong> Sequences both DNA strands for 99.9% accuracy (Q30+)</li>
                      <li>• <strong>GPU Acceleration:</strong> NVIDIA T4 GPU on g4dn.xlarge instance</li>
                      <li>• <strong>Real-time Processing:</strong> Can process data as it's generated</li>
                      <li>• <strong>Quality Filtering:</strong> Automatically filters reads below Q9</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Script Example</h4>
                    <CodeBlock
                      code={`#!/usr/bin/env bash
# Basecalling with Dorado Duplex

DORADO_BIN=/opt/dorado/bin/dorado
MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg
INPUT_DIR=/data/fast5
OUTPUT_DIR=/data/fastq

# Run duplex basecalling
$DORADO_BIN duplex \\
  --device cuda:0 \\
  $MODEL \\
  $INPUT_DIR > $OUTPUT_DIR/basecalled.fastq

# Generate sequencing summary
python3 generate_summary.py \\
  --input $OUTPUT_DIR/basecalled.fastq \\
  --output $OUTPUT_DIR/sequencing_summary.txt`}
                      language="bash"
                      filename="basecall_duplex.sh"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Output Metrics</h4>
                    <div className="grid grid-cols-2 gap-4 text-sm">
                      <div className="bg-muted p-3 rounded">
                        <div className="text-muted-foreground mb-1">Total Reads</div>
                        <div className="text-2xl font-bold">50,000+</div>
                      </div>
                      <div className="bg-muted p-3 rounded">
                        <div className="text-muted-foreground mb-1">Mean Quality</div>
                        <div className="text-2xl font-bold">Q10.5</div>
                      </div>
                      <div className="bg-muted p-3 rounded">
                        <div className="text-muted-foreground mb-1">Total Bases</div>
                        <div className="text-2xl font-bold">150 Mb</div>
                      </div>
                      <div className="bg-muted p-3 rounded">
                        <div className="text-muted-foreground mb-1">N50 Length</div>
                        <div className="text-2xl font-bold">3.5 kb</div>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Phase 2: QC */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 2: Quality Control</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>PycoQC & NanoPlot</CardTitle>
                  <CardDescription>Read quality assessment</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Comprehensive quality control analysis to ensure data meets minimum standards before proceeding to analysis phases.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Quality Thresholds</h4>
                    <ul className="text-sm space-y-1">
                      <li>✓ Minimum reads: <strong>10,000</strong></li>
                      <li>✓ Mean quality score: <strong>Q9+</strong></li>
                      <li>✓ Q30 reads: <strong>&gt;10%</strong></li>
                      <li>✓ Read length N50: <strong>&gt;2 kb</strong></li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">QC Reports Generated</h4>
                    <CodeBlock
                      code={`# PycoQC HTML report
pycoQC -f sequencing_summary.txt -o pycoQC_report.html

# NanoPlot visualization
NanoPlot --fastq basecalled.fastq --plots kde --legacy hex dot

# Output files:
# - NanoPlot-report.html
# - Read length distribution
# - Quality score distribution
# - Yield over time`}
                      language="bash"
                    />
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Phase 3: Host Removal */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 3: Host Genome Removal</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Minimap2 Alignment</CardTitle>
                  <CardDescription>Sus scrofa genome depletion</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Aligns reads to the porcine reference genome (Sus scrofa 11.1) and removes host DNA to enrich for pathogen sequences.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Alignment & Filtering</h4>
                    <CodeBlock
                      code={`# Align to host genome
minimap2 -ax map-ont \\
  /data/references/sus_scrofa_11.1.mmi \\
  basecalled.fastq > aligned.sam

# Extract unmapped reads (potential pathogens)
samtools view -f 4 aligned.sam | \\
  samtools fastq - > unmapped.fastq

# Calculate depletion statistics
python3 calculate_depletion_stats.py \\
  --total $(wc -l < basecalled.fastq) \\
  --unmapped $(wc -l < unmapped.fastq)`}
                      language="bash"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Expected Depletion</h4>
                    <Alert variant="success">
                      <AlertDescription>
                        Typical depletion efficiency: <strong>90-99%</strong> of reads should map to host genome for blood samples.
                        Remaining <strong>1-10%</strong> unmapped reads proceed to pathogen detection.
                      </AlertDescription>
                    </Alert>
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Phase 4: Pathogen Detection */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 4: Pathogen Detection</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Multi-Database Search</CardTitle>
                  <CardDescription>Kraken2, BLAST, Diamond, PERV-specific</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Comprehensive pathogen screening using multiple complementary databases and detection methods.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Detection Pipeline</h4>
                    <CodeBlock
                      code={`# 1. Kraken2 classification (rapid screening)
kraken2 --db /data/kraken2_db \\
  --threads 16 \\
  --report kraken_report.txt \\
  unmapped.fastq > kraken_output.txt

# 2. BLAST search against PMDA database
blastn -query unmapped.fastq \\
  -db /data/pmda_pathogens \\
  -num_threads 16 \\
  -outfmt 6 -out blast_results.txt

# 3. Diamond viral protein search
diamond blastx \\
  --query unmapped.fastq \\
  --db /data/rvdb.dmnd \\
  --threads 16 \\
  --outfmt 6 -out diamond_results.txt

# 4. PERV-specific detection
bash perv_analysis.sh \\
  --input unmapped.fastq \\
  --output perv_results.json`}
                      language="bash"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">PERV Critical Detection</h4>
                    <Alert variant="error">
                      <AlertTitle>Critical Alert Trigger</AlertTitle>
                      <AlertDescription>
                        Any PERV detection (PERV-A, PERV-B, PERV-C) triggers immediate SNS notification to alert recipients. This is the highest priority pathogen for xenotransplantation safety.
                      </AlertDescription>
                    </Alert>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Result Integration</h4>
                    <CodeBlock
                      code={`# Integrate results from all methods
python3 integrate_results.py \\
  --kraken kraken_output.txt \\
  --blast blast_results.txt \\
  --diamond diamond_results.txt \\
  --perv perv_results.json \\
  --output integrated_pathogens.json

# Validate against PMDA checklist
python3 pmda_check.py \\
  --results integrated_pathogens.json \\
  --checklist /data/pmda_91_pathogens.json \\
  --output pmda_compliance.json`}
                      language="python"
                    />
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Phase 5: Quantification */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 5: Quantification</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Spike-in Normalization</CardTitle>
                  <CardDescription>Absolute copy number calculation</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Converts read counts to absolute pathogen copy numbers using PhiX174 spike-in as internal standard.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Quantification Formula</h4>
                    <div className="bg-muted p-4 rounded-lg">
                      <code className="text-sm">
                        Pathogen copies/mL = (Pathogen reads / Spike-in reads) × Spike-in copies/mL
                      </code>
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Calculation Script</h4>
                    <CodeBlock
                      code={`# Absolute quantification with confidence intervals
python3 absolute_quantification.py \\
  --pathogen-reads 1250 \\
  --spike-in-reads 1000 \\
  --spike-in-copies 1000000 \\
  --confidence 0.95 \\
  --output quantification.json

# Output:
# {
#   "copies_per_ml": 1250000,
#   "log10_copies": 6.097,
#   "ci_lower": 1180000,
#   "ci_upper": 1320000,
#   "confidence_level": 0.95
# }`}
                      language="python"
                    />
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Phase 6: Reporting */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Phase 6: Report Generation</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>PMDA-Compliant Reports</CardTitle>
                  <CardDescription>PDF, JSON, HTML formats</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    Generates comprehensive analysis reports in multiple formats for different audiences.
                  </p>

                  <div>
                    <h4 className="font-semibold mb-2">Report Types</h4>
                    <div className="grid md:grid-cols-3 gap-4">
                      <div className="bg-muted p-3 rounded">
                        <div className="font-semibold mb-1">PDF Report</div>
                        <div className="text-xs text-muted-foreground">
                          Comprehensive report with visualizations for human review
                        </div>
                      </div>
                      <div className="bg-muted p-3 rounded">
                        <div className="font-semibold mb-1">JSON Report</div>
                        <div className="text-xs text-muted-foreground">
                          Machine-readable PMDA 91 pathogen checklist
                        </div>
                      </div>
                      <div className="bg-muted p-3 rounded">
                        <div className="font-semibold mb-1">HTML Report</div>
                        <div className="text-xs text-muted-foreground">
                          Interactive web-based report with search
                        </div>
                      </div>
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Report Generation</h4>
                    <CodeBlock
                      code={`# Generate all report formats
python3 generate_reports.py \\
  --run-id RUN-2024-001 \\
  --results integrated_pathogens.json \\
  --quantification quantification.json \\
  --pmda pmda_compliance.json \\
  --output-dir /data/reports

# Send notifications
python3 send_notifications.py \\
  --run-id RUN-2024-001 \\
  --reports /data/reports \\
  --recipients alerts@example.com`}
                      language="bash"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Report Contents</h4>
                    <ul className="text-sm space-y-1">
                      <li>✓ Executive summary</li>
                      <li>✓ PMDA 91 pathogen checklist</li>
                      <li>✓ Detailed pathogen detection results</li>
                      <li>✓ Quantification data (copies/mL)</li>
                      <li>✓ Quality control metrics</li>
                      <li>✓ PERV analysis section</li>
                      <li>✓ Pipeline execution log</li>
                    </ul>
                  </div>
                </CardContent>
              </Card>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
