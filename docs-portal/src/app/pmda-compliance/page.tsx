import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { CodeBlock } from "@/components/ui/CodeBlock";
import {
  ShieldCheckIcon,
  AlertTriangleIcon,
  CheckCircleIcon,
  VibrateIcon,
  FileTextIcon,
} from "lucide-react";

const pathogenCategories = [
  { name: "Viruses", count: 41, detected: 5, color: "text-destructive dark:text-red-400", bgColor: "bg-red-50 dark:bg-destructive/10" },
  { name: "Bacteria", count: 27, detected: 2, color: "text-primary", bgColor: "bg-primary/10" },
  { name: "Parasites", count: 19, detected: 0, color: "text-primary", bgColor: "bg-primary/10" },
  { name: "Fungi", count: 2, detected: 0, color: "text-secondary", bgColor: "bg-secondary/10" },
  { name: "Special Mgmt", count: 5, detected: 0, color: "text-purple-600 dark:text-purple-400", bgColor: "bg-purple-50 dark:bg-purple-900/10" },
];

const criticalPathogens = [
  { code: "PERV-A", name: "Porcine Endogenous Retrovirus A", risk: "CRITICAL" },
  { code: "PERV-B", name: "Porcine Endogenous Retrovirus B", risk: "CRITICAL" },
  { code: "PERV-C", name: "Porcine Endogenous Retrovirus C", risk: "CRITICAL" },
  { code: "ASFV", name: "African Swine Fever Virus", risk: "CRITICAL" },
  { code: "CSFV", name: "Classical Swine Fever Virus", risk: "CRITICAL" },
  { code: "FMDV", name: "Foot-and-Mouth Disease Virus", risk: "CRITICAL" },
  { code: "PRION", name: "Transmissible Spongiform Encephalopathy", risk: "CRITICAL" },
];

const detectedExample = [
  { code: "HEV", name: "Hepatitis E Virus", reads: 450, confidence: 0.96, copies: 8.5e4 },
  { code: "PCV2", name: "Porcine Circovirus 2", reads: 1250, confidence: 0.98, copies: 5.2e4 },
  { code: "SS", name: "Streptococcus suis", reads: 320, confidence: 0.92, copies: 2.1e3 },
  { code: "EC", name: "Escherichia coli", reads: 180, confidence: 0.88, copies: 1.5e3 },
  { code: "CA", name: "Candida albicans", reads: 95, confidence: 0.85, copies: 8.2e2 },
];

export default function PMDACompliancePage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Compliance</Badge>
              <h1 className="text-4xl font-bold mb-4">PMDA Compliance</h1>
              <p className="text-lg text-muted-foreground">
                PMDA guidelines for xenotransplantation pathogen screening and the 91 designated pathogens.
              </p>
            </div>

            <Alert className="mb-8">
              <ShieldCheckIcon className="h-4 w-4" />
              <AlertTitle>Regulatory Framework</AlertTitle>
              <AlertDescription>
                This pipeline is designed to meet PMDA (Pharmaceuticals and Medical Devices Agency, Japan) guidelines for donor pig screening in xenotransplantation research.
              </AlertDescription>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">91 Pathogen Overview</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>PMDA Designated Pathogens</CardTitle>
                  <CardDescription>Complete coverage across 5 categories</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="grid md:grid-cols-5 gap-4 mb-6">
                    {pathogenCategories.map((category) => (
                      <div key={category.name} className={`${category.bgColor} p-4 rounded-lg text-center`}>
                        <div className={`text-3xl font-bold ${category.color}`}>{category.count}</div>
                        <div className="text-sm font-medium mt-1">{category.name}</div>
                      </div>
                    ))}
                  </div>

                  <div className="flex items-center justify-center gap-2">
                    <div className="text-4xl font-bold text-primary">91</div>
                    <div className="text-sm text-muted-foreground">Total Designated Pathogens</div>
                  </div>
                </CardContent>
              </Card>

              <div className="grid md:grid-cols-3 gap-4">
                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Detection Methods</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm">
                    <ul className="space-y-1">
                      <li>• Kraken2 classification</li>
                      <li>• BLAST sequence alignment</li>
                      <li>• Diamond protein search</li>
                      <li>• PERV-specific PCR primers</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Reporting Requirements</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm">
                    <ul className="space-y-1">
                      <li>• All 91 pathogens tested</li>
                      <li>• PERV-specific analysis</li>
                      <li>• Quantitative results</li>
                      <li>• Quality control metrics</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Validation Standards</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm">
                    <ul className="space-y-1">
                      <li>• Sensitivity: &gt;95%</li>
                      <li>• Specificity: &gt;98%</li>
                      <li>• LOD documented</li>
                      <li>• Reproducibility &gt;90%</li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Sample Preparation Protocol</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Protocol 12 v2.1: Unified Workflow (Recommended)</CardTitle>
                  <CardDescription>TRUE 100% coverage - includes circular and ssDNA viruses</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <Alert>
                    <ShieldCheckIcon className="h-4 w-4" />
                    <AlertTitle>Protocol 12 v2.1 Update (2025-11-13)</AlertTitle>
                    <AlertDescription>
                      Critical update adds Step 2.5 for circular/ssDNA virus support, achieving TRUE 91/91 pathogen coverage. Previous v1.0 could only detect 87/91 pathogens due to Oxford Nanopore LSK114 ligation kit limitations with circular DNA and single-stranded DNA.
                    </AlertDescription>
                  </Alert>

                  <div className="grid md:grid-cols-2 gap-4">
                    <div className="bg-muted p-4 rounded-lg">
                      <h4 className="font-semibold mb-2">Workflow Specifications</h4>
                      <div className="text-sm space-y-1">
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Workflows:</span>
                          <span className="font-medium">2 (DNA + RNA)</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Hands-on time:</span>
                          <span className="font-medium">15.5 hours</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Cost/sample:</span>
                          <span className="font-medium">¥162,000</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Coverage:</span>
                          <span className="font-medium text-primary">91/91 (TRUE 100%)</span>
                        </div>
                      </div>
                    </div>

                    <div className="bg-muted p-4 rounded-lg">
                      <h4 className="font-semibold mb-2">Key Features</h4>
                      <ul className="text-sm space-y-1">
                        <li className="flex items-start gap-2">
                          <CheckCircleIcon className="h-4 w-4 text-primary mt-0.5 shrink-0" />
                          <span>Universal poly(A) selection for all RNA viruses</span>
                        </li>
                        <li className="flex items-start gap-2">
                          <CheckCircleIcon className="h-4 w-4 text-primary mt-0.5 shrink-0" />
                          <span>CpG methylation-based host depletion</span>
                        </li>
                        <li className="flex items-start gap-2">
                          <CheckCircleIcon className="h-4 w-4 text-primary mt-0.5 shrink-0" />
                          <span>Step 2.5: Circular DNA linearization + ssDNA conversion</span>
                        </li>
                        <li className="flex items-start gap-2">
                          <CheckCircleIcon className="h-4 w-4 text-primary mt-0.5 shrink-0" />
                          <span>LOD: 100-500 copies/mL (screening-sufficient)</span>
                        </li>
                      </ul>
                    </div>
                  </div>

                  <div className="bg-purple-50 dark:bg-purple-900/10 p-4 rounded-lg border border-purple-200 dark:border-purple-800">
                    <h4 className="font-semibold mb-2 text-sm">Step 2.5: Circular & ssDNA Virus Support</h4>
                    <p className="text-sm text-muted-foreground mb-2">
                      Added to enable detection of 4 critical PMDA pathogens that were incompatible with standard ligation-based library preparation:
                    </p>
                    <div className="grid md:grid-cols-2 gap-3 text-sm">
                      <div>
                        <div className="font-semibold mb-1">Affected Pathogens:</div>
                        <ul className="text-xs space-y-0.5">
                          <li>• <strong>PCV2</strong> - Porcine Circovirus 2 (circular ssDNA, 1.7kb) - Special Mgmt</li>
                          <li>• <strong>PCV3</strong> - Porcine Circovirus 3 (circular ssDNA, 2.0kb) - Special Mgmt</li>
                          <li>• <strong>TTV</strong> - Torque teno virus (circular ssDNA, 3.8kb)</li>
                          <li>• <strong>PPV</strong> - Porcine parvovirus (linear ssDNA, 5.0kb)</li>
                        </ul>
                      </div>
                      <div>
                        <div className="font-semibold mb-1">Process Steps:</div>
                        <ul className="text-xs space-y-0.5">
                          <li>• <strong>Sub-step 2.5.1:</strong> DNase I linearization (0.005 U, 30 min)</li>
                          <li>• <strong>Sub-step 2.5.2:</strong> Klenow Fragment second-strand synthesis (2 hours)</li>
                          <li>• <strong>Time added:</strong> +2.5 hours (+19%)</li>
                          <li>• <strong>Cost added:</strong> +¥5,000 (+3.2%)</li>
                        </ul>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Conditional Screening Strategy</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Two-Tier Detection Approach</CardTitle>
                  <CardDescription>Universal screening + conditional follow-up for low-probability pathogens</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    While most PMDA pathogens are screened universally (100% of samples), certain low-probability pathogens trigger conditional follow-up testing only when initial metagenomic screening detects relevant signatures.
                  </p>

                  <div className="grid md:grid-cols-2 gap-4">
                    <div className="bg-primary/5 dark:bg-primary/10 p-4 rounded-lg border border-primary/20">
                      <h4 className="font-semibold mb-2 flex items-center gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary" />
                        Universal Screening (90 pathogens)
                      </h4>
                      <p className="text-sm text-muted-foreground mb-2">
                        All samples undergo metagenomic analysis via Protocol 12
                      </p>
                      <ul className="text-sm space-y-1">
                        <li>• PERV (A/B/C subtypes)</li>
                        <li>• All 41 DNA/RNA viruses</li>
                        <li>• All 27 bacterial pathogens</li>
                        <li>• All 19 parasites</li>
                        <li>• All 2 fungi</li>
                      </ul>
                    </div>

                    <div className="bg-purple-50 dark:bg-purple-900/10 p-4 rounded-lg border border-purple-200 dark:border-purple-800">
                      <h4 className="font-semibold mb-2 flex items-center gap-2">
                        <AlertTriangleIcon className="h-5 w-5 text-purple-600 dark:text-purple-400" />
                        Conditional Testing (1 pathogen)
                      </h4>
                      <p className="text-sm text-muted-foreground mb-2">
                        Triggered only when metagenomic data shows relevant signatures
                      </p>
                      <div className="bg-white dark:bg-purple-950/30 p-3 rounded">
                        <div className="font-medium text-sm mb-1">Porcine Spumavirus</div>
                        <div className="text-xs text-muted-foreground mb-2">
                          Protocol 13: Nested PCR + Phylogenetic analysis
                        </div>
                        <div className="text-xs space-y-1">
                          <div><strong>Trigger rate:</strong> 5-10% of samples</div>
                          <div><strong>Detection probability:</strong> &lt;0.1%</div>
                          <div><strong>Reason:</strong> 0 detections in 70 years, tested precautionarily</div>
                        </div>
                      </div>
                    </div>
                  </div>

                  <Alert>
                    <AlertDescription>
                      <strong>Precautionary Testing:</strong> Spumavirus has never been detected in pigs despite 70 years of virology research. It remains on the PMDA list based on the precautionary principle: phylogenetic proximity to foamy virus-positive species (cattle, horses) and potential for cross-species transmission justifies continued surveillance despite zero historical detections.
                    </AlertDescription>
                  </Alert>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Critical Pathogens</h2>

              <Alert variant="destructive" className="mb-6">
                <AlertTriangleIcon className="h-4 w-4" />
                <AlertTitle>Immediate Alert Pathogens</AlertTitle>
                <AlertDescription>
                  Detection of these pathogens triggers immediate SNS notification and quarantine procedures.
                </AlertDescription>
              </Alert>

              <Card>
                <CardHeader>
                  <CardTitle>High-Priority Screening</CardTitle>
                  <CardDescription>7 pathogens requiring immediate action</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-3">
                    {criticalPathogens.map((pathogen) => (
                      <div key={pathogen.code} className="flex items-center justify-between p-3 bg-red-50 dark:bg-red-900/10 rounded-lg border border-red-200 dark:border-red-800">
                        <div className="flex items-center gap-3">
                          {pathogen.code.startsWith('PERV') && (
                            <VibrateIcon className="h-5 w-5 text-red-700 dark:text-red-400" />
                          )}
                          <div>
                            <div className="font-semibold">{pathogen.code}</div>
                            <div className="text-sm text-muted-foreground">{pathogen.name}</div>
                          </div>
                        </div>
                        <Badge variant="error">{pathogen.risk}</Badge>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">PERV Detection</h2>

              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-3">
                    <VibrateIcon className="h-8 w-8 text-red-600" />
                    <div>
                      <CardTitle>Porcine Endogenous Retrovirus (PERV)</CardTitle>
                      <CardDescription>Most critical pathogen for xenotransplantation</CardDescription>
                    </div>
                  </div>
                </CardHeader>
                <CardContent className="space-y-4">
                  <Alert variant="destructive">
                    <AlertDescription>
                      PERV is the highest priority pathogen. Any detection of PERV-A, PERV-B, or PERV-C results in immediate rejection of the donor pig and triggering of quarantine protocols.
                    </AlertDescription>
                  </Alert>

                  <div>
                    <h4 className="font-semibold mb-2">Detection Strategy</h4>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Targeted primers:</strong> PERV-specific PCR primers for env gene
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>NGS screening:</strong> Metagenomic reads aligned to PERV reference sequences
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Typing:</strong> Differentiation between PERV-A, PERV-B, and PERV-C subtypes
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Quantification:</strong> Copy number per mL plasma with 95% confidence intervals
                        </div>
                      </li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Action Thresholds</h4>
                    <div className="bg-red-50 dark:bg-red-900/10 p-4 rounded-lg border border-red-200 dark:border-red-800">
                      <div className="grid grid-cols-2 gap-4 text-sm">
                        <div>
                          <div className="text-muted-foreground mb-1">Detection Threshold</div>
                          <div className="font-bold">1 read</div>
                        </div>
                        <div>
                          <div className="text-muted-foreground mb-1">Action</div>
                          <div className="font-bold">Immediate Quarantine</div>
                        </div>
                        <div>
                          <div className="text-muted-foreground mb-1">Notification</div>
                          <div className="font-bold">Immediate SNS Alert</div>
                        </div>
                        <div>
                          <div className="text-muted-foreground mb-1">Follow-up</div>
                          <div className="font-bold">PCR Confirmation</div>
                        </div>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Example Detection Report</h2>

              <Card>
                <CardHeader>
                  <FileTextIcon className="h-8 w-8 text-primary mb-2" />
                  <CardTitle>Sample PMDA Report</CardTitle>
                  <CardDescription>Run ID: RUN-2024-001 | Date: 2024-01-15</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div className="bg-muted p-4 rounded-lg">
                    <h4 className="font-semibold mb-3">Detection Summary</h4>
                    <div className="grid grid-cols-2 gap-4 text-sm">
                      <div>
                        <div className="text-muted-foreground mb-1">Total PMDA Pathogens</div>
                        <div className="text-2xl font-bold">91</div>
                      </div>
                      <div>
                        <div className="text-muted-foreground mb-1">Detected</div>
                        <div className="text-2xl font-bold text-secondary">5</div>
                      </div>
                      <div>
                        <div className="text-muted-foreground mb-1">Not Detected</div>
                        <div className="text-2xl font-bold text-primary">86</div>
                      </div>
                      <div>
                        <div className="text-muted-foreground mb-1">Critical Detected</div>
                        <div className="text-2xl font-bold text-primary">0</div>
                      </div>
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-3">PERV Analysis</h4>
                    <Alert variant="success">
                      <CheckCircleIcon className="h-4 w-4" />
                      <AlertTitle>PERV Not Detected</AlertTitle>
                      <AlertDescription>
                        No PERV-A, PERV-B, or PERV-C sequences detected. Donor pig passes PERV screening.
                      </AlertDescription>
                    </Alert>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-3">Detected Pathogens (Non-Critical)</h4>
                    <div className="space-y-2">
                      {detectedExample.map((pathogen) => (
                        <div key={pathogen.code} className="flex items-center justify-between p-3 bg-secondary/5 dark:bg-secondary/10 rounded border border-secondary/30">
                          <div>
                            <div className="font-semibold">{pathogen.code} - {pathogen.name}</div>
                            <div className="text-sm text-muted-foreground">
                              {pathogen.reads} reads | Confidence: {(pathogen.confidence * 100).toFixed(1)}% | {pathogen.copies.toExponential(1)} copies/mL
                            </div>
                          </div>
                          <Badge variant="warning">MEDIUM</Badge>
                        </div>
                      ))}
                    </div>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">JSON Report Format</h4>
                    <CodeBlock
                      code={`{
  "run_id": "RUN-2024-001",
  "analysis_date": "2024-01-15T18:00:00Z",
  "pmda_compliance": {
    "total_pathogens": 91,
    "detected": 5,
    "not_detected": 86,
    "critical_detected": 0,
    "detection_rate": 0.055
  },
  "perv_analysis": {
    "detected": false,
    "perv_a": false,
    "perv_b": false,
    "perv_c": false,
    "pass": true
  },
  "detected_pathogens": [
    {
      "code": "PCV2",
      "name": "Porcine circovirus 2",
      "category": "virus",
      "reads": 1250,
      "confidence": 0.98,
      "risk_level": "MEDIUM",
      "is_pmda": true,
      "quantification": {
        "copies_per_ml": 52000,
        "log10_copies": 4.72
      }
    }
  ],
  "quality_metrics": {
    "total_reads": 50000,
    "mean_quality": 10.5,
    "depletion_efficiency": 0.95
  },
  "compliance_status": "PASS"
}`}
                      language="json"
                      filename="pmda_report.json"
                    />
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Validation Requirements</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Analytical Validation</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Limit of Detection (LOD):</strong> Documented for each pathogen
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Reproducibility:</strong> CV &lt;20% across replicates
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Specificity:</strong> No cross-reactivity with host genome
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Linearity:</strong> R² &gt; 0.90 for quantification
                        </div>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Quality Control</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Positive Control:</strong> Spiked pathogen standards
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Negative Control:</strong> Pathogen-free plasma
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Spike-in Control:</strong> PhiX174 for normalization
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <div>
                          <strong>Data Integrity:</strong> ALCOA+ principles
                        </div>
                      </li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-6">Regulatory Resources</h2>

              <Card>
                <CardHeader>
                  <CardTitle>PMDA Guidelines</CardTitle>
                </CardHeader>
                <CardContent>
                  <ul className="text-sm space-y-2">
                    <li>
                      <strong>Xenotransplantation Guidelines (2024):</strong>
                      <div className="text-muted-foreground">
                        Comprehensive guidelines for donor animal screening and pathogen testing
                      </div>
                    </li>
                    <li>
                      <strong>91 Designated Pathogens List:</strong>
                      <div className="text-muted-foreground">
                        Complete list of pathogens requiring screening for xenotransplantation
                      </div>
                    </li>
                    <li>
                      <strong>NGS Validation Guidance:</strong>
                      <div className="text-muted-foreground">
                        Technical requirements for NGS-based pathogen detection methods
                      </div>
                    </li>
                  </ul>
                </CardContent>
              </Card>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
