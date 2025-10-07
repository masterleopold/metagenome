import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { ShieldCheckIcon, AlertTriangleIcon } from "lucide-react";

export default function PMDACompliancePage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Compliance</Badge>
              <h1 className="text-4xl font-bold mb-4">PMDA Guidelines</h1>
              <p className="text-lg text-muted-foreground">
                Regulatory compliance for xenotransplantation pathogen screening in Japan.
              </p>
            </div>

            <Alert className="mb-8" style={{
              backgroundColor: 'hsl(var(--primary) / 0.1)',
              borderColor: 'hsl(var(--primary) / 0.3)'
            }}>
              <ShieldCheckIcon className="h-4 w-4" style={{ color: 'hsl(var(--primary))' }} />
              <div>
                <AlertTitle>PMDA Compliant</AlertTitle>
                <AlertDescription>
                  This pipeline implements all requirements specified by the Pharmaceuticals and
                  Medical Devices Agency (PMDA) for xenotransplantation donor screening.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Regulatory Overview</h2>
              <p className="text-muted-foreground mb-4">
                The PMDA (Pharmaceuticals and Medical Devices Agency) and MHLW (Ministry of Health,
                Labour and Welfare) have established comprehensive guidelines for xenotransplantation
                safety in Japan. This pipeline addresses the requirement for detection of both
                "known and unknown infectious diseases" in donor animals.
              </p>

              <div className="border-l-2 border-primary pl-6 mb-6">
                <h3 className="text-lg font-semibold mb-2">Key Requirements</h3>
                <ul className="list-disc pl-6 space-y-2 text-muted-foreground">
                  <li>Detection of 91 PMDA-designated pathogens</li>
                  <li>Screening during minimum 3-week quarantine period</li>
                  <li>Detection of pathogens during latency periods</li>
                  <li>Donor animals from closed, documented environments</li>
                  <li>Comprehensive documentation and traceability</li>
                </ul>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">91 Designated Pathogens</h2>
              <p className="text-muted-foreground mb-4">
                The PMDA requires screening for 91 specific pathogens across multiple categories:
              </p>

              <div className="grid md:grid-cols-2 gap-4 mb-6">
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-3">Viruses</h3>
                  <ul className="text-sm text-muted-foreground space-y-1">
                    <li>• Porcine Endogenous Retroviruses (PERV-A, B, C)</li>
                    <li>• Porcine Cytomegalovirus (PCMV)</li>
                    <li>• Porcine Circovirus (PCV1, PCV2)</li>
                    <li>• Hepatitis E Virus (HEV)</li>
                    <li>• Japanese Encephalitis Virus</li>
                    <li>• And 40+ additional viral pathogens</li>
                  </ul>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-3">Bacteria</h3>
                  <ul className="text-sm text-muted-foreground space-y-1">
                    <li>• Mycoplasma species</li>
                    <li>• Leptospira interrogans</li>
                    <li>• Salmonella species</li>
                    <li>• Brucella suis</li>
                    <li>• Streptococcus suis</li>
                    <li>• And 30+ additional bacterial pathogens</li>
                  </ul>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-3">Parasites</h3>
                  <ul className="text-sm text-muted-foreground space-y-1">
                    <li>• Toxoplasma gondii</li>
                    <li>• Trichinella spiralis</li>
                    <li>• And additional parasitic pathogens</li>
                  </ul>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-3">Fungi</h3>
                  <ul className="text-sm text-muted-foreground space-y-1">
                    <li>• Candida species</li>
                    <li>• Aspergillus species</li>
                    <li>• And additional fungal pathogens</li>
                  </ul>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">PERV Detection (Critical)</h2>
              <Alert className="mb-6">
                <AlertTriangleIcon className="h-4 w-4" />
                <div>
                  <AlertTitle>High Priority</AlertTitle>
                  <AlertDescription>
                    PERV (Porcine Endogenous Retrovirus) detection is the highest priority for
                    PMDA compliance. All detections trigger immediate SNS alerts.
                  </AlertDescription>
                </div>
              </Alert>

              <div className="space-y-4">
                <div className="border-l-2 pl-6" style={{ borderColor: 'hsl(var(--destructive))' }}>
                  <h3 className="text-lg font-semibold mb-2">PERV-A</h3>
                  <p className="text-sm text-muted-foreground">
                    Most common subtype, capable of infecting human cells in vitro
                  </p>
                </div>
                <div className="border-l-2 pl-6" style={{ borderColor: 'hsl(var(--destructive))' }}>
                  <h3 className="text-lg font-semibold mb-2">PERV-B</h3>
                  <p className="text-sm text-muted-foreground">
                    Also capable of human cell infection, variant envelope genes
                  </p>
                </div>
                <div className="border-l-2 pl-6" style={{ borderColor: 'hsl(var(--destructive))' }}>
                  <h3 className="text-lg font-semibold mb-2">PERV-C</h3>
                  <p className="text-sm text-muted-foreground">
                    Pig-tropic only, but recombination risk with PERV-A
                  </p>
                </div>
              </div>

              <div className="mt-6 bg-muted p-4 rounded-lg">
                <h4 className="font-semibold mb-2">Pipeline PERV Detection Features:</h4>
                <ul className="text-sm text-muted-foreground space-y-1">
                  <li>• Subtype-specific identification (PERV-A/B/C)</li>
                  <li>• Recombinant detection (PERV-A/C hybrids)</li>
                  <li>• Phylogenetic analysis for strain identification</li>
                  <li>• Quantification in copies/mL plasma</li>
                  <li>• Immediate SNS notification on detection</li>
                </ul>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Analytical Validation</h2>
              <p className="text-muted-foreground mb-4">
                The pipeline implements comprehensive validation following FDA and BloodPAC guidance:
              </p>

              <div className="space-y-4">
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Limit of Detection (LOD)</h3>
                  <p className="text-sm text-muted-foreground mb-2">
                    Determined using serially diluted spiked samples
                  </p>
                  <div className="bg-muted p-3 rounded text-sm">
                    Target: {'<'}10 copies/mL for priority pathogens
                  </div>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Reproducibility</h3>
                  <p className="text-sm text-muted-foreground mb-2">
                    Tested across operators, instruments, and dates
                  </p>
                  <div className="bg-muted p-3 rounded text-sm">
                    Target: CV {'<'}20% within runs, CV {'<'}30% between runs
                  </div>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Specificity</h3>
                  <p className="text-sm text-muted-foreground mb-2">
                    Verification against host genome cross-reactivity
                  </p>
                  <div className="bg-muted p-3 rounded text-sm">
                    Target: {'<'}1% false positive rate
                  </div>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Accuracy Metrics</h3>
                  <p className="text-sm text-muted-foreground mb-2">
                    Clinical performance evaluation
                  </p>
                  <div className="bg-muted p-3 rounded text-sm space-y-1">
                    <div>• Positive Percent Agreement (PPA): {'>'}95%</div>
                    <div>• Negative Percent Agreement (NPA): {'>'}98%</div>
                    <div>• Quantification Correlation: R² {'>'}0.90</div>
                  </div>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Quality Management</h2>
              <p className="text-muted-foreground mb-4">
                The pipeline implements ALCOA+ principles for data integrity:
              </p>

              <div className="grid md:grid-cols-2 gap-4">
                {[
                  { letter: "A", term: "Attributable", desc: "All data traced to original source" },
                  { letter: "L", term: "Legible", desc: "Clear, permanent records" },
                  { letter: "C", term: "Contemporaneous", desc: "Recorded at time of activity" },
                  { letter: "O", term: "Original", desc: "First recording or certified copy" },
                  { letter: "A", term: "Accurate", desc: "Error-free, validated data" },
                  { letter: "+C", term: "Complete", desc: "All relevant information captured" },
                  { letter: "+C", term: "Consistent", desc: "Standardized format and timing" },
                  { letter: "+E", term: "Enduring", desc: "Retained for required period" },
                  { letter: "+A", term: "Available", desc: "Accessible for review/audit" }
                ].map((principle) => (
                  <div key={principle.term} className="flex items-start gap-3">
                    <div className="w-12 h-12 rounded-lg flex items-center justify-center font-bold text-sm"
                         style={{ backgroundColor: 'hsl(var(--primary) / 0.1)', color: 'hsl(var(--primary))' }}>
                      {principle.letter}
                    </div>
                    <div>
                      <div className="font-semibold text-sm">{principle.term}</div>
                      <div className="text-xs text-muted-foreground">{principle.desc}</div>
                    </div>
                  </div>
                ))}
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Compliance Reporting</h2>
              <p className="text-muted-foreground mb-4">
                The pipeline automatically generates PMDA-compliant reports in multiple formats:
              </p>

              <div className="space-y-3">
                <div className="border-l-2 border-primary pl-6">
                  <h3 className="font-semibold mb-1">PDF Reports</h3>
                  <p className="text-sm text-muted-foreground">
                    Comprehensive human-readable reports with executive summary
                  </p>
                </div>
                <div className="border-l-2 border-primary pl-6">
                  <h3 className="font-semibold mb-1">JSON Data</h3>
                  <p className="text-sm text-muted-foreground">
                    Structured machine-readable data for downstream analysis
                  </p>
                </div>
                <div className="border-l-2 border-primary pl-6">
                  <h3 className="font-semibold mb-1">HTML Dashboard</h3>
                  <p className="text-sm text-muted-foreground">
                    Interactive results visualization with drill-down capabilities
                  </p>
                </div>
                <div className="border-l-2 border-primary pl-6">
                  <h3 className="font-semibold mb-1">91-Pathogen Checklist</h3>
                  <p className="text-sm text-muted-foreground">
                    PMDA compliance checklist with detection status for each pathogen
                  </p>
                </div>
              </div>
            </section>

            <Alert style={{
              backgroundColor: 'hsl(var(--muted))',
              borderColor: 'hsl(var(--border))'
            }}>
              <ShieldCheckIcon className="h-4 w-4" />
              <div>
                <AlertTitle>Regulatory Compliance</AlertTitle>
                <AlertDescription>
                  For complete regulatory documentation and validation protocols, refer to the
                  strategic planning documents in the <code className="bg-background px-2 py-1 rounded">md/</code> directory
                  of the repository.
                </AlertDescription>
              </div>
            </Alert>
          </div>
        </main>
      </div>
    </div>
  );
}
