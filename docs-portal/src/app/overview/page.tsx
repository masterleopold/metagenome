import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import {
  CheckCircleIcon,
  ZapIcon,
  ShieldCheckIcon,
  CloudIcon,
  DatabaseIcon,
  GitBranchIcon,
} from "lucide-react";

export default function OverviewPage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Introduction</Badge>
              <h1 className="text-4xl font-bold mb-4">Overview</h1>
              <p className="text-lg text-muted-foreground">
                Comprehensive introduction to the MinION pathogen screening pipeline for xenotransplantation.
              </p>
            </div>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">What is the MinION Pipeline?</h2>

              <Card className="mb-6">
                <CardContent className="pt-6">
                  <p className="text-base leading-relaxed mb-4">
                    The MinION Pathogen Screening Pipeline is a PMDA-compliant metagenomic analysis system designed for comprehensive pathogen detection in donor pigs intended for xenotransplantation. The system uses Oxford Nanopore MinION Mk1D long-read sequencing technology combined with AWS cloud infrastructure for scalable, cost-effective analysis.
                  </p>
                  <p className="text-base leading-relaxed text-muted-foreground">
                    This pipeline provides end-to-end automation from basecalling to reporting, ensuring complete coverage of all 91 PMDA-designated pathogens with special emphasis on critical PERV (Porcine Endogenous Retrovirus) detection.
                  </p>
                </CardContent>
              </Card>

              <div className="grid md:grid-cols-3 gap-6">
                <Card>
                  <CardHeader>
                    <ShieldCheckIcon className="h-8 w-8 text-primary mb-2" />
                    <CardTitle>PMDA Compliant</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm text-muted-foreground">
                    Full coverage of 91 designated pathogens for xenotransplantation safety
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CloudIcon className="h-8 w-8 text-blue-600 mb-2" />
                    <CardTitle>Cloud-Native</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm text-muted-foreground">
                    Serverless architecture on AWS with Lambda + EC2 on-demand
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <ZapIcon className="h-8 w-8 text-amber-600 mb-2" />
                    <CardTitle>Real-time Analysis</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm text-muted-foreground">
                    Streaming analysis capability with MinION sequencing
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Key Features</h2>

              <div className="space-y-4">
                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <CheckCircleIcon className="h-6 w-6 text-green-600" />
                      <CardTitle>Q30+ Accuracy</CardTitle>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Duplex basecalling with Dorado achieves 99.9% accuracy (Q30+) by sequencing both DNA strands.
                    </p>
                    <div className="bg-muted p-3 rounded text-sm">
                      <strong>Simplex mode:</strong> 95% accuracy (Q12-15)<br />
                      <strong>Duplex mode:</strong> 99.9% accuracy (Q30+)
                    </div>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <DatabaseIcon className="h-6 w-6 text-purple-600" />
                      <CardTitle>Multi-Database Search</CardTitle>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Comprehensive pathogen detection using multiple complementary databases:
                    </p>
                    <ul className="text-sm space-y-1">
                      <li>• <strong>Kraken2:</strong> Rapid taxonomic classification (Standard DB ~50GB)</li>
                      <li>• <strong>RVDB v30.0:</strong> Manually curated viral sequences (~10GB)</li>
                      <li>• <strong>BLAST:</strong> PMDA custom database (91 pathogens, ~5GB)</li>
                      <li>• <strong>Diamond:</strong> Protein-level viral search</li>
                      <li>• <strong>PERV-specific:</strong> Targeted detection and typing</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <GitBranchIcon className="h-6 w-6 text-teal-600" />
                      <CardTitle>Automated Workflow</CardTitle>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Six automated phases from raw data to PMDA-compliant reports:
                    </p>
                    <div className="space-y-2">
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>1</Badge>
                        <span>Basecalling (2-4 hours)</span>
                      </div>
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>2</Badge>
                        <span>Quality Control (10-15 min)</span>
                      </div>
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>3</Badge>
                        <span>Host Removal (30-60 min)</span>
                      </div>
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>4</Badge>
                        <span>Pathogen Detection (1-2 hours)</span>
                      </div>
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>5</Badge>
                        <span>Quantification (15-30 min)</span>
                      </div>
                      <div className="flex items-center gap-2 text-sm">
                        <Badge>6</Badge>
                        <span>Reporting (10-15 min)</span>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Use Cases</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Xenotransplantation Research</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Primary use case for screening donor pigs before organ transplantation to humans.
                    </p>
                    <ul className="text-sm space-y-1">
                      <li>✓ Pre-transplant donor screening</li>
                      <li>✓ Quarterly monitoring of breeding colony</li>
                      <li>✓ Clinical trial donor qualification</li>
                      <li>✓ Post-transplant monitoring</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Pathogen Surveillance</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Broader applications in veterinary and agricultural settings.
                    </p>
                    <ul className="text-sm space-y-1">
                      <li>✓ Outbreak investigation</li>
                      <li>✓ Breeding facility monitoring</li>
                      <li>✓ Biosecurity assessment</li>
                      <li>✓ Research animal screening</li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">System Requirements</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Hardware</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="grid md:grid-cols-2 gap-6">
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">Sequencing</h4>
                        <ul className="text-sm space-y-1 text-muted-foreground">
                          <li>• Oxford Nanopore MinION Mk1D</li>
                          <li>• R10.4.1 flow cells</li>
                          <li>• Laptop/desktop for MinKNOW</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">AWS Resources</h4>
                        <ul className="text-sm space-y-1 text-muted-foreground">
                          <li>• g4dn.xlarge (GPU basecalling)</li>
                          <li>• r5.4xlarge (pathogen detection)</li>
                          <li>• t3.large/xlarge (general analysis)</li>
                        </ul>
                      </div>
                    </div>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Software Dependencies</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="grid md:grid-cols-3 gap-4 text-sm">
                      <div>
                        <h4 className="font-semibold mb-2">Basecalling</h4>
                        <ul className="space-y-1 text-muted-foreground">
                          <li>• Dorado v0.5+</li>
                          <li>• CUDA 11.8+</li>
                          <li>• PycoQC</li>
                          <li>• NanoPlot</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2">Alignment</h4>
                        <ul className="space-y-1 text-muted-foreground">
                          <li>• Minimap2</li>
                          <li>• SAMtools</li>
                          <li>• BWA-MEM</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2">Detection</h4>
                        <ul className="space-y-1 text-muted-foreground">
                          <li>• Kraken2/Bracken</li>
                          <li>• BLAST+</li>
                          <li>• Diamond</li>
                        </ul>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Performance Metrics</h2>

              <div className="grid md:grid-cols-4 gap-6">
                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription className="text-xs">Total Duration</CardDescription>
                    <CardTitle className="text-3xl">4-8 hrs</CardTitle>
                  </CardHeader>
                  <CardContent className="text-xs text-muted-foreground">
                    Complete analysis pipeline
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription className="text-xs">Accuracy</CardDescription>
                    <CardTitle className="text-3xl">99.9%</CardTitle>
                  </CardHeader>
                  <CardContent className="text-xs text-muted-foreground">
                    Duplex basecalling (Q30+)
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription className="text-xs">Sensitivity</CardDescription>
                    <CardTitle className="text-3xl">&gt;95%</CardTitle>
                  </CardHeader>
                  <CardContent className="text-xs text-muted-foreground">
                    For PMDA pathogens
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription className="text-xs">Cost per Run</CardDescription>
                    <CardTitle className="text-3xl">~$50</CardTitle>
                  </CardHeader>
                  <CardContent className="text-xs text-muted-foreground">
                    AWS compute costs
                  </CardContent>
                </Card>
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-6">Next Steps</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>Quick Start</CardTitle>
                    <CardDescription>
                      Get started with deployment and configuration
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground">
                      Follow our step-by-step guide to deploy the pipeline to AWS and run your first analysis.
                    </p>
                  </CardContent>
                </Card>

                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>Architecture Deep Dive</CardTitle>
                    <CardDescription>
                      Learn about system design and components
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground">
                      Explore the cloud-native architecture and understand how each component works together.
                    </p>
                  </CardContent>
                </Card>

                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>API Integration</CardTitle>
                    <CardDescription>
                      Integrate with your existing systems
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground">
                      Use our REST API to automate workflow execution and retrieve results programmatically.
                    </p>
                  </CardContent>
                </Card>

                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>PMDA Compliance</CardTitle>
                    <CardDescription>
                      Understand regulatory requirements
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground">
                      Review the 91 designated pathogens and PERV-specific detection requirements.
                    </p>
                  </CardContent>
                </Card>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
