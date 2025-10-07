import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";

export default function OverviewPage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Introduction</Badge>
              <h1 className="text-4xl font-bold mb-4">Overview</h1>
              <p className="text-lg text-muted-foreground">
                High-level overview of the MinION pathogen screening pipeline for xenotransplantation.
              </p>
            </div>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Project Overview</h2>
              <p className="text-muted-foreground mb-4">
                This repository implements a production-ready MinION pathogen screening pipeline for
                xenotransplantation donor pig screening, with comprehensive strategic planning documentation.
                The system detects 91 PMDA-designated pathogens using Oxford Nanopore MinION long-read
                sequencing and AWS cloud infrastructure.
              </p>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Research Background</h2>
              <div className="space-y-4 text-muted-foreground">
                <div>
                  <h3 className="text-lg font-medium text-foreground mb-2">Target Animal</h3>
                  <p>
                    Yucatan miniature pig (genetically modified with 69 genomic edits: 3KO, 7TG,
                    59 PERV inactivations)
                  </p>
                </div>
                <div>
                  <h3 className="text-lg font-medium text-foreground mb-2">Sample Type</h3>
                  <p>Plasma samples for cfDNA/cfRNA extraction</p>
                </div>
                <div>
                  <h3 className="text-lg font-medium text-foreground mb-2">Regulatory Framework</h3>
                  <p>PMDA guidelines for xenotransplantation safety</p>
                </div>
                <div>
                  <h3 className="text-lg font-medium text-foreground mb-2">Primary Goal</h3>
                  <p>
                    Comprehensive detection of both known (91 PMDA-designated pathogens) and
                    unknown pathogens
                  </p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Three-Phase Workflow</h2>
              <div className="space-y-6">
                <div className="border-l-2 border-primary pl-6">
                  <h3 className="text-lg font-semibold mb-2">1. Lab Work Phase</h3>
                  <ul className="list-disc pl-6 space-y-2 text-muted-foreground">
                    <li>Blood plasma separation and preservation</li>
                    <li>cfDNA/cfRNA extraction using magnetic bead methods</li>
                    <li>Host DNA depletion using endonuclease-based methods</li>
                  </ul>
                </div>

                <div className="border-l-2 border-secondary pl-6">
                  <h3 className="text-lg font-semibold mb-2">2. Sequencing Phase</h3>
                  <ul className="list-disc pl-6 space-y-2 text-muted-foreground">
                    <li>Oxford Nanopore MinION long-read sequencing</li>
                    <li>Duplex mode for Q30+ accuracy</li>
                    <li>Real-time analysis capability</li>
                    <li>Lower initial investment compared to Illumina platforms</li>
                  </ul>
                </div>

                <div className="border-l-2 pl-6" style={{ borderColor: 'hsl(var(--primary))' }}>
                  <h3 className="text-lg font-semibold mb-2">3. Bioinformatics Phase</h3>
                  <ul className="list-disc pl-6 space-y-2 text-muted-foreground">
                    <li>Quality control (FastQC, PycoQC)</li>
                    <li>Host genome removal (Minimap2)</li>
                    <li>Pathogen detection via multi-database search</li>
                    <li>Quantitative analysis and reporting</li>
                    <li>AWS cloud-based pipeline for scalable processing</li>
                  </ul>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Key Features</h2>
              <div className="grid md:grid-cols-2 gap-6">
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">PMDA Compliance</h3>
                  <p className="text-sm text-muted-foreground">
                    Full coverage of 91 designated pathogens with regulatory-compliant reporting
                  </p>
                </div>
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">PERV Detection</h3>
                  <p className="text-sm text-muted-foreground">
                    Critical detection and typing of Porcine Endogenous Retroviruses (PERV-A, B, C)
                  </p>
                </div>
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Cloud-Native</h3>
                  <p className="text-sm text-muted-foreground">
                    Serverless AWS architecture with Lambda orchestration and EC2 compute
                  </p>
                </div>
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Quality Assured</h3>
                  <p className="text-sm text-muted-foreground">
                    Q30+ accuracy with duplex basecalling and comprehensive validation
                  </p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">System Architecture</h2>
              <p className="text-muted-foreground mb-4">
                The pipeline uses a serverless AWS architecture designed for scalability,
                cost-efficiency, and PMDA compliance:
              </p>
              <div className="bg-muted p-6 rounded-lg font-mono text-sm space-y-2">
                <div>MinION Sequencer → S3 Upload → Lambda Trigger → Step Functions</div>
                <div className="pl-8">↓</div>
                <div className="pl-4">[6 Sequential EC2 Phases]</div>
                <div className="pl-8">↓</div>
                <div className="pl-4">Reports + SNS Notifications</div>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
