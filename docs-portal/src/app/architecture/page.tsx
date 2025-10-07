import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { InfoIcon } from "lucide-react";

export default function ArchitecturePage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Core Concepts</Badge>
              <h1 className="text-4xl font-bold mb-4">System Architecture</h1>
              <p className="text-lg text-muted-foreground">
                AWS serverless architecture for scalable MinION pathogen screening.
              </p>
            </div>

            <Alert className="mb-8">
              <InfoIcon className="h-4 w-4" />
              <div>
                <AlertTitle>Documentation in progress</AlertTitle>
                <AlertDescription>
                  This page is being actively developed. Check back soon for complete documentation.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Overview</h2>
              <p className="text-muted-foreground">
                The MinION Pipeline uses a serverless AWS architecture designed for scalability,
                cost-efficiency, and PMDA compliance. The system processes Oxford Nanopore sequencing
                data through six sequential phases, with automated orchestration and monitoring.
              </p>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">High-Level Architecture</h2>
              <div className="bg-muted p-6 rounded-lg font-mono text-sm space-y-2">
                <div>MinION Sequencer → S3 Upload → Lambda Trigger → Step Functions</div>
                <div className="pl-8">↓</div>
                <div className="pl-4">[6 Sequential EC2 Phases]</div>
                <div className="pl-8">↓</div>
                <div className="pl-4">Reports + SNS Notifications</div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Core Components</h2>
              <div className="grid gap-4 md:grid-cols-2">
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon S3</h3>
                  <p className="text-sm text-muted-foreground">
                    Raw data storage, intermediate results, and final reports with lifecycle policies
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">AWS Lambda</h3>
                  <p className="text-sm text-muted-foreground">
                    Orchestration logic, API endpoints, and phase coordination
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">AWS Step Functions</h3>
                  <p className="text-sm text-muted-foreground">
                    Workflow state machine managing the 6-phase pipeline execution
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon EC2</h3>
                  <p className="text-sm text-muted-foreground">
                    GPU-enabled compute instances for basecalling and bioinformatics analysis
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon EFS</h3>
                  <p className="text-sm text-muted-foreground">
                    Shared reference database storage (Kraken2, RVDB, PMDA custom DB)
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon RDS Aurora</h3>
                  <p className="text-sm text-muted-foreground">
                    Serverless PostgreSQL for workflow metadata and run tracking
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon CloudWatch</h3>
                  <p className="text-sm text-muted-foreground">
                    Metrics, logs, dashboards, and alerting for pipeline monitoring
                  </p>
                </div>
                <div className="border rounded-lg p-6 space-y-2">
                  <h3 className="font-semibold text-lg">Amazon SNS</h3>
                  <p className="text-sm text-muted-foreground">
                    Email/SMS notifications for critical events (PERV detection, failures)
                  </p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Data Flow</h2>
              <ol className="list-decimal list-inside space-y-2 text-muted-foreground leading-relaxed">
                <li>Raw FAST5/POD5 files uploaded to S3 bucket</li>
                <li>S3 event triggers Lambda orchestrator</li>
                <li>Lambda starts Step Functions state machine</li>
                <li>Each phase launches EC2 instance with pre-configured AMI</li>
                <li>EC2 processes data, uploads results to S3, updates DynamoDB</li>
                <li>Instance auto-terminates on completion</li>
                <li>Final reports generated and SNS notifications sent</li>
              </ol>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Infrastructure as Code</h2>
              <p className="text-muted-foreground mb-4">
                All infrastructure is defined using Terraform, enabling version control,
                reproducible deployments, and environment isolation (development, staging, production).
              </p>
              <div className="bg-muted p-4 rounded-lg">
                <p className="text-sm font-mono">infrastructure/terraform/</p>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  )
}
