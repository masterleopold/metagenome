import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { AlertTriangleIcon, CheckCircleIcon, InfoIcon } from "lucide-react";

export default function DeploymentPage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Development</Badge>
              <h1 className="text-4xl font-bold mb-4">Deployment</h1>
              <p className="text-lg text-muted-foreground">
                Deploy and manage the MinION pipeline infrastructure on AWS.
              </p>
            </div>

            <Alert className="mb-8" style={{
              backgroundColor: 'hsl(var(--muted))',
              borderColor: 'hsl(var(--border))'
            }}>
              <InfoIcon className="h-4 w-4" />
              <div>
                <AlertTitle>Prerequisites</AlertTitle>
                <AlertDescription>
                  Make sure you've completed the Getting Started guide before deploying to production.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Deployment Environments</h2>
              <p className="text-muted-foreground mb-6">
                The pipeline supports multiple deployment environments with isolated resources:
              </p>

              <div className="space-y-4">
                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Development</h3>
                  <p className="text-sm text-muted-foreground mb-3">
                    For testing and development with minimal resources and cost
                  </p>
                  <code className="bg-muted px-3 py-1 rounded text-sm">
                    export ENVIRONMENT=development
                  </code>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Staging</h3>
                  <p className="text-sm text-muted-foreground mb-3">
                    Pre-production environment for validation testing
                  </p>
                  <code className="bg-muted px-3 py-1 rounded text-sm">
                    export ENVIRONMENT=staging
                  </code>
                </div>

                <div className="border rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-2">Production</h3>
                  <p className="text-sm text-muted-foreground mb-3">
                    Production environment with full resources and redundancy
                  </p>
                  <code className="bg-muted px-3 py-1 rounded text-sm">
                    export ENVIRONMENT=production
                  </code>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Infrastructure Components</h2>
              <p className="text-muted-foreground mb-6">
                The deployment creates the following AWS resources:
              </p>

              <div className="space-y-3">
                {[
                  { name: "S3 Buckets", desc: "Data storage with lifecycle policies" },
                  { name: "Lambda Functions", desc: "Orchestration and API endpoints" },
                  { name: "Step Functions", desc: "Workflow state machine" },
                  { name: "EC2 AMIs", desc: "Pre-configured compute images for each phase" },
                  { name: "RDS Aurora", desc: "Serverless PostgreSQL database" },
                  { name: "EFS", desc: "Reference database storage" },
                  { name: "CloudWatch", desc: "Monitoring dashboards and alarms" },
                  { name: "SNS Topics", desc: "Email/SMS notifications" },
                  { name: "IAM Roles", desc: "Service permissions and policies" },
                  { name: "VPC & Networking", desc: "Isolated network configuration" }
                ].map((component) => (
                  <div key={component.name} className="flex items-start gap-3 border-l-2 pl-4">
                    <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                    <div>
                      <strong className="text-sm">{component.name}</strong>
                      <p className="text-sm text-muted-foreground">{component.desc}</p>
                    </div>
                  </div>
                ))}
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Deployment Process</h2>

              <div className="space-y-8">
                <div>
                  <h3 className="text-lg font-semibold mb-3">1. Initialize Terraform Backend</h3>
                  <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                    cd infrastructure/terraform<br />
                    terraform init
                  </div>
                </div>

                <div>
                  <h3 className="text-lg font-semibold mb-3">2. Configure Variables</h3>
                  <div className="bg-muted p-4 rounded-lg font-mono text-sm space-y-1">
                    <div>export AWS_REGION=ap-northeast-1</div>
                    <div>export ENVIRONMENT=production</div>
                    <div>export PROJECT_NAME=minion-pipeline</div>
                  </div>
                </div>

                <div>
                  <h3 className="text-lg font-semibold mb-3">3. Plan Infrastructure</h3>
                  <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                    terraform plan -out=tfplan
                  </div>
                  <p className="text-sm text-muted-foreground mt-2">
                    Review the planned changes carefully before applying.
                  </p>
                </div>

                <div>
                  <h3 className="text-lg font-semibold mb-3">4. Deploy Infrastructure</h3>
                  <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                    terraform apply tfplan
                  </div>
                  <p className="text-sm text-muted-foreground mt-2">
                    Deployment takes approximately 10-15 minutes.
                  </p>
                </div>

                <div>
                  <h3 className="text-lg font-semibold mb-3">5. Setup Reference Databases</h3>
                  <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                    ./tools/database_setup.sh --all
                  </div>
                  <p className="text-sm text-muted-foreground mt-2">
                    Downloads and installs Kraken2, RVDB, and PMDA databases (~150GB, 2-4 hours).
                  </p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Post-Deployment Validation</h2>
              <div className="bg-muted p-4 rounded-lg font-mono text-sm space-y-1 mb-4">
                <div># Validate deployment</div>
                <div>./tools/deployment_script.sh validate</div>
                <div className="pt-2"># Check services status</div>
                <div>./tools/deployment_script.sh status</div>
              </div>
              <p className="text-sm text-muted-foreground">
                All services should report as healthy before processing production workloads.
              </p>
            </section>

            <Alert className="mb-8" style={{
              backgroundColor: 'hsl(var(--primary) / 0.1)',
              borderColor: 'hsl(var(--primary) / 0.3)'
            }}>
              <CheckCircleIcon className="h-4 w-4" style={{ color: 'hsl(var(--primary))' }} />
              <div>
                <AlertTitle>Deployment Complete</AlertTitle>
                <AlertDescription>
                  Your MinION pipeline is ready to process sequencing data. See the API Reference
                  for integration options.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Cost Estimation</h2>
              <p className="text-muted-foreground mb-4">
                Estimated monthly costs for production environment:
              </p>
              <div className="border rounded-lg p-6 space-y-3">
                <div className="flex justify-between">
                  <span>Lambda + Step Functions</span>
                  <span className="font-semibold">~$50/month</span>
                </div>
                <div className="flex justify-between">
                  <span>EC2 (on-demand, per run)</span>
                  <span className="font-semibold">~$15-25/run</span>
                </div>
                <div className="flex justify-between">
                  <span>RDS Aurora Serverless</span>
                  <span className="font-semibold">~$30-50/month</span>
                </div>
                <div className="flex justify-between">
                  <span>S3 Storage</span>
                  <span className="font-semibold">~$20-40/month</span>
                </div>
                <div className="flex justify-between">
                  <span>EFS (reference databases)</span>
                  <span className="font-semibold">~$50/month</span>
                </div>
                <div className="border-t pt-3 flex justify-between font-semibold">
                  <span>Total (24 analyses/year)</span>
                  <span className="text-primary">~$500-800/month</span>
                </div>
              </div>
            </section>

            <Alert style={{
              backgroundColor: 'hsl(var(--destructive) / 0.1)',
              borderColor: 'hsl(var(--destructive) / 0.3)'
            }}>
              <AlertTriangleIcon className="h-4 w-4" style={{ color: 'hsl(var(--destructive))' }} />
              <div>
                <AlertTitle>Important</AlertTitle>
                <AlertDescription>
                  Remember to destroy development/staging environments when not in use to avoid
                  unnecessary costs: <code className="bg-muted px-2 py-1 rounded">terraform destroy</code>
                </AlertDescription>
              </div>
            </Alert>
          </div>
        </main>
      </div>
    </div>
  );
}
