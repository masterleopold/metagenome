import { Sidebar } from "@/components/layout/Sidebar";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { CodeBlock } from "@/components/ui/CodeBlock";
import { CheckCircleIcon, AlertTriangleIcon, TerminalIcon, ServerIcon } from "lucide-react";

export default function GettingStartedPage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Introduction</Badge>
              <h1 className="text-4xl font-bold mb-4">Getting Started</h1>
              <p className="text-lg text-muted-foreground">
                Set up your development environment and deploy the MinION pathogen screening pipeline.
              </p>
            </div>

            <Alert className="mb-8" style={{
              backgroundColor: 'hsl(var(--primary) / 0.1)',
              borderColor: 'hsl(var(--primary) / 0.3)'
            }}>
              <CheckCircleIcon className="h-4 w-4" style={{ color: 'hsl(var(--primary))' }} />
              <div>
                <AlertTitle>Prerequisites Checklist</AlertTitle>
                <AlertDescription>
                  Make sure you have all required tools and access before starting.
                </AlertDescription>
              </div>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Prerequisites</h2>

              <div className="grid md:grid-cols-2 gap-6 mb-6">
                <Card>
                  <CardHeader>
                    <TerminalIcon className="h-8 w-8 mb-2 text-primary" />
                    <CardTitle>Local Tools</CardTitle>
                    <CardDescription>Required software on your machine</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>AWS CLI 2.0+</strong>
                          <p className="text-sm text-muted-foreground">Command line interface for AWS</p>
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>Terraform 1.0+</strong>
                          <p className="text-sm text-muted-foreground">Infrastructure as Code tool</p>
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>Python 3.9+</strong>
                          <p className="text-sm text-muted-foreground">Programming language runtime</p>
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>Git</strong>
                          <p className="text-sm text-muted-foreground">Version control system</p>
                        </div>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <ServerIcon className="h-8 w-8 mb-2 text-secondary" />
                    <CardTitle>AWS Resources</CardTitle>
                    <CardDescription>Required AWS access and quotas</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>AWS Account</strong>
                          <p className="text-sm text-muted-foreground">With admin/PowerUser access</p>
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>Service Quotas</strong>
                          <p className="text-sm text-muted-foreground">GPU instances, Lambda concurrency</p>
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 mt-0.5" style={{ color: 'hsl(var(--primary))' }} />
                        <div>
                          <strong>Region</strong>
                          <p className="text-sm text-muted-foreground">Recommended: ap-northeast-1</p>
                        </div>
                      </li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Installation</h2>

              <div className="space-y-6">
                <div>
                  <h3 className="text-lg font-medium mb-3">1. Clone the Repository</h3>
                  <CodeBlock
                    code={`git clone https://github.com/your-org/minion-pipeline.git
cd minion-pipeline`}
                    language="bash"
                  />
                </div>

                <div>
                  <h3 className="text-lg font-medium mb-3">2. Install Python Dependencies</h3>
                  <CodeBlock
                    code={`pip install -r requirements.txt`}
                    language="bash"
                  />
                </div>

                <div>
                  <h3 className="text-lg font-medium mb-3">3. Configure AWS Credentials</h3>
                  <CodeBlock
                    code={`aws configure
# Enter your AWS Access Key ID
# Enter your AWS Secret Access Key
# Default region name: ap-northeast-1
# Default output format: json`}
                    language="bash"
                  />
                </div>

                <div>
                  <h3 className="text-lg font-medium mb-3">4. Set Environment Variables</h3>
                  <CodeBlock
                    code={`export AWS_REGION=ap-northeast-1
export ENVIRONMENT=production
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)`}
                    language="bash"
                  />
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Deploy Infrastructure</h2>

              <Alert variant="warning" className="mb-6">
                <AlertTriangleIcon className="h-4 w-4" />
                <div>
                  <AlertTitle>Cost Warning</AlertTitle>
                  <AlertDescription>
                    Deploying this infrastructure will incur AWS charges. Review the cost estimation in the deployment guide.
                  </AlertDescription>
                </div>
              </Alert>

              <div className="space-y-6">
                <div>
                  <h3 className="text-lg font-medium mb-3">Initialize Terraform</h3>
                  <CodeBlock
                    code={`cd infrastructure/terraform
terraform init`}
                    language="bash"
                  />
                </div>

                <div>
                  <h3 className="text-lg font-medium mb-3">Plan Deployment</h3>
                  <CodeBlock
                    code={`terraform plan \\
  -var="environment=$ENVIRONMENT" \\
  -var="region=$AWS_REGION" \\
  -out=tfplan`}
                    language="bash"
                  />
                </div>

                <div>
                  <h3 className="text-lg font-medium mb-3">Apply Infrastructure</h3>
                  <CodeBlock
                    code={`terraform apply tfplan

# Save outputs for reference
terraform output -json > outputs.json`}
                    language="bash"
                  />
                  <p className="text-sm text-muted-foreground mt-2">
                    This will take approximately 10-15 minutes to complete.
                  </p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Setup Databases</h2>

              <div className="space-y-4">
                <p className="text-muted-foreground">
                  After infrastructure deployment, set up reference databases on EFS:
                </p>
                <CodeBlock
                  code={`# Run database setup script
./tools/database_setup.sh --all

# Verify installation
./tools/database_setup.sh --check`}
                  language="bash"
                />
                <p className="text-sm text-muted-foreground">
                  This will download and install Kraken2, RVDB, BLAST databases, and PMDA pathogen sequences.
                  Total download size: ~150GB. Allow 2-4 hours for completion.
                </p>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Verify Deployment</h2>

              <div className="space-y-4">
                <CodeBlock
                  code={`# Run deployment validation
./tools/deployment_script.sh validate

# Check all services status
./tools/deployment_script.sh status`}
                  language="bash"
                />

                <Card style={{
                  backgroundColor: 'hsl(var(--primary) / 0.1)',
                  borderColor: 'hsl(var(--primary) / 0.3)'
                }}>
                  <CardHeader>
                    <CheckCircleIcon className="h-8 w-8 mb-2" style={{ color: 'hsl(var(--primary))' }} />
                    <CardTitle>Deployment Successful</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="mb-4" style={{ color: 'hsl(var(--foreground))' }}>
                      If all checks pass, your MinION pipeline is ready to use!
                    </p>
                    <p className="text-sm text-muted-foreground">
                      Next steps: Run your first workflow or explore the API documentation.
                    </p>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Run Your First Workflow</h2>

              <div className="space-y-4">
                <p className="text-muted-foreground">
                  Upload test FAST5 files and start an analysis workflow:
                </p>
                <CodeBlock
                  code={`# Upload FAST5 files to S3
aws s3 cp test-data/sample.fast5 s3://minion-data-production/runs/TEST-001/fast5/

# Start workflow via CLI
./tools/workflow_cli.py start \\
  --run-id TEST-001 \\
  --bucket minion-data-production \\
  --input-prefix runs/TEST-001/fast5/

# Monitor progress
./tools/workflow_cli.py status --run-id TEST-001 --watch`}
                  language="bash"
                />
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-4">Next Steps</h2>
              <div className="grid md:grid-cols-2 gap-4">
                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>Architecture Overview</CardTitle>
                    <CardDescription>
                      Learn about system architecture and data flow
                    </CardDescription>
                  </CardHeader>
                </Card>
                <Card className="hover:shadow-lg transition-shadow cursor-pointer">
                  <CardHeader>
                    <CardTitle>API Reference</CardTitle>
                    <CardDescription>
                      Explore API endpoints and integration options
                    </CardDescription>
                  </CardHeader>
                </Card>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
