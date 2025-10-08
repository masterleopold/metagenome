import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import {
  CloudIcon,
  DatabaseIcon,
  ServerIcon,
  WorkflowIcon,
  LayersIcon,
  ZapIcon,
  ShieldIcon,
  InfoIcon,
} from "lucide-react";

export default function ArchitecturePage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Core Concepts</Badge>
              <h1 className="text-4xl font-bold mb-4">System Architecture</h1>
              <p className="text-lg text-muted-foreground">
                Understand the cloud-native architecture and component interactions of the MinION pipeline.
              </p>
            </div>

            <Alert className="mb-8">
              <InfoIcon className="h-4 w-4" />
              <AlertTitle>Architecture Overview</AlertTitle>
              <AlertDescription>
                The MinION pipeline uses a serverless architecture on AWS, combining Lambda functions for orchestration and EC2 instances for compute-intensive analysis.
              </AlertDescription>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">High-Level Architecture</h2>

              <div className="bg-slate-50 dark:bg-slate-900 rounded-xl p-8 mb-6 border">
                <div className="space-y-4">
                  <div className="flex items-center gap-4">
                    <div className="w-32 bg-primary/10 border-2 border-primary rounded-lg p-3 text-center">
                      <div className="text-sm font-semibold">MinION</div>
                      <div className="text-xs text-muted-foreground">Sequencer</div>
                    </div>
                    <div className="flex-1 border-t-2 border-dashed border-primary"></div>
                    <div className="w-32 bg-blue-100 dark:bg-blue-900/30 border-2 border-blue-500 rounded-lg p-3 text-center">
                      <div className="text-sm font-semibold">S3</div>
                      <div className="text-xs text-muted-foreground">Upload</div>
                    </div>
                  </div>

                  <div className="flex items-center gap-4">
                    <div className="w-full bg-purple-100 dark:bg-purple-900/30 border-2 border-purple-500 rounded-lg p-3 text-center">
                      <div className="text-sm font-semibold">Lambda Orchestration</div>
                      <div className="text-xs text-muted-foreground">Event-driven workflow management</div>
                    </div>
                  </div>

                  <div className="flex items-center gap-4">
                    <div className="w-full bg-green-100 dark:bg-green-900/30 border-2 border-green-500 rounded-lg p-3 text-center">
                      <div className="text-sm font-semibold">Step Functions Workflow</div>
                      <div className="text-xs text-muted-foreground">6-phase analysis pipeline</div>
                    </div>
                  </div>

                  <div className="grid grid-cols-3 gap-2">
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">Basecalling</div>
                      <div className="text-[10px] text-muted-foreground">GPU EC2</div>
                    </div>
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">QC</div>
                      <div className="text-[10px] text-muted-foreground">EC2</div>
                    </div>
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">Host Removal</div>
                      <div className="text-[10px] text-muted-foreground">EC2</div>
                    </div>
                  </div>

                  <div className="grid grid-cols-3 gap-2">
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">Pathogen Detection</div>
                      <div className="text-[10px] text-muted-foreground">High-mem EC2</div>
                    </div>
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">Quantification</div>
                      <div className="text-[10px] text-muted-foreground">EC2</div>
                    </div>
                    <div className="bg-orange-100 dark:bg-orange-900/30 border border-orange-500 rounded p-2 text-center text-xs">
                      <div className="font-semibold">Reporting</div>
                      <div className="text-[10px] text-muted-foreground">EC2</div>
                    </div>
                  </div>

                  <div className="flex items-center gap-4">
                    <div className="flex-1 bg-secondary/10 border-2 border-secondary rounded-lg p-3 text-center">
                      <div className="text-sm font-semibold">Reports</div>
                      <div className="text-xs text-muted-foreground">PDF, JSON, HTML</div>
                    </div>
                  </div>
                </div>
              </div>

              <p className="text-sm text-muted-foreground">
                <strong>Data Flow:</strong> MinION Sequencer → S3 Upload → Lambda Trigger → Step Functions → EC2 Processing (6 phases) → Results Storage → Report Generation
              </p>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">AWS Services</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card>
                  <CardHeader>
                    <CloudIcon className="h-8 w-8 text-blue-600 mb-2" />
                    <CardTitle>Compute</CardTitle>
                    <CardDescription>Scalable processing power</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2 text-sm">
                      <li className="flex justify-between">
                        <span className="font-medium">Lambda</span>
                        <Badge variant="secondary">Orchestration</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">EC2 (g4dn.xlarge)</span>
                        <Badge variant="secondary">GPU Basecalling</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">EC2 (r5.4xlarge)</span>
                        <Badge variant="secondary">High-memory</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">Step Functions</span>
                        <Badge variant="secondary">Workflow</Badge>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <DatabaseIcon className="h-8 w-8 text-green-600 mb-2" />
                    <CardTitle>Storage</CardTitle>
                    <CardDescription>Data persistence and databases</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2 text-sm">
                      <li className="flex justify-between">
                        <span className="font-medium">S3</span>
                        <Badge variant="secondary">Object Storage</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">EFS</span>
                        <Badge variant="secondary">Reference DBs</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">RDS Aurora</span>
                        <Badge variant="secondary">Metadata</Badge>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <WorkflowIcon className="h-8 w-8 text-purple-600 mb-2" />
                    <CardTitle>Integration</CardTitle>
                    <CardDescription>Event-driven communication</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2 text-sm">
                      <li className="flex justify-between">
                        <span className="font-medium">EventBridge</span>
                        <Badge variant="secondary">Events</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">SNS</span>
                        <Badge variant="secondary">Notifications</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">API Gateway</span>
                        <Badge variant="secondary">REST API</Badge>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <ServerIcon className="h-8 w-8 text-orange-600 mb-2" />
                    <CardTitle>Monitoring</CardTitle>
                    <CardDescription>Observability and logging</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2 text-sm">
                      <li className="flex justify-between">
                        <span className="font-medium">CloudWatch</span>
                        <Badge variant="secondary">Metrics & Logs</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">CloudWatch Alarms</span>
                        <Badge variant="secondary">Alerts</Badge>
                      </li>
                      <li className="flex justify-between">
                        <span className="font-medium">X-Ray</span>
                        <Badge variant="secondary">Tracing</Badge>
                      </li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Component Details</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <LayersIcon className="h-6 w-6 text-primary" />
                      <div>
                        <CardTitle>Lambda Functions (16 functions)</CardTitle>
                        <CardDescription>Serverless orchestration and control</CardDescription>
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <div className="grid md:grid-cols-2 gap-4">
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">Orchestration</h4>
                        <ul className="space-y-1 text-sm text-muted-foreground">
                          <li>• Pipeline orchestrator</li>
                          <li>• Phase state manager</li>
                          <li>• Error handler</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">EC2 Management</h4>
                        <ul className="space-y-1 text-sm text-muted-foreground">
                          <li>• Instance launcher</li>
                          <li>• Instance monitor</li>
                          <li>• Instance terminator</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">Data Processing</h4>
                        <ul className="space-y-1 text-sm text-muted-foreground">
                          <li>• FASTQ validator</li>
                          <li>• Result aggregator</li>
                          <li>• Metric calculator</li>
                        </ul>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2 text-sm">Monitoring</h4>
                        <ul className="space-y-1 text-sm text-muted-foreground">
                          <li>• Alert handler</li>
                          <li>• Status updater</li>
                          <li>• Cost tracker</li>
                        </ul>
                      </div>
                    </div>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <ZapIcon className="h-6 w-6 text-orange-600" />
                      <div>
                        <CardTitle>EC2 Instances</CardTitle>
                        <CardDescription>On-demand compute for analysis phases</CardDescription>
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <div className="space-y-4">
                      <div className="flex items-start gap-3 pb-3 border-b">
                        <Badge variant="warning">GPU</Badge>
                        <div className="flex-1">
                          <div className="font-semibold">g4dn.xlarge</div>
                          <div className="text-sm text-muted-foreground">Basecalling with Dorado (NVIDIA T4 GPU)</div>
                          <div className="text-xs text-muted-foreground mt-1">4 vCPU, 16GB RAM, 1x T4 GPU</div>
                        </div>
                      </div>
                      <div className="flex items-start gap-3 pb-3 border-b">
                        <Badge variant="secondary">Memory</Badge>
                        <div className="flex-1">
                          <div className="font-semibold">r5.4xlarge</div>
                          <div className="text-sm text-muted-foreground">Pathogen detection (Kraken2, BLAST)</div>
                          <div className="text-xs text-muted-foreground mt-1">16 vCPU, 128GB RAM</div>
                        </div>
                      </div>
                      <div className="flex items-start gap-3">
                        <Badge>General</Badge>
                        <div className="flex-1">
                          <div className="font-semibold">t3.large / r5.xlarge</div>
                          <div className="text-sm text-muted-foreground">QC, Host Removal, Quantification, Reporting</div>
                          <div className="text-xs text-muted-foreground mt-1">2-4 vCPU, 8-32GB RAM</div>
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <div className="flex items-center gap-3">
                      <ShieldIcon className="h-6 w-6 text-green-600" />
                      <div>
                        <CardTitle>Security & IAM</CardTitle>
                        <CardDescription>Access control and data protection</CardDescription>
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent>
                    <ul className="space-y-2 text-sm">
                      <li className="flex items-start gap-2">
                        <span className="text-green-600 mt-1">✓</span>
                        <div>
                          <strong>VPC Isolation:</strong> Private subnets for EC2 instances
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <span className="text-green-600 mt-1">✓</span>
                        <div>
                          <strong>IAM Roles:</strong> Least-privilege access for Lambda and EC2
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <span className="text-green-600 mt-1">✓</span>
                        <div>
                          <strong>S3 Encryption:</strong> Server-side encryption (SSE-S3)
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <span className="text-green-600 mt-1">✓</span>
                        <div>
                          <strong>RDS Encryption:</strong> At-rest encryption with KMS
                        </div>
                      </li>
                      <li className="flex items-start gap-2">
                        <span className="text-green-600 mt-1">✓</span>
                        <div>
                          <strong>API Gateway:</strong> API key authentication
                        </div>
                      </li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Data Flow</h2>

              <div className="space-y-4">
                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">1</div>
                  <div>
                    <h3 className="font-semibold mb-1">Data Upload</h3>
                    <p className="text-sm text-muted-foreground">
                      FAST5/POD5 files uploaded from MinION to S3 bucket (<code className="bg-muted px-1 py-0.5 rounded">runs/&#123;run_id&#125;/fast5/</code>)
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">2</div>
                  <div>
                    <h3 className="font-semibold mb-1">Event Trigger</h3>
                    <p className="text-sm text-muted-foreground">
                      S3 ObjectCreated event triggers Lambda orchestrator function
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">3</div>
                  <div>
                    <h3 className="font-semibold mb-1">Workflow Initialization</h3>
                    <p className="text-sm text-muted-foreground">
                      Lambda creates workflow record in RDS and starts Step Functions execution
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">4</div>
                  <div>
                    <h3 className="font-semibold mb-1">Phase Execution</h3>
                    <p className="text-sm text-muted-foreground">
                      Each phase triggers Lambda to launch EC2 instance with appropriate configuration
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">5</div>
                  <div>
                    <h3 className="font-semibold mb-1">Analysis Processing</h3>
                    <p className="text-sm text-muted-foreground">
                      EC2 instance downloads data from S3, processes with analysis tools, uploads results back to S3
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">6</div>
                  <div>
                    <h3 className="font-semibold mb-1">Phase Completion</h3>
                    <p className="text-sm text-muted-foreground">
                      EC2 instance signals completion to Lambda, which updates RDS and terminates instance
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">7</div>
                  <div>
                    <h3 className="font-semibold mb-1">Critical Alerts</h3>
                    <p className="text-sm text-muted-foreground">
                      PERV detection triggers immediate SNS notification to alert recipients
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-4">
                  <div className="bg-primary text-white rounded-full w-8 h-8 flex items-center justify-center font-bold shrink-0">8</div>
                  <div>
                    <h3 className="font-semibold mb-1">Report Generation</h3>
                    <p className="text-sm text-muted-foreground">
                      Final phase generates PMDA-compliant reports (PDF, JSON, HTML) and stores in S3
                    </p>
                  </div>
                </div>
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-6">Cost Optimization</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Spot Instances</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      70% cost reduction for basecalling and analysis phases
                    </p>
                    <ul className="space-y-1 text-sm">
                      <li>• Automatic fallback to on-demand</li>
                      <li>• Checkpoint/resume capability</li>
                      <li>• Priority-based instance selection</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Auto-termination</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Instances automatically terminate after phase completion
                    </p>
                    <ul className="space-y-1 text-sm">
                      <li>• No idle instance costs</li>
                      <li>• Lambda-based monitoring</li>
                      <li>• Configurable timeout protection</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Serverless Services</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Pay-per-use pricing for orchestration
                    </p>
                    <ul className="space-y-1 text-sm">
                      <li>• Lambda (millisecond billing)</li>
                      <li>• RDS Aurora Serverless</li>
                      <li>• S3 Intelligent-Tiering</li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Data Lifecycle</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Automated data retention and archival
                    </p>
                    <ul className="space-y-1 text-sm">
                      <li>• Raw data: 30-day retention</li>
                      <li>• Results: 90-day retention</li>
                      <li>• Reports: 365-day retention</li>
                    </ul>
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
