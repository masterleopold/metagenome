import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { CodeBlock, InlineCode } from "@/components/ui/CodeBlock";
import { KeyIcon, InfoIcon } from "lucide-react";

export default function APIReferencePage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-8 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Development</Badge>
              <h1 className="text-4xl font-bold mb-4">API Reference</h1>
              <p className="text-lg text-muted-foreground">
                Complete REST API documentation for the MinION pathogen screening pipeline.
              </p>
            </div>

            <Alert className="mb-8">
              <InfoIcon className="h-4 w-4" />
              <AlertTitle>Base URL</AlertTitle>
              <AlertDescription>
                <code className="text-sm">https://api.minion-pipeline.com/&#123;environment&#125;/</code>
                <div className="mt-2 space-y-1 text-sm">
                  <div><strong>production:</strong> Production environment</div>
                  <div><strong>staging:</strong> Staging environment</div>
                  <div><strong>development:</strong> Development environment</div>
                </div>
              </AlertDescription>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Authentication</h2>

              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-3">
                    <KeyIcon className="h-6 w-6 text-primary" />
                    <div>
                      <CardTitle>API Key Authentication</CardTitle>
                      <CardDescription>All requests require an API key in the header</CardDescription>
                    </div>
                  </div>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    code={`# Include in request headers
x-api-key: YOUR_API_KEY`}
                    language="http"
                  />
                  <div className="mt-4 space-y-2 text-sm">
                    <p><strong>Obtaining an API Key:</strong></p>
                    <ol className="list-decimal list-inside space-y-1 text-muted-foreground ml-4">
                      <li>Contact your administrator for API key generation</li>
                      <li>Store the key securely in environment variables</li>
                      <li>Never commit API keys to version control</li>
                    </ol>
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Endpoints</h2>

              {/* POST /workflows */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge>POST</Badge>
                    <InlineCode>/workflows</InlineCode>
                  </div>
                  <CardTitle>Start Workflow</CardTitle>
                  <CardDescription>Start a new analysis workflow</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Request Body</h4>
                    <CodeBlock
                      code={`{
  "run_id": "RUN-2024-001",
  "bucket": "minion-data",
  "input_prefix": "runs/RUN-2024-001/fast5/",
  "config": {
    "phases": {
      "basecalling": {
        "enabled": true,
        "skip_duplex": false
      },
      "qc": {
        "enabled": true,
        "min_quality": 9
      },
      "host_removal": {
        "enabled": true,
        "reference": "sus_scrofa_11.1"
      },
      "pathogen_detection": {
        "enabled": true,
        "databases": ["kraken2", "rvdb", "pmda"]
      },
      "quantification": {
        "enabled": true,
        "spike_in": "PhiX174"
      },
      "reporting": {
        "enabled": true,
        "formats": ["pdf", "json"]
      }
    }
  }
}`}
                      language="json"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "run_id": "RUN-2024-001",
  "workflow_id": "123456",
  "execution_arn": "arn:aws:states:...",
  "status": "STARTED",
  "estimated_duration_hours": 8.5
}`}
                      language="json"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Example Request</h4>
                    <CodeBlock
                      code={`curl -X POST https://api.minion-pipeline.com/production/workflows \\
  -H "x-api-key: YOUR_API_KEY" \\
  -H "Content-Type: application/json" \\
  -d '{
    "run_id": "RUN-2024-001",
    "bucket": "minion-data",
    "input_prefix": "runs/RUN-2024-001/fast5/"
  }'`}
                      language="bash"
                    />
                  </div>
                </CardContent>
              </Card>

              {/* GET /workflows/{workflow_id} */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="secondary">GET</Badge>
                    <InlineCode>/workflows/&#123;workflow_id&#125;</InlineCode>
                  </div>
                  <CardTitle>Get Workflow Status</CardTitle>
                  <CardDescription>Get the status of a workflow execution</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Path Parameters</h4>
                    <ul className="text-sm space-y-1">
                      <li><InlineCode>workflow_id</InlineCode> - Unique workflow identifier</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "workflow_id": "123456",
  "run_id": "RUN-2024-001",
  "status": "RUNNING",
  "created_at": "2024-01-15T10:00:00Z",
  "updated_at": "2024-01-15T14:30:00Z",
  "phases": [
    {
      "name": "basecalling",
      "status": "COMPLETED",
      "started_at": "2024-01-15T10:05:00Z",
      "completed_at": "2024-01-15T12:00:00Z"
    },
    {
      "name": "qc",
      "status": "COMPLETED",
      "started_at": "2024-01-15T12:05:00Z",
      "completed_at": "2024-01-15T12:30:00Z"
    },
    {
      "name": "host_removal",
      "status": "RUNNING",
      "started_at": "2024-01-15T12:35:00Z"
    }
  ],
  "metrics": {
    "basecalling": {
      "total_reads": 50000,
      "mean_quality": 10.5,
      "total_bases": 150000000
    }
  },
  "pathogens": []
}`}
                      language="json"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Status Values</h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div className="flex items-center gap-2">
                        <Badge variant="outline">INITIATED</Badge>
                        <span className="text-muted-foreground">Created but not started</span>
                      </div>
                      <div className="flex items-center gap-2">
                        <Badge variant="warning">RUNNING</Badge>
                        <span className="text-muted-foreground">Currently executing</span>
                      </div>
                      <div className="flex items-center gap-2">
                        <Badge variant="success">COMPLETED</Badge>
                        <span className="text-muted-foreground">Successfully completed</span>
                      </div>
                      <div className="flex items-center gap-2">
                        <Badge variant="error">FAILED</Badge>
                        <span className="text-muted-foreground">Execution failed</span>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>

              {/* GET /workflows */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="secondary">GET</Badge>
                    <InlineCode>/workflows</InlineCode>
                  </div>
                  <CardTitle>List Workflows</CardTitle>
                  <CardDescription>List all workflow executions with optional filtering</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Query Parameters</h4>
                    <ul className="text-sm space-y-1">
                      <li><InlineCode>status</InlineCode> - Filter by status (RUNNING, COMPLETED, FAILED)</li>
                      <li><InlineCode>run_id</InlineCode> - Filter by run ID (partial match)</li>
                      <li><InlineCode>limit</InlineCode> - Maximum results (default: 20, max: 100)</li>
                      <li><InlineCode>offset</InlineCode> - Pagination offset</li>
                      <li><InlineCode>days</InlineCode> - Get workflows from last N days (default: 7)</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Example Request</h4>
                    <CodeBlock
                      code={`curl -X GET "https://api.minion-pipeline.com/production/workflows?status=COMPLETED&limit=10" \\
  -H "x-api-key: YOUR_API_KEY"`}
                      language="bash"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "workflows": [
    {
      "workflow_id": "123456",
      "run_id": "RUN-2024-001",
      "status": "COMPLETED",
      "created_at": "2024-01-15T10:00:00Z",
      "updated_at": "2024-01-15T18:00:00Z",
      "phases_completed": 6,
      "error_count": 0,
      "has_pathogens": true,
      "duration_seconds": 28800
    }
  ],
  "pagination": {
    "total": 45,
    "limit": 20,
    "offset": 0,
    "has_more": true
  }
}`}
                      language="json"
                    />
                  </div>
                </CardContent>
              </Card>

              {/* GET /pathogens/{run_id} */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="secondary">GET</Badge>
                    <InlineCode>/pathogens/&#123;run_id&#125;</InlineCode>
                  </div>
                  <CardTitle>Get Pathogen Results</CardTitle>
                  <CardDescription>Get pathogen detection results for a specific run</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Path Parameters</h4>
                    <ul className="text-sm space-y-1">
                      <li><InlineCode>run_id</InlineCode> - Run identifier</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Query Parameters</h4>
                    <ul className="text-sm space-y-1">
                      <li><InlineCode>pmda_only</InlineCode> - Show only PMDA 91 pathogens (default: false)</li>
                      <li><InlineCode>min_reads</InlineCode> - Minimum read count threshold (default: 0)</li>
                      <li><InlineCode>include_quantification</InlineCode> - Include absolute quantification (default: false)</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "run_id": "RUN-2024-001",
  "workflow_id": "123456",
  "total_pathogens": 8,
  "pmda_pathogens_detected": 5,
  "pmda_compliance": {
    "total": 91,
    "detected": 5,
    "not_detected": 86,
    "detection_rate": 0.055,
    "detected_list": ["HEV", "PCV2", "SS", "EC", "CA"],
    "critical_detected": false
  },
  "perv_analysis": {
    "detected": false
  },
  "pathogens": [
    {
      "code": "PCV2",
      "name": "Porcine circovirus 2",
      "read_count": 1250,
      "confidence": 0.98,
      "method": "kraken2",
      "is_pmda": true,
      "risk_level": "MEDIUM",
      "quantification": {
        "copies_per_ml": 5.2e4,
        "log10_copies": 4.72,
        "ci_lower": 4.8e4,
        "ci_upper": 5.6e4
      }
    }
  ]
}`}
                      language="json"
                    />
                  </div>
                </CardContent>
              </Card>

              {/* DELETE /workflows/{workflow_id} */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="error">DELETE</Badge>
                    <InlineCode>/workflows/&#123;workflow_id&#125;</InlineCode>
                  </div>
                  <CardTitle>Stop Workflow</CardTitle>
                  <CardDescription>Stop a running workflow execution</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Request Body (optional)</h4>
                    <CodeBlock
                      code={`{
  "reason": "User requested termination"
}`}
                      language="json"
                    />
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "workflow_id": "123456",
  "status": "STOPPING",
  "message": "Workflow termination initiated"
}`}
                      language="json"
                    />
                  </div>
                </CardContent>
              </Card>

              {/* GET /reports/{run_id} */}
              <Card className="mb-6">
                <CardHeader>
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="secondary">GET</Badge>
                    <InlineCode>/reports/&#123;run_id&#125;</InlineCode>
                  </div>
                  <CardTitle>Get Report URLs</CardTitle>
                  <CardDescription>Get analysis report download URLs</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <h4 className="font-semibold mb-2">Query Parameters</h4>
                    <ul className="text-sm space-y-1">
                      <li><InlineCode>format</InlineCode> - Report format (pdf, json, html)</li>
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-semibold mb-2">Response (200 OK)</h4>
                    <CodeBlock
                      code={`{
  "run_id": "RUN-2024-001",
  "reports": {
    "pdf": "https://s3.amazonaws.com/...",
    "json": "https://s3.amazonaws.com/...",
    "html": "https://s3.amazonaws.com/...",
    "checklist": "https://s3.amazonaws.com/..."
  },
  "expires_at": "2024-01-22T10:00:00Z"
}`}
                      language="json"
                    />
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Error Responses</h2>

              <Card className="mb-6">
                <CardHeader>
                  <CardTitle>Error Format</CardTitle>
                  <CardDescription>All error responses follow this format</CardDescription>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    code={`{
  "error": {
    "code": "VALIDATION_ERROR",
    "message": "Invalid input parameters",
    "details": {
      "field": "run_id",
      "issue": "Required field missing"
    }
  },
  "timestamp": "2024-01-15T10:00:00Z",
  "request_id": "abc123"
}`}
                    language="json"
                  />
                </CardContent>
              </Card>

              <div className="grid md:grid-cols-2 gap-4">
                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Error Codes</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-center gap-2">
                        <Badge variant="error">400</Badge>
                        <span>VALIDATION_ERROR</span>
                      </li>
                      <li className="flex items-center gap-2">
                        <Badge variant="error">401</Badge>
                        <span>AUTHENTICATION_ERROR</span>
                      </li>
                      <li className="flex items-center gap-2">
                        <Badge variant="error">403</Badge>
                        <span>AUTHORIZATION_ERROR</span>
                      </li>
                      <li className="flex items-center gap-2">
                        <Badge variant="error">404</Badge>
                        <span>NOT_FOUND</span>
                      </li>
                      <li className="flex items-center gap-2">
                        <Badge variant="error">429</Badge>
                        <span>RATE_LIMIT_EXCEEDED</span>
                      </li>
                      <li className="flex items-center gap-2">
                        <Badge variant="error">500</Badge>
                        <span>INTERNAL_ERROR</span>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Rate Limiting</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm space-y-2">
                    <p><strong>Default limit:</strong> 100 requests/minute</p>
                    <p><strong>Burst limit:</strong> 200 requests</p>
                    <div className="mt-3">
                      <p className="font-medium mb-1">Response Headers:</p>
                      <ul className="space-y-1 text-muted-foreground">
                        <li>• <InlineCode>X-RateLimit-Limit</InlineCode></li>
                        <li>• <InlineCode>X-RateLimit-Remaining</InlineCode></li>
                        <li>• <InlineCode>X-RateLimit-Reset</InlineCode></li>
                      </ul>
                    </div>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-6">SDK Examples</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Python</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`import requests

API_KEY = 'your-api-key'
BASE_URL = 'https://api.minion-pipeline.com/production'

headers = {'x-api-key': API_KEY}

# Start workflow
response = requests.post(
    f'{BASE_URL}/workflows',
    headers=headers,
    json={
        'run_id': 'RUN-2024-001',
        'bucket': 'minion-data',
        'input_prefix': 'runs/RUN-2024-001/fast5/'
    }
)

workflow = response.json()
print(f"Workflow started: {workflow['workflow_id']}")

# Check status
status = requests.get(
    f"{BASE_URL}/workflows/{workflow['workflow_id']}",
    headers=headers
).json()

print(f"Status: {status['status']}")`}
                      language="python"
                      filename="example.py"
                    />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Node.js</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`const axios = require('axios');

const API_KEY = 'your-api-key';
const BASE_URL = 'https://api.minion-pipeline.com/production';

const client = axios.create({
  baseURL: BASE_URL,
  headers: {'x-api-key': API_KEY}
});

// Start workflow
async function startWorkflow(runId, bucket, prefix) {
  const response = await client.post('/workflows', {
    run_id: runId,
    bucket: bucket,
    input_prefix: prefix
  });

  return response.data;
}

// Usage
startWorkflow('RUN-2024-001', 'minion-data', 'runs/RUN-2024-001/fast5/')
  .then(workflow => console.log(\`Started: \${workflow.workflow_id}\`))
  .catch(error => console.error(error));`}
                      language="javascript"
                      filename="example.js"
                    />
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
