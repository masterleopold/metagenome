import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { CodeBlock } from "@/components/ui/CodeBlock";
import {
  ActivityIcon,
  AlertTriangleIcon,
  CheckCircleIcon,
  DatabaseIcon,
  GlobeIcon,
  BellIcon,
  BarChartIcon,
  SearchIcon,
  ClockIcon,
} from "lucide-react";

const targetViruses = [
  {
    name: "Hantavirus",
    nameJa: "ハンタウイルス",
    priority: "HIGH",
    source: "Rodent-borne zoonotic",
    thresholdCritical: null,
    thresholdHigh: 100,
    status: "Monitored",
  },
  {
    name: "Polyomavirus",
    nameJa: "ポリオーマウイルス",
    priority: "MEDIUM",
    source: "Sus scrofa polyomavirus",
    thresholdCritical: null,
    thresholdHigh: 100,
    status: "Monitored",
  },
  {
    name: "Spumavirus",
    nameJa: "スピューマウイルス",
    priority: "CRITICAL",
    source: "MHLW Special Mgmt #5",
    thresholdCritical: 500,
    thresholdHigh: 100,
    status: "Monitored",
  },
  {
    name: "EEEV",
    nameJa: "東部ウマ脳炎ウイルス",
    priority: "CRITICAL",
    source: "NOT endemic to Japan",
    thresholdCritical: "ANY",
    thresholdHigh: null,
    status: "Monitored",
  },
];

const externalSources = [
  {
    name: "MAFF",
    fullName: "Ministry of Agriculture, Forestry and Fisheries",
    url: "www.maff.go.jp",
    frequency: "Daily",
    method: "Web Scraping",
    dataType: "Surveillance Reports",
  },
  {
    name: "E-Stat",
    fullName: "Government Statistics Portal",
    url: "www.e-stat.go.jp",
    frequency: "Daily",
    method: "REST API",
    dataType: "Livestock Statistics",
  },
  {
    name: "PubMed",
    fullName: "NCBI PubMed Database",
    url: "pubmed.ncbi.nlm.nih.gov",
    frequency: "Daily",
    method: "E-utilities API",
    dataType: "Research Publications",
  },
  {
    name: "J-STAGE",
    fullName: "Japan Science and Technology",
    url: "www.jstage.jst.go.jp",
    frequency: "Daily",
    method: "Web Scraping",
    dataType: "Japanese Publications",
  },
];

const severityLevels = [
  {
    level: "CRITICAL",
    color: "bg-red-50 dark:bg-red-900/10 text-red-600 dark:text-red-400",
    response: "< 5 min",
    actions: ["SNS Immediate Alert", "SMS to Key Personnel", "Dashboard Flashing", "Pipeline Pause"],
    criteria: "Spumavirus >500 copies/mL, ANY EEEV detection",
  },
  {
    level: "HIGH",
    color: "bg-orange-50 dark:bg-orange-900/10 text-orange-600 dark:text-orange-400",
    response: "< 30 min",
    actions: ["SNS Notification", "Email Alert", "Dashboard Warning"],
    criteria: "Hantavirus >100 copies/mL, Polyomavirus >100 copies/mL",
  },
  {
    level: "MEDIUM",
    color: "bg-blue-50 dark:bg-blue-900/10 text-blue-600 dark:text-blue-400",
    response: "< 2 hours",
    actions: ["Email Notification", "Dashboard Display"],
    criteria: "External keyword match, Low-level detections",
  },
  {
    level: "LOW",
    color: "bg-gray-50 dark:bg-gray-900/10 text-gray-600 dark:text-gray-400",
    response: "< 24 hours",
    actions: ["Dashboard Record", "Daily Summary"],
    criteria: "Academic publications, Informational",
  },
];

export default function SurveillancePage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            {/* Header */}
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Surveillance System</Badge>
              <h1 className="text-4xl font-bold mb-4">4-Virus Surveillance</h1>
              <p className="text-lg text-muted-foreground">
                Real-time monitoring system for Hantavirus, Polyomavirus, Spumavirus, and EEEV in Japanese pig populations.
              </p>
            </div>

            {/* System Overview Alert */}
            <Alert className="mb-8">
              <ActivityIcon className="h-4 w-4" />
              <AlertTitle>Integrated Monitoring</AlertTitle>
              <AlertDescription>
                This system combines external information sources (government, academic) with internal MinION pipeline detections to provide comprehensive surveillance of 4 target viruses.
              </AlertDescription>
            </Alert>

            {/* Target Viruses Section */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Target Viruses</h2>

              <div className="grid md:grid-cols-2 gap-4">
                {targetViruses.map((virus) => (
                  <Card key={virus.name}>
                    <CardHeader>
                      <div className="flex items-center justify-between">
                        <CardTitle className="text-lg">{virus.name}</CardTitle>
                        <Badge
                          variant={
                            virus.priority === "CRITICAL" ? "destructive" :
                            virus.priority === "HIGH" ? "default" :
                            "secondary"
                          }
                        >
                          {virus.priority}
                        </Badge>
                      </div>
                      <CardDescription>{virus.nameJa}</CardDescription>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-2 text-sm">
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Source:</span>
                          <span className="font-medium">{virus.source}</span>
                        </div>
                        {virus.thresholdCritical && (
                          <div className="flex justify-between">
                            <span className="text-muted-foreground">Critical:</span>
                            <span className="font-medium text-red-600 dark:text-red-400">
                              {typeof virus.thresholdCritical === 'number'
                                ? `>${virus.thresholdCritical} copies/mL`
                                : virus.thresholdCritical}
                            </span>
                          </div>
                        )}
                        {virus.thresholdHigh && (
                          <div className="flex justify-between">
                            <span className="text-muted-foreground">High:</span>
                            <span className="font-medium text-orange-600 dark:text-orange-400">
                              >{virus.thresholdHigh} copies/mL
                            </span>
                          </div>
                        )}
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Status:</span>
                          <span className="flex items-center gap-1">
                            <CheckCircleIcon className="h-3 w-3 text-green-600" />
                            {virus.status}
                          </span>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                ))}
              </div>
            </section>

            {/* Data Sources Section */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">External Information Sources</h2>

              <Card>
                <CardHeader>
                  <CardTitle>Daily Monitoring Sources</CardTitle>
                  <CardDescription>Automated collection at 11:00 JST (02:00 UTC)</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-4">
                    {externalSources.map((source) => (
                      <div key={source.name} className="flex items-start gap-4 p-4 border rounded-lg">
                        <div className="flex-shrink-0">
                          {source.method === "REST API" ? (
                            <DatabaseIcon className="h-5 w-5 text-primary" />
                          ) : source.method === "Web Scraping" ? (
                            <GlobeIcon className="h-5 w-5 text-blue-600 dark:text-blue-400" />
                          ) : (
                            <SearchIcon className="h-5 w-5 text-green-600 dark:text-green-400" />
                          )}
                        </div>
                        <div className="flex-1">
                          <div className="flex items-center justify-between mb-1">
                            <h3 className="font-semibold">{source.name}</h3>
                            <Badge variant="outline">
                              <ClockIcon className="h-3 w-3 mr-1" />
                              {source.frequency}
                            </Badge>
                          </div>
                          <p className="text-sm text-muted-foreground mb-2">{source.fullName}</p>
                          <div className="grid grid-cols-3 gap-2 text-xs">
                            <div>
                              <span className="text-muted-foreground">URL:</span>{" "}
                              <span className="font-mono">{source.url}</span>
                            </div>
                            <div>
                              <span className="text-muted-foreground">Method:</span>{" "}
                              <span>{source.method}</span>
                            </div>
                            <div>
                              <span className="text-muted-foreground">Data:</span>{" "}
                              <span>{source.dataType}</span>
                            </div>
                          </div>
                        </div>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Severity Classification */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Severity Classification</h2>

              <div className="space-y-3">
                {severityLevels.map((severity) => (
                  <Card key={severity.level}>
                    <CardHeader className={`${severity.color} rounded-t-lg`}>
                      <div className="flex items-center justify-between">
                        <CardTitle className="text-lg">{severity.level}</CardTitle>
                        <Badge variant="outline" className="bg-background">
                          Response: {severity.response}
                        </Badge>
                      </div>
                    </CardHeader>
                    <CardContent className="pt-4">
                      <div className="space-y-3">
                        <div>
                          <p className="text-sm font-medium mb-1">Criteria:</p>
                          <p className="text-sm text-muted-foreground">{severity.criteria}</p>
                        </div>
                        <div>
                          <p className="text-sm font-medium mb-2">Actions:</p>
                          <div className="flex flex-wrap gap-2">
                            {severity.actions.map((action) => (
                              <Badge key={action} variant="secondary" className="text-xs">
                                {action}
                              </Badge>
                            ))}
                          </div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                ))}
              </div>
            </section>

            {/* System Architecture */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">System Architecture</h2>

              <Card>
                <CardHeader>
                  <CardTitle>Data Flow</CardTitle>
                  <CardDescription>Dual-source monitoring with severity-based alerting</CardDescription>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    language="plaintext"
                    code={`External Sources (Daily)        Internal Pipeline (Real-time)
      ↓                               ↓
  Lambda Collector              S3 Event Trigger
      ↓                               ↓
  DynamoDB ←──────────────────── Lambda Listener
      ↓
Severity Engine
      ↓
Notification Router
   ├─ SNS/SES (Email/SMS)
   ├─ Streamlit Dashboard
   └─ REST API`}
                  />
                </CardContent>
              </Card>
            </section>

            {/* Components */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">System Components</h2>

              <div className="grid md:grid-cols-3 gap-4">
                <Card>
                  <CardHeader>
                    <BellIcon className="h-8 w-8 text-primary mb-2" />
                    <CardTitle className="text-base">Alerting</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm space-y-1">
                    <p>• AWS SNS Topics</p>
                    <p>• SES Email Templates</p>
                    <p>• SMS for Critical</p>
                    <p>• Deduplication (1h)</p>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <BarChartIcon className="h-8 w-8 text-primary mb-2" />
                    <CardTitle className="text-base">Dashboard</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm space-y-1">
                    <p>• Streamlit UI</p>
                    <p>• Real-time Updates (30s)</p>
                    <p>• Plotly Charts</p>
                    <p>• Trend Analysis</p>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <DatabaseIcon className="h-8 w-8 text-primary mb-2" />
                    <CardTitle className="text-base">Storage</CardTitle>
                  </CardHeader>
                  <CardContent className="text-sm space-y-1">
                    <p>• DynamoDB Tables (3)</p>
                    <p>• S3 Data Lake</p>
                    <p>• TTL: 90 days</p>
                    <p>• Point-in-Time Recovery</p>
                  </CardContent>
                </Card>
              </div>
            </section>

            {/* Usage Examples */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Quick Start</h2>

              <Card className="mb-4">
                <CardHeader>
                  <CardTitle>Launch Dashboard</CardTitle>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    language="bash"
                    code={`streamlit run surveillance/dashboard/app.py
# Access: http://localhost:8501`}
                  />
                </CardContent>
              </Card>

              <Card className="mb-4">
                <CardHeader>
                  <CardTitle>Start REST API</CardTitle>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    language="bash"
                    code={`uvicorn surveillance.api.main:app --reload --port 8000
# API Docs: http://localhost:8000/docs`}
                  />
                </CardContent>
              </Card>

              <Card>
                <CardHeader>
                  <CardTitle>Manual Collection Test</CardTitle>
                </CardHeader>
                <CardContent>
                  <CodeBlock
                    language="bash"
                    code={`# Test PubMed + J-STAGE search
python surveillance/external/academic_monitor.py

# Test MAFF scraping
python -m surveillance.external.maff_scraper

# Test E-Stat API
python -m surveillance.external.estat_client`}
                  />
                </CardContent>
              </Card>
            </section>

            {/* API Reference */}
            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">API Endpoints</h2>

              <Card>
                <CardHeader>
                  <CardTitle>REST API (FastAPI)</CardTitle>
                  <CardDescription>Programmatic access to surveillance data</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-3 text-sm font-mono">
                    <div className="flex items-start gap-3">
                      <Badge variant="outline" className="bg-blue-50 dark:bg-blue-900/10">GET</Badge>
                      <div>
                        <div className="font-semibold">/api/v1/detections</div>
                        <div className="text-muted-foreground text-xs">Get virus detections (filterable)</div>
                      </div>
                    </div>
                    <div className="flex items-start gap-3">
                      <Badge variant="outline" className="bg-blue-50 dark:bg-blue-900/10">GET</Badge>
                      <div>
                        <div className="font-semibold">/api/v1/alerts/active</div>
                        <div className="text-muted-foreground text-xs">Get active alerts summary</div>
                      </div>
                    </div>
                    <div className="flex items-start gap-3">
                      <Badge variant="outline" className="bg-blue-50 dark:bg-blue-900/10">GET</Badge>
                      <div>
                        <div className="font-semibold">/api/v1/external/daily-updates</div>
                        <div className="text-muted-foreground text-xs">Get external source updates</div>
                      </div>
                    </div>
                    <div className="flex items-start gap-3">
                      <Badge variant="outline" className="bg-blue-50 dark:bg-blue-900/10">GET</Badge>
                      <div>
                        <div className="font-semibold">/api/v1/statistics/trends</div>
                        <div className="text-muted-foreground text-xs">Get detection trends</div>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>
            </section>

            {/* Important Notes */}
            <section>
              <h2 className="text-2xl font-semibold mb-6">Important Notes</h2>

              <Alert variant="destructive" className="mb-4">
                <AlertTriangleIcon className="h-4 w-4" />
                <AlertTitle>EEEV Detection</AlertTitle>
                <AlertDescription>
                  Eastern Equine Encephalitis Virus is NOT endemic to Japan/Asia. Any detection requires immediate sequence confirmation and verification of sample origin.
                </AlertDescription>
              </Alert>

              <Alert className="mb-4">
                <AlertTriangleIcon className="h-4 w-4" />
                <AlertTitle>Hantavirus Context</AlertTitle>
                <AlertDescription>
                  Hantavirus is primarily rodent-borne. Pig detection suggests facility rodent contamination rather than pig infection. Focus on biosecurity and rodent control.
                </AlertDescription>
              </Alert>

              <Alert>
                <CheckCircleIcon className="h-4 w-4" />
                <AlertTitle>Spumavirus Priority</AlertTitle>
                <AlertDescription>
                  Spumavirus is designated as MHLW Special Management Pathogen #5. No porcine reference genome exists - detection uses cross-genus comparison.
                </AlertDescription>
              </Alert>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
