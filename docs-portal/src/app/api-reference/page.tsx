import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import { InfoIcon } from "lucide-react";

export default function ApiReferencePage() {
  return (
    <div className="container flex-1">
      <div className="flex gap-6">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Development</Badge>
              <h1 className="text-4xl font-bold mb-4">API Reference</h1>
              <p className="text-lg text-muted-foreground">
                Complete API documentation for the MinION pathogen screening pipeline.
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
                The MinION Pipeline API provides REST endpoints for workflow management,
                status monitoring, and result retrieval. All endpoints require authentication
                via the <code className="bg-muted px-2 py-1 rounded text-sm">x-api-key</code> header.
              </p>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Base URL</h2>
              <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                https://api.minion-pipeline.com/{'{environment}'}/
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Endpoints</h2>
              <div className="space-y-6">
                <div className="border rounded-lg p-6 space-y-3">
                  <div className="flex items-center gap-3">
                    <span className="px-3 py-1 rounded font-mono text-sm font-semibold" style={{
                      backgroundColor: 'hsl(var(--primary) / 0.1)',
                      color: 'hsl(var(--primary))'
                    }}>
                      POST
                </span>
                <code className="text-lg font-mono">/workflows</code>
              </div>
                  <p className="text-muted-foreground">Start a new workflow execution</p>
                </div>

                <div className="border rounded-lg p-6 space-y-3">
                  <div className="flex items-center gap-3">
                    <span className="px-3 py-1 rounded font-mono text-sm font-semibold" style={{
                      backgroundColor: 'hsl(var(--primary) / 0.1)',
                      color: 'hsl(var(--primary))'
                    }}>
                      GET
                    </span>
                    <code className="text-lg font-mono">/workflows/{'{workflow_id}'}</code>
                  </div>
                  <p className="text-muted-foreground">Get workflow status and details</p>
                </div>

                <div className="border rounded-lg p-6 space-y-3">
                  <div className="flex items-center gap-3">
                    <span className="px-3 py-1 rounded font-mono text-sm font-semibold" style={{
                      backgroundColor: 'hsl(var(--primary) / 0.1)',
                      color: 'hsl(var(--primary))'
                    }}>
                      GET
                    </span>
                    <code className="text-lg font-mono">/workflows/{'{workflow_id}'}/metrics</code>
                  </div>
                  <p className="text-muted-foreground">Get performance metrics for a workflow</p>
                </div>

                <div className="border rounded-lg p-6 space-y-3">
                  <div className="flex items-center gap-3">
                    <span className="px-3 py-1 rounded font-mono text-sm font-semibold" style={{
                      backgroundColor: 'hsl(var(--primary) / 0.1)',
                      color: 'hsl(var(--primary))'
                    }}>
                      GET
                    </span>
                    <code className="text-lg font-mono">/workflows/{'{workflow_id}'}/results</code>
                  </div>
                  <p className="text-muted-foreground">Get analysis results and pathogen detections</p>
                </div>

                <div className="border rounded-lg p-6 space-y-3">
                  <div className="flex items-center gap-3">
                    <span className="px-3 py-1 rounded font-mono text-sm font-semibold" style={{
                      backgroundColor: 'hsl(var(--destructive) / 0.1)',
                      color: 'hsl(var(--destructive))'
                    }}>
                      DELETE
                    </span>
                    <code className="text-lg font-mono">/workflows/{'{workflow_id}'}</code>
                  </div>
                  <p className="text-muted-foreground">Cancel a running workflow</p>
                </div>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-4">Authentication</h2>
              <p className="text-muted-foreground mb-4">
                All API requests require an API key passed in the request header:
              </p>
              <div className="bg-muted p-4 rounded-lg font-mono text-sm">
                x-api-key: YOUR_API_KEY
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  )
}
