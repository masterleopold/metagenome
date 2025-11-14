import Link from "next/link";
import {
  ShieldCheckIcon,
  VibrateIcon,
  ActivityIcon,
  CloudIcon,
  WorkflowIcon,
  CheckCircleIcon,
  ArrowRightIcon,
  DatabaseIcon,
  ServerIcon,
  ZapIcon,
  RadarIcon,
} from "lucide-react";
import { Button } from "@/components/ui/Button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";

const features = [
  {
    title: "PMDA Compliance",
    description: "Full coverage of 91 designated pathogens for xenotransplantation safety",
    icon: ShieldCheckIcon,
    color: "text-primary",
  },
  {
    title: "PERV Detection",
    description: "Critical detection of Porcine Endogenous Retroviruses (PERV-A, B, C)",
    icon: VibrateIcon,
    color: "text-red-600",
  },
  {
    title: "Real-time Analysis",
    description: "Streaming analysis capability with MinION sequencing",
    icon: ActivityIcon,
    color: "text-secondary",
  },
  {
    title: "Cloud-Native",
    description: "Serverless architecture on AWS (Lambda + EC2 on-demand)",
    icon: CloudIcon,
    color: "text-primary",
  },
  {
    title: "Automated Workflow",
    description: "End-to-end automation from basecalling to reporting",
    icon: WorkflowIcon,
    color: "text-secondary",
  },
  {
    title: "Quality Assured",
    description: "Q30+ accuracy with duplex basecalling",
    icon: CheckCircleIcon,
    color: "text-green-600",
  },
];

const quickLinks = [
  {
    title: "Getting Started",
    description: "Set up your development environment and run your first workflow",
    href: "/getting-started" as const,
    icon: ZapIcon,
  },
  {
    title: "Architecture",
    description: "Understand the system architecture and component interactions",
    href: "/architecture" as const,
    icon: ServerIcon,
  },
  {
    title: "4-Virus Surveillance",
    description: "Real-time monitoring of Hantavirus, Polyomavirus, Spumavirus, and EEEV",
    href: "/surveillance" as const,
    icon: RadarIcon,
  },
  {
    title: "API Reference",
    description: "Complete API documentation with request/response examples",
    href: "/api-reference" as const,
    icon: DatabaseIcon,
  },
] as const;

const stats = [
  { label: "PMDA Pathogens", value: "91" },
  { label: "Surveillance Viruses", value: "4" },
  { label: "Pipeline Phases", value: "6" },
  { label: "AWS Services", value: "10+" },
  { label: "Accuracy (Q30+)", value: "99.9%" },
];

export default function Home() {
  return (
    <div className="flex flex-col">
      {/* Hero Section */}
      <section className="relative overflow-hidden border-b bg-gradient-to-b from-background to-muted/20">
        <div className="container mx-auto max-w-7xl px-4 py-24 md:py-32">
          <div className="mx-auto max-w-4xl text-center">
            <Badge variant="secondary" className="mb-4">
              Version 1.0.0
            </Badge>
            <h1 className="text-4xl font-bold tracking-tight sm:text-6xl mb-6 bg-gradient-to-r from-primary via-secondary to-primary bg-clip-text text-transparent">
              MinION Pathogen Screening Pipeline
            </h1>
            <p className="text-xl text-muted-foreground mb-8 leading-relaxed">
              PMDA-compliant metagenomic analysis pipeline for xenotransplantation donor pig screening using Oxford Nanopore MinION Mk1D
            </p>
            <div className="flex flex-wrap gap-4 justify-center">
              <Button size="lg" asChild>
                <Link href="/getting-started" className="inline-flex items-center">
                  <span>Get Started</span>
                  <ArrowRightIcon className="ml-2 h-4 w-4" />
                </Link>
              </Button>
              <Button size="lg" variant="outline" asChild>
                <Link href="/api-reference">View API Docs</Link>
              </Button>
            </div>
          </div>
        </div>

        {/* Animated gradient blob */}
        <div className="absolute top-0 left-1/2 -translate-x-1/2 -translate-y-1/2 w-[800px] h-[800px] bg-gradient-to-r from-primary/20 to-secondary/20 rounded-full blur-3xl opacity-20 pointer-events-none" />
      </section>

      {/* Stats Section */}
      <section className="border-b bg-muted/30">
        <div className="container mx-auto max-w-7xl px-4 py-12">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-8">
            {stats.map((stat) => (
              <div key={stat.label} className="text-center">
                <div className="text-4xl font-bold text-primary mb-2">{stat.value}</div>
                <div className="text-sm text-muted-foreground">{stat.label}</div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section className="container mx-auto max-w-7xl px-4 py-24">
        <div className="text-center mb-12">
          <h2 className="text-3xl font-bold mb-4">Key Features</h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Comprehensive pathogen detection with enterprise-grade reliability and compliance
          </p>
        </div>
        <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-6">
          {features.map((feature) => {
            const Icon = feature.icon;
            return (
              <Card key={feature.title} className="hover:shadow-lg transition-shadow">
                <CardHeader>
                  <div className={`w-12 h-12 rounded-lg bg-gradient-to-br from-primary/10 to-secondary/10 flex items-center justify-center mb-4`}>
                    <Icon className={`h-6 w-6 ${feature.color}`} />
                  </div>
                  <CardTitle>{feature.title}</CardTitle>
                </CardHeader>
                <CardContent>
                  <CardDescription className="text-base">
                    {feature.description}
                  </CardDescription>
                </CardContent>
              </Card>
            );
          })}
        </div>
      </section>

      {/* Quick Links Section */}
      <section className="bg-muted/30 border-y">
        <div className="container mx-auto max-w-7xl px-4 py-24">
          <div className="text-center mb-12">
            <h2 className="text-3xl font-bold mb-4">Quick Links</h2>
            <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
              Jump right into the documentation
            </p>
          </div>
          <div className="grid md:grid-cols-3 gap-6">
            {quickLinks.map((link) => {
              const Icon = link.icon;
              return (
                <Link key={link.href} href={link.href} className="group">
                  <Card className="h-full hover:shadow-lg transition-all hover:border-primary/50">
                    <CardHeader>
                      <div className="flex items-center justify-between mb-2">
                        <Icon className="h-8 w-8 text-primary group-hover:scale-110 transition-transform" />
                        <ArrowRightIcon className="h-5 w-5 text-muted-foreground group-hover:translate-x-1 transition-transform" />
                      </div>
                      <CardTitle>{link.title}</CardTitle>
                    </CardHeader>
                    <CardContent>
                      <CardDescription className="text-base">
                        {link.description}
                      </CardDescription>
                    </CardContent>
                  </Card>
                </Link>
              );
            })}
          </div>
        </div>
      </section>

      {/* Tech Stack Section */}
      <section className="container mx-auto max-w-7xl px-4 py-24">
        <div className="text-center mb-12">
          <h2 className="text-3xl font-bold mb-4">Technology Stack</h2>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Built with modern, production-ready technologies
          </p>
        </div>
        <div className="max-w-4xl mx-auto">
          <div className="grid md:grid-cols-2 gap-8">
            <div>
              <h3 className="text-xl font-semibold mb-4">Cloud & Infrastructure</h3>
              <ul className="space-y-2">
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>AWS Lambda + EC2 on-demand</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>Step Functions workflow orchestration</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>RDS Aurora Serverless PostgreSQL</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>S3 + EFS for storage</span>
                </li>
              </ul>
            </div>
            <div>
              <h3 className="text-xl font-semibold mb-4">Analysis Tools</h3>
              <ul className="space-y-2">
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>Oxford Nanopore Dorado (Duplex)</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>Kraken2, BLAST, Diamond</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>Minimap2, SAMtools</span>
                </li>
                <li className="flex items-center gap-2">
                  <CheckCircleIcon className="h-5 w-5 text-green-600" />
                  <span>Python 3.9+ analysis scripts</span>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="border-t">
        <div className="container mx-auto max-w-7xl px-4 py-24">
          <div className="max-w-3xl mx-auto text-center">
            <h2 className="text-3xl font-bold mb-4">Ready to get started?</h2>
            <p className="text-lg text-muted-foreground mb-8">
              Follow our comprehensive guide to deploy your own MinION pathogen screening pipeline
            </p>
            <div className="flex flex-wrap gap-4 justify-center">
              <Button size="lg" asChild>
                <Link href="/getting-started" className="inline-flex items-center">
                  <span>Read the Documentation</span>
                  <ArrowRightIcon className="ml-2 h-4 w-4" />
                </Link>
              </Button>
              <Button size="lg" variant="outline" asChild>
                <a href="https://github.com/masterleopold/metagenome" target="_blank" rel="noreferrer">
                  View on GitHub
                </a>
              </Button>
            </div>
          </div>
        </div>
      </section>
    </div>
  );
}
