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
} from "lucide-react";
import { Button } from "@/components/ui/Button";

const features = [
  {
    title: "PMDA Compliance",
    description: "Full coverage of 91 designated pathogens for xenotransplantation safety",
    icon: ShieldCheckIcon,
  },
  {
    title: "PERV Detection",
    description: "Critical detection of Porcine Endogenous Retroviruses (PERV-A, B, C)",
    icon: VibrateIcon,
  },
  {
    title: "Real-time Analysis",
    description: "Streaming analysis capability with MinION sequencing",
    icon: ActivityIcon,
  },
  {
    title: "Cloud-Native",
    description: "Serverless architecture on AWS (Lambda + EC2 on-demand)",
    icon: CloudIcon,
  },
  {
    title: "Automated Workflow",
    description: "End-to-end automation from basecalling to reporting",
    icon: WorkflowIcon,
  },
  {
    title: "Quality Assured",
    description: "Q30+ accuracy with duplex basecalling",
    icon: CheckCircleIcon,
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
    title: "API Reference",
    description: "Complete API documentation with request/response examples",
    href: "/api-reference" as const,
    icon: DatabaseIcon,
  },
];

const stats = [
  { label: "PMDA Pathogens", value: "91" },
  { label: "Pipeline Phases", value: "6" },
  { label: "AWS Services", value: "10+" },
  { label: "Accuracy (Q30+)", value: "99.9%" },
];

export default function Home() {
  return (
    <div className="flex flex-col">
      {/* Hero Section */}
      <section className="relative overflow-hidden">
        <div className="max-w-6xl mx-auto px-6 py-20 md:py-28 lg:py-32">
          <div className="mx-auto max-w-3xl text-center">
            <div className="inline-flex items-center rounded-full border px-3 py-1 text-sm mb-8">
              Version 1.0.0
            </div>
            <h1 className="text-5xl md:text-6xl lg:text-7xl font-bold tracking-tight mb-6" style={{ letterSpacing: '-0.03em' }}>
              MinION Pipeline
            </h1>
            <p className="text-lg md:text-xl leading-relaxed mb-10" style={{ color: 'hsl(var(--text-secondary))' }}>
              PMDA-compliant metagenomic analysis pipeline for xenotransplantation donor pig screening using Oxford Nanopore MinION Mk1D
            </p>
            <div className="flex flex-wrap gap-3 justify-center">
              <Button size="lg" asChild>
                <Link href="/getting-started" className="gap-2">
                  Get Started
                  <ArrowRightIcon className="h-4 w-4" />
                </Link>
              </Button>
              <Button size="lg" variant="outline" asChild>
                <Link href="/api-reference">View API Docs</Link>
              </Button>
            </div>
          </div>
        </div>
      </section>

      {/* Stats Section */}
      <section className="border-y" style={{ backgroundColor: 'hsl(var(--muted) / 0.3)' }}>
        <div className="max-w-6xl mx-auto px-6 py-16">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-8 md:gap-12">
            {stats.map((stat) => (
              <div key={stat.label} className="text-center">
                <div className="text-4xl md:text-5xl font-bold mb-2">{stat.value}</div>
                <div className="text-sm md:text-base" style={{ color: 'hsl(var(--text-tertiary))' }}>{stat.label}</div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section className="max-w-6xl mx-auto px-6 py-20 md:py-28">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold mb-4" style={{ letterSpacing: '-0.02em' }}>
            Key Features
          </h2>
          <p className="text-lg max-w-2xl mx-auto" style={{ color: 'hsl(var(--text-secondary))' }}>
            Comprehensive pathogen detection with enterprise-grade reliability
          </p>
        </div>
        <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-8 max-w-6xl mx-auto">
          {features.map((feature) => {
            const Icon = feature.icon;
            return (
              <div key={feature.title} className="group">
                <div className="mb-4">
                  <div className="inline-flex items-center justify-center w-10 h-10 rounded-lg transition-colors"
                       style={{ backgroundColor: 'hsl(var(--muted))' }}>
                    <Icon className="h-5 w-5" style={{ color: 'hsl(var(--foreground))' }} />
                  </div>
                </div>
                <h3 className="text-lg font-semibold mb-2">{feature.title}</h3>
                <p className="text-sm leading-relaxed" style={{ color: 'hsl(var(--text-secondary))' }}>
                  {feature.description}
                </p>
              </div>
            );
          })}
        </div>
      </section>

      {/* Quick Links Section */}
      <section className="border-y" style={{ backgroundColor: 'hsl(var(--muted) / 0.3)' }}>
        <div className="max-w-6xl mx-auto px-6 py-20 md:py-28">
          <div className="text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4" style={{ letterSpacing: '-0.02em' }}>
              Documentation
            </h2>
            <p className="text-lg max-w-2xl mx-auto" style={{ color: 'hsl(var(--text-secondary))' }}>
              Everything you need to get up and running
            </p>
          </div>
          <div className="grid md:grid-cols-3 gap-6 max-w-5xl mx-auto">
            {quickLinks.map((link) => {
              const Icon = link.icon;
              return (
                <Link
                  key={link.href}
                  href={link.href}
                  className="group block p-6 rounded-lg border transition-all hover:border-[hsl(var(--foreground))] hover:shadow-sm"
                  style={{ backgroundColor: 'hsl(var(--background))' }}
                >
                  <div className="flex items-center justify-between mb-4">
                    <Icon className="h-6 w-6" />
                    <ArrowRightIcon
                      className="h-4 w-4 transition-transform group-hover:translate-x-1"
                      style={{ color: 'hsl(var(--text-tertiary))' }}
                    />
                  </div>
                  <h3 className="text-lg font-semibold mb-2">{link.title}</h3>
                  <p className="text-sm leading-relaxed" style={{ color: 'hsl(var(--text-secondary))' }}>
                    {link.description}
                  </p>
                </Link>
              );
            })}
          </div>
        </div>
      </section>

      {/* Tech Stack Section */}
      <section className="max-w-6xl mx-auto px-6 py-20 md:py-28">
        <div className="text-center mb-16">
          <h2 className="text-3xl md:text-4xl font-bold mb-4" style={{ letterSpacing: '-0.02em' }}>
            Technology Stack
          </h2>
          <p className="text-lg max-w-2xl mx-auto" style={{ color: 'hsl(var(--text-secondary))' }}>
            Built with modern, production-ready technologies
          </p>
        </div>
        <div className="max-w-4xl mx-auto">
          <div className="grid md:grid-cols-2 gap-12">
            <div>
              <h3 className="text-xl font-semibold mb-6">Cloud & Infrastructure</h3>
              <ul className="space-y-3">
                {[
                  "AWS Lambda + EC2 on-demand",
                  "Step Functions workflow orchestration",
                  "RDS Aurora Serverless PostgreSQL",
                  "S3 + EFS for storage"
                ].map((item) => (
                  <li key={item} className="flex items-start gap-3">
                    <CheckCircleIcon className="h-5 w-5 mt-0.5 flex-shrink-0" style={{ color: 'hsl(var(--primary))' }} />
                    <span style={{ color: 'hsl(var(--text-secondary))' }}>{item}</span>
                  </li>
                ))}
              </ul>
            </div>
            <div>
              <h3 className="text-xl font-semibold mb-6">Analysis Tools</h3>
              <ul className="space-y-3">
                {[
                  "Oxford Nanopore Dorado (Duplex)",
                  "Kraken2, BLAST, Diamond",
                  "Minimap2, SAMtools",
                  "Python 3.9+ analysis scripts"
                ].map((item) => (
                  <li key={item} className="flex items-start gap-3">
                    <CheckCircleIcon className="h-5 w-5 mt-0.5 flex-shrink-0" style={{ color: 'hsl(var(--primary))' }} />
                    <span style={{ color: 'hsl(var(--text-secondary))' }}>{item}</span>
                  </li>
                ))}
              </ul>
            </div>
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="border-t">
        <div className="max-w-6xl mx-auto px-6 py-20 md:py-28">
          <div className="max-w-2xl mx-auto text-center">
            <h2 className="text-3xl md:text-4xl font-bold mb-4" style={{ letterSpacing: '-0.02em' }}>
              Ready to get started?
            </h2>
            <p className="text-lg mb-10" style={{ color: 'hsl(var(--text-secondary))' }}>
              Follow our comprehensive guide to deploy your own MinION pathogen screening pipeline
            </p>
            <div className="flex flex-wrap gap-3 justify-center">
              <Button size="lg" asChild>
                <Link href="/getting-started" className="gap-2">
                  Read the Documentation
                  <ArrowRightIcon className="h-4 w-4" />
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
