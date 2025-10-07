import { type ClassValue, clsx } from "clsx"
import { twMerge } from "tailwind-merge"

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs))
}

export function formatDate(date: string | Date): string {
  const d = typeof date === 'string' ? new Date(date) : date
  return d.toLocaleDateString('en-US', {
    year: 'numeric',
    month: 'long',
    day: 'numeric'
  })
}

export function copyToClipboard(text: string): Promise<void> {
  if (navigator.clipboard) {
    return navigator.clipboard.writeText(text)
  }
  // Fallback for older browsers
  const textarea = document.createElement('textarea')
  textarea.value = text
  textarea.style.position = 'fixed'
  textarea.style.opacity = '0'
  document.body.appendChild(textarea)
  textarea.select()
  document.execCommand('copy')
  document.body.removeChild(textarea)
  return Promise.resolve()
}

export function slugify(text: string): string {
  return text
    .toLowerCase()
    .replace(/[^\w\s-]/g, '')
    .replace(/[\s_-]+/g, '-')
    .replace(/^-+|-+$/g, '')
}

export const navigation = [
  {
    title: 'Getting Started',
    href: '/getting-started',
    icon: 'RocketLaunch',
  },
  {
    title: 'Architecture',
    href: '/architecture',
    icon: 'Layers',
  },
  {
    title: 'API Reference',
    href: '/api-reference',
    icon: 'Code',
  },
  {
    title: 'Pipeline Phases',
    href: '/pipeline-phases',
    icon: 'GitBranch',
  },
  {
    title: 'Deployment',
    href: '/deployment',
    icon: 'Cloud',
  },
  {
    title: 'PMDA Compliance',
    href: '/pmda-compliance',
    icon: 'ShieldCheck',
  },
]

export const features = [
  {
    title: 'PMDA Compliance',
    description: 'Full coverage of 91 designated pathogens for xenotransplantation safety',
    icon: 'ShieldCheck',
    color: 'text-primary-600',
  },
  {
    title: 'PERV Detection',
    description: 'Critical detection of Porcine Endogenous Retroviruses (PERV-A, B, C)',
    icon: 'Virus',
    color: 'text-red-600',
  },
  {
    title: 'Real-time Analysis',
    description: 'Streaming analysis capability with MinION sequencing',
    icon: 'Activity',
    color: 'text-secondary-600',
  },
  {
    title: 'Cloud-Native',
    description: 'Serverless architecture on AWS (Lambda + EC2 on-demand)',
    icon: 'Cloud',
    color: 'text-blue-600',
  },
  {
    title: 'Automated Workflow',
    description: 'End-to-end automation from basecalling to reporting',
    icon: 'Workflow',
    color: 'text-purple-600',
  },
  {
    title: 'Quality Assured',
    description: 'Q30+ accuracy with duplex basecalling',
    icon: 'CheckCircle',
    color: 'text-green-600',
  },
]

export const pipelinePhases = [
  {
    id: 1,
    name: 'Basecalling',
    description: 'Convert raw signal (FAST5/POD5) to sequences (FASTQ)',
    tool: 'Dorado Duplex',
    duration: '2-4 hours',
    instance: 'g4dn.xlarge (GPU)',
  },
  {
    id: 2,
    name: 'Quality Control',
    description: 'Assess read quality metrics and filter low-quality reads',
    tool: 'PycoQC, NanoPlot',
    duration: '10-15 minutes',
    instance: 't3.large',
  },
  {
    id: 3,
    name: 'Host Removal',
    description: 'Align reads to Sus scrofa genome and remove host DNA',
    tool: 'Minimap2, SAMtools',
    duration: '30-60 minutes',
    instance: 'r5.xlarge',
  },
  {
    id: 4,
    name: 'Pathogen Detection',
    description: 'Multi-database screening for PMDA pathogens',
    tool: 'Kraken2, BLAST, Diamond',
    duration: '1-2 hours',
    instance: 'r5.4xlarge',
  },
  {
    id: 5,
    name: 'Quantification',
    description: 'Absolute copy number calculation with spike-in normalization',
    tool: 'Custom Python scripts',
    duration: '15-30 minutes',
    instance: 't3.large',
  },
  {
    id: 6,
    name: 'Reporting',
    description: 'Generate PMDA-compliant reports in PDF, JSON, and HTML',
    tool: 'ReportLab, WeasyPrint',
    duration: '10-15 minutes',
    instance: 't3.medium',
  },
]

export const apiEndpoints = [
  {
    method: 'POST',
    path: '/workflows',
    description: 'Start a new analysis workflow',
    requiresAuth: true,
  },
  {
    method: 'GET',
    path: '/workflows/{workflow_id}',
    description: 'Get workflow execution status',
    requiresAuth: true,
  },
  {
    method: 'GET',
    path: '/workflows',
    description: 'List all workflow executions',
    requiresAuth: true,
  },
  {
    method: 'DELETE',
    path: '/workflows/{workflow_id}',
    description: 'Stop a running workflow',
    requiresAuth: true,
  },
  {
    method: 'GET',
    path: '/pathogens/{run_id}',
    description: 'Get pathogen detection results',
    requiresAuth: true,
  },
  {
    method: 'GET',
    path: '/reports/{run_id}',
    description: 'Get analysis report URLs',
    requiresAuth: true,
  },
]
