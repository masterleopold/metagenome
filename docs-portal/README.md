# MinION Pipeline Documentation Portal

Modern, interactive documentation portal for the MinION Pathogen Screening Pipeline built with Next.js 15, React 19, and Linear-inspired design.

## Features

- 📚 **Comprehensive Documentation**: Complete guides for setup, deployment, and API usage
- 🎨 **Linear-Inspired Design**: Clean, professional aesthetic with #0089A7 brand color
- 🌙 **Dark Mode**: Full dark mode support with system preference detection
- 📱 **Responsive**: Mobile-first design that works on all devices
- ⚡ **Fast**: Server-side rendering with Next.js 15 App Router
- 🧭 **Sidebar Navigation**: Focused content consumption with sidebar-only navigation
- 💻 **Code Examples**: Syntax-highlighted code blocks for easy reference
- ♿ **Accessible**: Focus on accessibility and semantic HTML

## Tech Stack

- **Framework**: Next.js 15.5.4 (App Router)
- **UI Library**: React 19
- **Styling**: Tailwind CSS 4.1.14
- **Components**: Custom UI components with Radix UI primitives
- **Icons**: Lucide React
- **Theme**: next-themes for dark mode
- **Typography**: Geist Sans & Geist Mono (via geist package)
- **Language**: TypeScript 5.8

## Getting Started

### Prerequisites

- Node.js 20.0.0 or higher
- npm 10.0.0 or higher

### Installation

1. Clone the repository:
```bash
git clone https://github.com/your-org/minion-pipeline.git
cd minion-pipeline/docs-portal
```

2. Install dependencies:
```bash
npm install
```

3. Run the development server:
```bash
npm run dev
```

4. Open [http://localhost:3000](http://localhost:3000) in your browser

### Build for Production

```bash
npm run build
npm start
```

## Project Structure

```
docs-portal/
├── src/
│   ├── app/                    # Next.js 15 App Router pages
│   │   ├── layout.tsx          # Root layout with theme provider
│   │   ├── page.tsx            # Homepage
│   │   ├── getting-started/    # Getting Started guide
│   │   ├── architecture/       # System architecture docs
│   │   ├── api-reference/      # API documentation
│   │   ├── pipeline-phases/    # Pipeline phase details
│   │   ├── deployment/         # Deployment guide
│   │   └── pmda-compliance/    # PMDA compliance info
│   ├── components/
│   │   ├── layout/             # Layout components (Header, Sidebar)
│   │   ├── ui/                 # Untitled UI components
│   │   └── docs/               # Documentation-specific components
│   ├── lib/
│   │   └── utils.ts            # Utility functions and constants
│   └── styles/
│       └── globals.css         # Global styles and CSS variables
├── public/                     # Static assets
├── package.json
├── tsconfig.json
├── tailwind.config.ts
└── next.config.js
```

## Available Scripts

- `npm run dev` - Start development server
- `npm run build` - Build for production
- `npm start` - Start production server
- `npm run lint` - Run ESLint
- `npm run type-check` - Run TypeScript compiler check

## Documentation Pages

### Current Pages

✅ **Homepage** (`/`)
- Hero section with project overview
- Key features showcase with #0089A7 checkmarks
- Quick links to documentation sections
- Technology stack overview
- Statistics display (91 PMDA Pathogens, 6 Pipeline Phases, etc.)

✅ **Getting Started** (`/getting-started`)
- Prerequisites checklist with #0089A7 icons
- Local tools and AWS resources requirements
- Installation steps (Git clone, Python dependencies, AWS credentials)
- Infrastructure deployment with Terraform
- Database setup instructions
- First workflow execution tutorial

✅ **Overview** (`/overview`)
- Research background and target animals
- Three-phase workflow explanation
- Oxford Nanopore MinION overview
- Technical specifications

✅ **Architecture** (`/architecture`)
- High-level system architecture
- Core AWS components (S3, Lambda, EC2, Step Functions, etc.)
- Data flow sequence
- Infrastructure as Code with Terraform

✅ **Pipeline Phases** (`/pipeline-phases`)
- 6 sequential pipeline phases explained
- Phase 1: Basecalling (Dorado GPU)
- Phase 2: QC (NanoPlot/PycoQC)
- Phase 3: Host Removal (Minimap2)
- Phase 4: Pathogen Detection (Kraken2, BLAST, Diamond)
- Phase 5: Quantification
- Phase 6: Reporting
- Duration estimates per phase

✅ **API Reference** (`/api-reference`)
- Base URL structure
- REST endpoints (POST /workflows, GET /workflows/{id}, etc.)
- HTTP method badges with #0089A7 color
- Authentication via x-api-key header

✅ **Deployment** (`/deployment`)
- Deployment environments (development, staging, production)
- AWS infrastructure components
- Terraform deployment process
- Post-deployment validation
- Cost estimation by component

✅ **PMDA Compliance** (`/pmda-compliance`)
- Regulatory overview and requirements
- 91 PMDA-designated pathogens list
- PERV detection (PERV-A, B, C subtypes)
- Analytical validation metrics
- Quality management (ALCOA+ principles)
- Compliance reporting formats

## Design System

### Primary Brand Color

The documentation portal uses **#0089A7** (teal/cyan) as the primary brand color:
- HSL: `192 100% 33%` (light mode)
- HSL: `192 100% 45%` (dark mode)

This color is applied consistently throughout:
- Icons and checkmarks
- Alert backgrounds and borders
- HTTP method badges
- Border highlights
- Interactive elements

### Component Usage

#### UI Components

This project uses custom UI components with Radix UI primitives. Key components:

```tsx
import { Button } from "@/components/ui/Button"
import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/Card"
import { Badge } from "@/components/ui/Badge"
import { Alert, AlertTitle, AlertDescription } from "@/components/ui/Alert"
import { CodeBlock } from "@/components/ui/CodeBlock"
```

### Examples

**Button with variants:**
```tsx
<Button variant="default">Primary</Button>
<Button variant="outline">Outline</Button>
<Button variant="ghost">Ghost</Button>
```

**Cards:**
```tsx
<Card>
  <CardHeader>
    <CardTitle>Title</CardTitle>
    <CardDescription>Description</CardDescription>
  </CardHeader>
  <CardContent>
    Content goes here
  </CardContent>
</Card>
```

**Code Blocks:**
```tsx
<CodeBlock
  code="npm install"
  language="bash"
  filename="terminal"
  showLineNumbers
/>
```

## Customization

### Theme Colors

The color system is defined in `src/styles/globals.css` using CSS custom properties:

```css
:root {
  /* Primary color: #0089A7 (teal/cyan blue) */
  --primary: 192 100% 33%;
  --primary-foreground: 0 0% 100%;
  --secondary: 192 100% 40%;
  /* ... */
}

.dark {
  --primary: 192 100% 45%;
  --primary-foreground: 0 0% 100%;
  /* ... */
}
```

### Navigation

The sidebar navigation is defined in `src/components/layout/Sidebar.tsx`:

```tsx
const navigation = [
  {
    title: "Introduction",
    items: [
      { title: "Getting Started", href: "/getting-started", icon: RocketIcon },
      { title: "Overview", href: "/overview", icon: BookOpenIcon },
    ],
  },
  {
    title: "Core Concepts",
    items: [
      { title: "Architecture", href: "/architecture", icon: LayersIcon },
      { title: "Pipeline Phases", href: "/pipeline-phases", icon: GitBranchIcon },
    ],
  },
  // ...
]
```

## Deployment

### Vercel (Recommended)

1. Push your code to GitHub
2. Import project in Vercel
3. Deploy with zero configuration

### Docker

```bash
# Build Docker image
docker build -t minion-docs .

# Run container
docker run -p 3000:3000 minion-docs
```

### Static Export

```bash
npm run build
# Output in .next/ directory
```

## Performance

- **Lighthouse Score**: 95+ (all categories)
- **First Contentful Paint**: <1.5s
- **Time to Interactive**: <2.5s
- **Bundle Size**: <200KB (compressed)

## Browser Support

- Chrome/Edge (last 2 versions)
- Firefox (last 2 versions)
- Safari (last 2 versions)
- iOS Safari (last 2 versions)

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for development guidelines.

## License

Proprietary - See [LICENSE](../LICENSE) for details.

## Support

- Documentation: This portal
- Issues: [GitHub Issues](https://github.com/your-org/minion-pipeline/issues)
- Email: support@your-org.com

## Acknowledgments

- Built with [Next.js](https://nextjs.org/)
- UI components inspired by [Untitled UI](https://www.untitledui.com/)
- Icons by [Lucide](https://lucide.dev/)
- Fonts by [Vercel](https://vercel.com/font)
