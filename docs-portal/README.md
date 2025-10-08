# MinION Pipeline Documentation Portal

Modern, interactive documentation portal for the MinION Pathogen Screening Pipeline built with Next.js 15, React 19, and Untitled UI components.

## Features

- 📚 **Comprehensive Documentation**: Complete guides for setup, deployment, and API usage
- 🎨 **Modern UI**: Built with Untitled UI React components and Tailwind CSS 4
- 🌙 **Dark Mode**: Full dark mode support with system preference detection
- 📱 **Responsive**: Mobile-first design that works on all devices
- ⚡ **Fast**: Server-side rendering with Next.js 15 App Router
- 🔍 **Search**: Fast documentation search (coming soon)
- 💻 **Interactive Code**: Copy-to-clipboard code examples
- ♿ **Accessible**: WCAG 2.1 AA compliant with React Aria

## Tech Stack

- **Framework**: Next.js 15 (App Router)
- **UI Library**: React 19
- **Styling**: Tailwind CSS 4
- **Components**: Untitled UI React + Radix UI
- **Icons**: Lucide React
- **Theme**: next-themes
- **Typography**: Geist Sans & Geist Mono
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
- Key features showcase
- Quick links to documentation sections
- Technology stack overview

✅ **Getting Started** (`/getting-started`)
- Prerequisites checklist
- Installation steps
- Infrastructure deployment guide
- First workflow tutorial

✅ **Overview** (`/overview`)
- Research background and context
- Sample processing workflow
- Regulatory framework
- Technical approach

✅ **Architecture** (`/architecture`)
- System architecture diagram
- Component breakdown
- Data flow visualization
- AWS services overview

✅ **Pipeline Phases** (`/pipeline-phases`)
- 6 pipeline phases explained
- Tool documentation
- Performance benchmarks
- Best practices

✅ **API Reference** (`/api-reference`)
- All API endpoints
- Request/response examples
- Authentication guide
- Error codes reference

✅ **Deployment** (`/deployment`)
- Step-by-step deployment
- Configuration options
- Troubleshooting guide
- Cost optimization tips

✅ **PMDA Compliance** (`/pmda-compliance`)
- 91 pathogen list
- PERV detection details
- Validation requirements
- Reporting formats

## Component Usage

### Untitled UI Components

This project uses Untitled UI React components. Key components:

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

The portal uses `#0089A7` (HSL: 192 100% 33%) as the primary brand color. Edit `src/styles/globals.css` to customize:

```css
@theme {
  --color-primary: #0089A7;
  --color-secondary: oklch(0.65 0.15 180);
}
```

### Navigation

Edit navigation in `src/components/layout/Sidebar.tsx` and `src/components/layout/Header.tsx`:

**Important**: When using Next.js `typedRoutes`, all href values must be const-asserted:

```ts
const navigation = [
  {
    title: 'Getting Started',
    href: '/getting-started' as const,  // Required for typed routes
    icon: RocketIcon,
  },
  // ...
] as const  // Required for typed routes
```

## Deployment

### Vercel (Recommended)

1. Push your code to GitHub
2. Import project in Vercel
3. Deploy with zero configuration

**Note**: The project is configured for Vercel deployment with:
- Next.js 15.5.4 with App Router
- TypeScript typed routes enabled (`typedRoutes: true`)
- Automatic static page optimization
- Edge runtime support

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
