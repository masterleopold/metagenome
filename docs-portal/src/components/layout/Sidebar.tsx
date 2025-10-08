'use client'

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import {
  RocketIcon,
  LayersIcon,
  Code2Icon,
  GitBranchIcon,
  CloudIcon,
  ShieldCheckIcon,
  BookOpenIcon,
} from "lucide-react"

const navigation = [
  {
    title: "Introduction",
    items: [
      {
        title: "Getting Started",
        href: "/getting-started",
        icon: RocketIcon,
      },
      {
        title: "Overview",
        href: "/overview",
        icon: BookOpenIcon,
      },
    ],
  },
  {
    title: "Core Concepts",
    items: [
      {
        title: "Architecture",
        href: "/architecture",
        icon: LayersIcon,
      },
      {
        title: "Pipeline Phases",
        href: "/pipeline-phases",
        icon: GitBranchIcon,
      },
    ],
  },
  {
    title: "Development",
    items: [
      {
        title: "API Reference",
        href: "/api-reference",
        icon: Code2Icon,
      },
      {
        title: "Deployment",
        href: "/deployment",
        icon: CloudIcon,
      },
    ],
  },
  {
    title: "Compliance",
    items: [
      {
        title: "PMDA Guidelines",
        href: "/pmda-compliance",
        icon: ShieldCheckIcon,
      },
    ],
  },
]

export function Sidebar() {
  const pathname = usePathname()

  return (
    <aside className="fixed top-16 z-30 hidden h-[calc(100vh-4rem)] w-64 shrink-0 border-r bg-background md:sticky md:block">
      <div className="h-full overflow-y-auto py-8 px-4">
        <nav className="space-y-6">
          {navigation.map((section) => (
            <div key={section.title}>
              <h4 className="mb-2 px-3 text-xs font-semibold text-foreground/60 uppercase tracking-wider">
                {section.title}
              </h4>
              <div className="space-y-1">
                {section.items.map((item) => {
                  const Icon = item.icon
                  const isActive = pathname === item.href
                  return (
                    <Link
                      key={item.href}
                      href={item.href}
                      className={cn(
                        "flex items-center gap-3 rounded-lg px-3 py-2 text-sm font-medium transition-colors",
                        isActive
                          ? "bg-primary/10 text-primary"
                          : "text-foreground/60 hover:bg-accent hover:text-foreground"
                      )}
                    >
                      <Icon className="h-4 w-4" />
                      {item.title}
                    </Link>
                  )
                })}
              </div>
            </div>
          ))}
        </nav>
      </div>
    </aside>
  )
}
