'use client'

import Link from "next/link"
import { usePathname } from "next/navigation"
import { Menu, X, Github, Moon, Sun } from "lucide-react"
import { Button } from "@/components/ui/Button"
import { cn } from "@/lib/utils"
import { useState } from "react"
import { useTheme } from "next-themes"

export function Header() {
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)
  const pathname = usePathname()
  const { theme, setTheme } = useTheme()

  const navigation = [
    { name: "Docs", href: "/getting-started" },
    { name: "API", href: "/api-reference" },
    { name: "Architecture", href: "/architecture" },
    { name: "Pipeline", href: "/pipeline-phases" },
  ]

  return (
    <header className="sticky top-0 z-50 w-full border-b bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60">
      <div className="container mx-auto max-w-7xl">
        <nav className="flex h-16 items-center justify-between px-4">
          <div className="flex items-center gap-6 md:gap-10">
            <Link href="/" className="flex items-center space-x-2">
              <div className="h-8 w-8 rounded-lg bg-gradient-to-br from-primary to-secondary flex items-center justify-center">
                <span className="text-white font-bold text-lg">M</span>
              </div>
              <span className="hidden font-bold sm:inline-block">
                MinION Pipeline
              </span>
            </Link>
            <div className="hidden md:flex gap-6">
              {navigation.map((item) => (
                <Link
                  key={item.href}
                  href={item.href}
                  className={cn(
                    "text-sm font-medium transition-colors hover:text-primary",
                    pathname?.startsWith(item.href)
                      ? "text-foreground"
                      : "text-foreground/60"
                  )}
                >
                  {item.name}
                </Link>
              ))}
            </div>
          </div>
          <div className="flex items-center gap-1">
            <Button
              variant="ghost"
              size="icon"
              onClick={() => {
                if (theme === "dark") {
                  setTheme("light")
                } else if (theme === "light") {
                  setTheme("dark")
                } else {
                  // If system, set to the opposite of current appearance
                  setTheme("dark")
                }
              }}
            >
              {theme === "dark" ? (
                <Sun className="h-5 w-5 transition-all" />
              ) : (
                <Moon className="h-5 w-5 transition-all" />
              )}
              <span className="sr-only">Toggle theme</span>
            </Button>
            <Button variant="ghost" size="icon" asChild className="hidden md:flex">
              <a
                href="https://github.com/masterleopold/metagenome"
                target="_blank"
                rel="noreferrer"
              >
                <Github className="h-5 w-5" />
                <span className="sr-only">GitHub</span>
              </a>
            </Button>
            <Button
              variant="ghost"
              size="icon"
              className="md:hidden"
              onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
            >
              {mobileMenuOpen ? (
                <X className="h-5 w-5" />
              ) : (
                <Menu className="h-5 w-5" />
              )}
            </Button>
          </div>
        </nav>
      </div>
      {mobileMenuOpen && (
        <div className="md:hidden border-t">
          <div className="container mx-auto max-w-7xl px-4">
            <div className="space-y-1 py-4">
              {navigation.map((item) => (
                <Link
                  key={item.href}
                  href={item.href}
                  className={cn(
                    "block px-3 py-2 rounded-md text-base font-medium",
                    pathname?.startsWith(item.href)
                      ? "bg-primary/10 text-primary"
                      : "text-foreground/60 hover:bg-accent"
                  )}
                  onClick={() => setMobileMenuOpen(false)}
                >
                  {item.name}
                </Link>
              ))}
            </div>
          </div>
        </div>
      )}
    </header>
  )
}
