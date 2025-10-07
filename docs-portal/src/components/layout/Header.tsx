'use client'

import Link from "next/link"
import { Github, Moon, Sun } from "lucide-react"
import { Button } from "@/components/ui/Button"
import { useTheme } from "next-themes"

export function Header() {
  const { theme, setTheme } = useTheme()

  return (
    <header className="sticky top-0 z-50 w-full border-b bg-background/80 backdrop-blur-md">
      <nav className="flex h-14 items-center justify-between px-4">
        <Link href="/" className="flex items-center space-x-2">
          <div className="h-7 w-7 rounded-md bg-foreground flex items-center justify-center">
            <span className="text-background font-bold text-base">M</span>
          </div>
          <span className="font-semibold text-sm">
            MinION Pipeline
          </span>
        </Link>
        <div className="flex items-center gap-2">
          <Button
            variant="ghost"
            size="icon"
            onClick={() => setTheme(theme === "dark" ? "light" : "dark")}
          >
            <Sun className="h-5 w-5 rotate-0 scale-100 transition-all dark:-rotate-90 dark:scale-0" />
            <Moon className="absolute h-5 w-5 rotate-90 scale-0 transition-all dark:rotate-0 dark:scale-100" />
            <span className="sr-only">Toggle theme</span>
          </Button>
          <Button variant="ghost" size="icon" asChild>
            <a
              href="https://github.com/masterleopold/metagenome"
              target="_blank"
              rel="noreferrer"
            >
              <Github className="h-5 w-5" />
              <span className="sr-only">GitHub</span>
            </a>
          </Button>
        </div>
      </nav>
    </header>
  )
}
