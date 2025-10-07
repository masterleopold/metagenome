'use client'

import * as React from "react"
import { cn } from "@/lib/utils"
import { Check, Copy } from "lucide-react"
import { Button } from "./Button"

interface CodeBlockProps extends React.HTMLAttributes<HTMLDivElement> {
  code: string
  language?: string
  filename?: string
  showLineNumbers?: boolean
}

export function CodeBlock({
  code,
  language = "typescript",
  filename,
  showLineNumbers = false,
  className,
  ...props
}: CodeBlockProps) {
  const [copied, setCopied] = React.useState(false)

  const copyToClipboard = async () => {
    await navigator.clipboard.writeText(code)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }

  const lines = code.split('\n')

  return (
    <div
      className={cn(
        "relative rounded-lg border bg-slate-950 dark:bg-slate-900",
        className
      )}
      {...props}
    >
      {filename && (
        <div className="flex items-center justify-between border-b border-slate-800 px-4 py-2">
          <span className="text-sm font-mono text-slate-400">{filename}</span>
          <span className="text-xs text-slate-500 uppercase">{language}</span>
        </div>
      )}
      <div className="relative">
        <Button
          size="icon"
          variant="ghost"
          className="absolute right-2 top-2 h-8 w-8 text-slate-400 hover:text-slate-100 hover:bg-slate-800"
          onClick={copyToClipboard}
        >
          {copied ? (
            <Check className="h-4 w-4" />
          ) : (
            <Copy className="h-4 w-4" />
          )}
        </Button>
        <pre className="overflow-x-auto p-4 text-sm">
          <code className="font-mono text-slate-100">
            {showLineNumbers ? (
              lines.map((line, i) => (
                <div key={i} className="flex">
                  <span className="mr-4 w-8 text-right text-slate-600 select-none">
                    {i + 1}
                  </span>
                  <span>{line}</span>
                </div>
              ))
            ) : (
              code
            )}
          </code>
        </pre>
      </div>
    </div>
  )
}

interface InlineCodeProps extends React.HTMLAttributes<HTMLElement> {}

export function InlineCode({ className, ...props }: InlineCodeProps) {
  return (
    <code
      className={cn(
        "relative rounded bg-muted px-[0.3rem] py-[0.2rem] font-mono text-sm font-semibold",
        className
      )}
      {...props}
    />
  )
}
