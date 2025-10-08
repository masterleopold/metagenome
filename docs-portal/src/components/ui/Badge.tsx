import * as React from "react"
import { cn } from "@/lib/utils"

export interface BadgeProps extends React.HTMLAttributes<HTMLDivElement> {
  variant?: "default" | "secondary" | "success" | "warning" | "error" | "outline"
}

function Badge({ className, variant = "default", ...props }: BadgeProps) {
  return (
    <div
      className={cn(
        "inline-flex items-center rounded-full border px-2.5 py-0.5 text-xs font-semibold transition-colors focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2",
        {
          "border-transparent bg-primary/20 text-primary":
            variant === "default",
          "border-transparent bg-secondary/20 text-secondary":
            variant === "secondary",
          "border-transparent bg-primary/10 text-primary dark:bg-primary/20":
            variant === "success",
          "border-transparent bg-amber-100 text-amber-900 dark:bg-amber-900/30 dark:text-amber-300":
            variant === "warning",
          "border-transparent bg-red-100 text-red-900 dark:bg-red-900/30 dark:text-red-300":
            variant === "error",
          "text-foreground": variant === "outline",
        },
        className
      )}
      {...props}
    />
  )
}

export { Badge }
