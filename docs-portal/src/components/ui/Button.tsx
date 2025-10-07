import * as React from "react"
import { cn } from "@/lib/utils"

export interface ButtonProps
  extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: "default" | "destructive" | "outline" | "secondary" | "ghost" | "link"
  size?: "default" | "sm" | "lg" | "icon"
  asChild?: boolean
}

const Button = React.forwardRef<HTMLButtonElement, ButtonProps>(
  ({ className, variant = "default", size = "default", asChild, ...props }, ref) => {
    return (
      <button
        className={cn(
          "inline-flex items-center justify-center whitespace-nowrap rounded-md text-sm font-medium transition-all duration-150 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:pointer-events-none disabled:opacity-50",
          {
            // Variants - more subtle, Linear-style
            "bg-foreground text-background hover:opacity-90":
              variant === "default",
            "bg-destructive text-destructive-foreground hover:opacity-90":
              variant === "destructive",
            "border border-border bg-background hover:bg-muted":
              variant === "outline",
            "bg-muted text-foreground hover:bg-muted/80":
              variant === "secondary",
            "hover:bg-muted": variant === "ghost",
            "text-foreground underline-offset-4 hover:underline": variant === "link",
            // Sizes
            "h-9 px-4 py-2 text-sm": size === "default",
            "h-8 px-3 text-xs": size === "sm",
            "h-10 px-6 text-base": size === "lg",
            "h-9 w-9": size === "icon",
          },
          className
        )}
        ref={ref}
        {...props}
      />
    )
  }
)
Button.displayName = "Button"

export { Button }
