import * as React from "react"
import { cn } from "@/lib/utils"

const Alert = React.forwardRef<
  HTMLDivElement,
  React.HTMLAttributes<HTMLDivElement> & {
    variant?: "default" | "destructive" | "warning" | "success"
  }
>(({ className, variant = "default", ...props }, ref) => (
  <div
    ref={ref}
    role="alert"
    className={cn(
      "relative w-full rounded-lg border p-4",
      {
        "bg-background text-foreground": variant === "default",
        "border-red-300 bg-red-50 text-red-950 dark:bg-destructive/10 dark:border-destructive dark:text-red-200 [&>svg]:text-red-800 dark:[&>svg]:text-red-300":
          variant === "destructive",
        "border-amber-300 bg-amber-50 text-amber-950 dark:bg-amber-900/10 dark:border-amber-800 dark:text-amber-200 [&>svg]:text-amber-800 dark:[&>svg]:text-amber-300":
          variant === "warning",
        "border-primary/30 bg-primary/5 text-primary dark:bg-primary/10 dark:text-primary [&>svg]:text-primary":
          variant === "success",
      },
      className
    )}
    {...props}
  />
))
Alert.displayName = "Alert"

const AlertTitle = React.forwardRef<
  HTMLParagraphElement,
  React.HTMLAttributes<HTMLHeadingElement>
>(({ className, ...props }, ref) => (
  <h5
    ref={ref}
    className={cn("my-4 font-medium leading-none tracking-tight", className)}
    {...props}
  />
))
AlertTitle.displayName = "AlertTitle"

const AlertDescription = React.forwardRef<
  HTMLParagraphElement,
  React.HTMLAttributes<HTMLParagraphElement>
>(({ className, ...props }, ref) => (
  <div
    ref={ref}
    className={cn("text-sm [&_p]:leading-relaxed", className)}
    {...props}
  />
))
AlertDescription.displayName = "AlertDescription"

export { Alert, AlertTitle, AlertDescription }
