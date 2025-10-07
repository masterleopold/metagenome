# MinION Metagenomics Pipeline - SNS Configuration
# SNS topics for notifications and alerts

# ===== Main Alert Topic =====
resource "aws_sns_topic" "alerts" {
  name         = "${var.project_name}-alerts"
  display_name = "MinION Pipeline Alerts"

  kms_master_key_id = "alias/aws/sns"  # Enable encryption

  tags = {
    Name        = "${var.project_name}-alerts"
    Environment = var.environment
    Purpose     = "Pipeline alerts and notifications"
  }
}

# ===== Critical Alerts Topic (PERV, BSL-3) =====
resource "aws_sns_topic" "critical_alerts" {
  name         = "${var.project_name}-critical-alerts"
  display_name = "MinION CRITICAL Alerts"

  kms_master_key_id = "alias/aws/sns"

  tags = {
    Name        = "${var.project_name}-critical-alerts"
    Environment = var.environment
    Purpose     = "Critical pathogen detection alerts"
    Priority    = "CRITICAL"
  }
}

# ===== Pipeline Status Topic =====
resource "aws_sns_topic" "pipeline_status" {
  name         = "${var.project_name}-pipeline-status"
  display_name = "MinION Pipeline Status"

  kms_master_key_id = "alias/aws/sns"

  tags = {
    Name        = "${var.project_name}-pipeline-status"
    Environment = var.environment
    Purpose     = "Pipeline completion status"
  }
}

# ===== Email Subscriptions =====
resource "aws_sns_topic_subscription" "alerts_email" {
  count = var.notification_email != "" ? 1 : 0

  topic_arn = aws_sns_topic.alerts.arn
  protocol  = "email"
  endpoint  = var.notification_email

  depends_on = [aws_sns_topic.alerts]
}

resource "aws_sns_topic_subscription" "critical_alerts_email" {
  count = var.notification_email != "" ? 1 : 0

  topic_arn = aws_sns_topic.critical_alerts.arn
  protocol  = "email"
  endpoint  = var.notification_email

  depends_on = [aws_sns_topic.critical_alerts]
}

resource "aws_sns_topic_subscription" "pipeline_status_email" {
  count = var.notification_email != "" ? 1 : 0

  topic_arn = aws_sns_topic.pipeline_status.arn
  protocol  = "email"
  endpoint  = var.notification_email

  depends_on = [aws_sns_topic.pipeline_status]
}

# ===== SMS Subscriptions (Optional) =====
variable "sms_phone_number" {
  description = "Phone number for SMS alerts (optional)"
  type        = string
  default     = ""
  sensitive   = true
}

resource "aws_sns_topic_subscription" "critical_alerts_sms" {
  count = var.sms_phone_number != "" ? 1 : 0

  topic_arn = aws_sns_topic.critical_alerts.arn
  protocol  = "sms"
  endpoint  = var.sms_phone_number

  depends_on = [aws_sns_topic.critical_alerts]
}

# ===== Lambda Subscriptions =====
resource "aws_sns_topic_subscription" "alerts_to_notify_lambda" {
  topic_arn = aws_sns_topic.alerts.arn
  protocol  = "lambda"
  endpoint  = aws_lambda_function.notify.arn

  depends_on = [
    aws_sns_topic.alerts,
    aws_lambda_function.notify
  ]
}

# Lambda permission for SNS to invoke
resource "aws_lambda_permission" "sns_invoke_notify" {
  statement_id  = "AllowSNSInvokeNotifyLambda"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.notify.function_name
  principal     = "sns.amazonaws.com"
  source_arn    = aws_sns_topic.alerts.arn
}

# ===== SNS Topic Policies =====
resource "aws_sns_topic_policy" "alerts" {
  arn = aws_sns_topic.alerts.arn

  policy = jsonencode({
    Version = "2012-10-17"
    Id      = "alerts-topic-policy"
    Statement = [
      {
        Sid    = "AllowCloudWatchPublish"
        Effect = "Allow"
        Principal = {
          Service = "cloudwatch.amazonaws.com"
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.alerts.arn
        Condition = {
          StringEquals = {
            "aws:SourceAccount" = data.aws_caller_identity.current.account_id
          }
        }
      },
      {
        Sid    = "AllowLambdaPublish"
        Effect = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.alerts.arn
        Condition = {
          StringEquals = {
            "aws:SourceAccount" = data.aws_caller_identity.current.account_id
          }
        }
      },
      {
        Sid    = "AllowEC2Publish"
        Effect = "Allow"
        Principal = {
          AWS = aws_iam_role.ec2_instance.arn
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.alerts.arn
      }
    ]
  })
}

resource "aws_sns_topic_policy" "critical_alerts" {
  arn = aws_sns_topic.critical_alerts.arn

  policy = jsonencode({
    Version = "2012-10-17"
    Id      = "critical-alerts-topic-policy"
    Statement = [
      {
        Sid    = "AllowLambdaPublish"
        Effect = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.critical_alerts.arn
        Condition = {
          StringEquals = {
            "aws:SourceAccount" = data.aws_caller_identity.current.account_id
          }
        }
      },
      {
        Sid    = "AllowEC2Publish"
        Effect = "Allow"
        Principal = {
          AWS = aws_iam_role.ec2_instance.arn
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.critical_alerts.arn
      }
    ]
  })
}

resource "aws_sns_topic_policy" "pipeline_status" {
  arn = aws_sns_topic.pipeline_status.arn

  policy = jsonencode({
    Version = "2012-10-17"
    Id      = "pipeline-status-topic-policy"
    Statement = [
      {
        Sid    = "AllowLambdaPublish"
        Effect = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
        Action = [
          "SNS:Publish"
        ]
        Resource = aws_sns_topic.pipeline_status.arn
        Condition = {
          StringEquals = {
            "aws:SourceAccount" = data.aws_caller_identity.current.account_id
          }
        }
      }
    ]
  })
}

# ===== SNS Data Protection Policy (Sensitive Data Filtering) =====
resource "aws_sns_data_protection_policy" "alerts" {
  arn = aws_sns_topic.alerts.arn

  policy = jsonencode({
    Name        = "AlertsDataProtection"
    Description = "Redact sensitive information from alerts"
    Version     = "2021-06-01"
    Statement = [
      {
        DataIdentifier = [
          "arn:aws:dataprotection::aws:data-identifier/EmailAddress",
          "arn:aws:dataprotection::aws:data-identifier/IpAddress"
        ]
        Operation = {
          Deidentify = {
            MaskConfig = {
              MaskWithCharacter = "*"
            }
          }
        }
      }
    ]
  })
}

# ===== CloudWatch Log Group for SNS Failed Deliveries =====
resource "aws_cloudwatch_log_group" "sns_failures" {
  name              = "/aws/sns/${var.project_name}/failures"
  retention_in_days = 7

  tags = {
    Name        = "${var.project_name}-sns-failures"
    Environment = var.environment
  }
}

# ===== SNS Delivery Status Logging =====
resource "aws_sns_topic_subscription" "alerts_failure_logging" {
  topic_arn = aws_sns_topic.alerts.arn
  protocol  = "https"
  endpoint  = "https://logs.${var.aws_region}.amazonaws.com/sns/deliver"

  raw_message_delivery = true

  subscription_role_arn = aws_iam_role.sns_delivery_status.arn

  depends_on = [
    aws_sns_topic.alerts,
    aws_iam_role.sns_delivery_status
  ]
}

# IAM Role for SNS Delivery Status Logging
resource "aws_iam_role" "sns_delivery_status" {
  name = "${var.project_name}-sns-delivery-status-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "sns.amazonaws.com"
        }
      }
    ]
  })

  tags = {
    Name        = "${var.project_name}-sns-delivery-status-role"
    Environment = var.environment
  }
}

resource "aws_iam_role_policy" "sns_delivery_status" {
  name = "${var.project_name}-sns-delivery-status-policy"
  role = aws_iam_role.sns_delivery_status.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "logs:CreateLogGroup",
          "logs:CreateLogStream",
          "logs:PutLogEvents",
          "logs:PutMetricFilter",
          "logs:PutRetentionPolicy"
        ]
        Resource = [
          "arn:aws:logs:${var.aws_region}:${data.aws_caller_identity.current.account_id}:log-group:/aws/sns/*"
        ]
      }
    ]
  })
}

# ===== Outputs =====
output "sns_alerts_topic_arn" {
  description = "ARN of the main alerts SNS topic"
  value       = aws_sns_topic.alerts.arn
}

output "sns_critical_alerts_topic_arn" {
  description = "ARN of the critical alerts SNS topic"
  value       = aws_sns_topic.critical_alerts.arn
}

output "sns_pipeline_status_topic_arn" {
  description = "ARN of the pipeline status SNS topic"
  value       = aws_sns_topic.pipeline_status.arn
}

output "sns_alerts_topic_name" {
  description = "Name of the main alerts SNS topic"
  value       = aws_sns_topic.alerts.name
}

output "sns_critical_alerts_topic_name" {
  description = "Name of the critical alerts SNS topic"
  value       = aws_sns_topic.critical_alerts.name
}

output "sns_pipeline_status_topic_name" {
  description = "Name of the pipeline status SNS topic"
  value       = aws_sns_topic.pipeline_status.name
}