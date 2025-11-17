# Terraform Configuration for 4-Virus Surveillance System
# Version: 1.0
# Region: ap-northeast-1 (Tokyo)

terraform {
  required_version = ">= 1.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }

  backend "s3" {
    bucket = "minion-terraform-state"
    key    = "surveillance/terraform.tfstate"
    region = "ap-northeast-1"
  }
}

provider "aws" {
  region = var.aws_region

  default_tags {
    tags = {
      Project    = "MinION-Surveillance"
      Component  = "4-Virus-Detection"
      ManagedBy  = "Terraform"
      Environment = var.environment
    }
  }
}

# Variables
variable "aws_region" {
  description = "AWS Region"
  type        = string
  default     = "ap-northeast-1"
}

variable "environment" {
  description = "Environment (dev, staging, prod)"
  type        = string
  default     = "prod"
}

variable "surveillance_bucket_name" {
  description = "S3 bucket for surveillance data"
  type        = string
  default     = "surveillance-data"
}

variable "estat_app_id" {
  description = "E-Stat API Application ID"
  type        = string
  sensitive   = true
}

variable "pubmed_email" {
  description = "Email for PubMed API"
  type        = string
}

# ============================================================================
# S3 Bucket for Surveillance Data
# ============================================================================

resource "aws_s3_bucket" "surveillance_data" {
  bucket = var.surveillance_bucket_name
}

resource "aws_s3_bucket_versioning" "surveillance_data" {
  bucket = aws_s3_bucket.surveillance_data.id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_lifecycle_configuration" "surveillance_data" {
  bucket = aws_s3_bucket.surveillance_data.id

  # J-STAGE ToS compliance: Article 3, Clause 5
  # Academic data (includes J-STAGE) must be deleted after 24 hours
  rule {
    id     = "expire-academic-after-24h"
    status = "Enabled"

    filter {
      prefix = "external/academic/"
    }

    expiration {
      days = 1  # 24 hours (J-STAGE Terms of Service compliance)
    }
  }

  # Other external sources (MAFF, E-Stat) - separate retention policy
  rule {
    id     = "expire-other-external-after-1-year"
    status = "Enabled"

    filter {
      prefix = "external/maff/"
    }

    expiration {
      days = 365
    }
  }

  rule {
    id     = "expire-estat-after-1-year"
    status = "Enabled"

    filter {
      prefix = "external/estat/"
    }

    expiration {
      days = 365
    }
  }

  rule {
    id     = "expire-internal-after-2-years"
    status = "Enabled"

    filter {
      prefix = "internal/"
    }

    expiration {
      days = 730
    }
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "surveillance_data" {
  bucket = aws_s3_bucket.surveillance_data.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

# ============================================================================
# DynamoDB Tables
# ============================================================================

# Detections Table
resource "aws_dynamodb_table" "surveillance_detections" {
  name         = "surveillance-detections"
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "detection_id"
  range_key    = "timestamp"

  attribute {
    name = "detection_id"
    type = "S"
  }

  attribute {
    name = "timestamp"
    type = "S"
  }

  attribute {
    name = "virus_type"
    type = "S"
  }

  attribute {
    name = "severity"
    type = "S"
  }

  attribute {
    name = "timestamp_sort"
    type = "N"
  }

  global_secondary_index {
    name            = "virus-type-index"
    hash_key        = "virus_type"
    range_key       = "timestamp_sort"
    projection_type = "ALL"
  }

  global_secondary_index {
    name            = "severity-index"
    hash_key        = "severity"
    range_key       = "timestamp_sort"
    projection_type = "ALL"
  }

  stream_enabled   = true
  stream_view_type = "NEW_AND_OLD_IMAGES"

  point_in_time_recovery {
    enabled = true
  }

  tags = {
    Name = "surveillance-detections"
  }
}

# External Updates Table
resource "aws_dynamodb_table" "surveillance_external_updates" {
  name         = "surveillance-external-updates"
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "source#date"
  range_key    = "update_id"

  attribute {
    name = "source#date"
    type = "S"
  }

  attribute {
    name = "update_id"
    type = "S"
  }

  attribute {
    name = "source"
    type = "S"
  }

  attribute {
    name = "timestamp"
    type = "N"
  }

  global_secondary_index {
    name            = "source-timestamp-index"
    hash_key        = "source"
    range_key       = "timestamp"
    projection_type = "ALL"
  }

  ttl {
    attribute_name = "ttl"
    enabled        = true
  }

  tags = {
    Name = "surveillance-external-updates"
  }
}

# Notifications Table
resource "aws_dynamodb_table" "surveillance_notifications" {
  name         = "surveillance-notifications"
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "notification_id"
  range_key    = "timestamp"

  attribute {
    name = "notification_id"
    type = "S"
  }

  attribute {
    name = "timestamp"
    type = "S"
  }

  ttl {
    attribute_name = "ttl"
    enabled        = true
  }

  tags = {
    Name = "surveillance-notifications"
  }
}

# ============================================================================
# SNS Topics for Alerts
# ============================================================================

resource "aws_sns_topic" "critical_alerts" {
  name         = "4virus-critical-alerts"
  display_name = "Critical 4-Virus Detections"
}

resource "aws_sns_topic" "high_alerts" {
  name         = "4virus-high-alerts"
  display_name = "High Priority 4-Virus Detections"
}

resource "aws_sns_topic" "daily_summary" {
  name         = "4virus-daily-summary"
  display_name = "Daily Surveillance Summary"
}

# ============================================================================
# Lambda Function: External Collector
# ============================================================================

resource "aws_iam_role" "external_collector_role" {
  name = "surveillance-external-collector-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
      }
    ]
  })
}

resource "aws_iam_role_policy" "external_collector_policy" {
  name = "surveillance-external-collector-policy"
  role = aws_iam_role.external_collector_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:PutObject",
          "s3:GetObject",
          "s3:ListBucket"
        ]
        Resource = [
          aws_s3_bucket.surveillance_data.arn,
          "${aws_s3_bucket.surveillance_data.arn}/*"
        ]
      },
      {
        Effect = "Allow"
        Action = [
          "dynamodb:PutItem",
          "dynamodb:GetItem",
          "dynamodb:Query",
          "dynamodb:Scan"
        ]
        Resource = [
          aws_dynamodb_table.surveillance_detections.arn,
          aws_dynamodb_table.surveillance_external_updates.arn
        ]
      },
      {
        Effect = "Allow"
        Action = [
          "sns:Publish"
        ]
        Resource = [
          aws_sns_topic.critical_alerts.arn,
          aws_sns_topic.high_alerts.arn
        ]
      },
      {
        Effect = "Allow"
        Action = [
          "logs:CreateLogGroup",
          "logs:CreateLogStream",
          "logs:PutLogEvents"
        ]
        Resource = "arn:aws:logs:${var.aws_region}:*:*"
      }
    ]
  })
}

resource "aws_lambda_function" "external_collector" {
  filename      = "lambda_packages/external_collector.zip"
  function_name = "surveillance-external-collector"
  role          = aws_iam_role.external_collector_role.arn
  handler       = "handler.lambda_handler"
  runtime       = "python3.11"
  timeout       = 900 # 15 minutes
  memory_size   = 512

  environment {
    variables = {
      SURVEILLANCE_BUCKET = aws_s3_bucket.surveillance_data.bucket
      E_STAT_APP_ID       = var.estat_app_id
      PUBMED_EMAIL        = var.pubmed_email
      AWS_REGION          = var.aws_region
    }
  }

  depends_on = [
    aws_iam_role_policy.external_collector_policy
  ]
}

# CloudWatch Log Group
resource "aws_cloudwatch_log_group" "external_collector" {
  name              = "/aws/lambda/${aws_lambda_function.external_collector.function_name}"
  retention_in_days = 30
}

# ============================================================================
# EventBridge Schedule for Daily Collection
# ============================================================================

resource "aws_cloudwatch_event_rule" "daily_collection" {
  name                = "surveillance-daily-collection"
  description         = "Trigger daily external source collection at 11:00 JST (02:00 UTC)"
  schedule_expression = "cron(0 2 * * ? *)" # 11:00 JST
}

resource "aws_cloudwatch_event_target" "daily_collection_target" {
  rule      = aws_cloudwatch_event_rule.daily_collection.name
  target_id = "ExternalCollector"
  arn       = aws_lambda_function.external_collector.arn
}

resource "aws_lambda_permission" "allow_eventbridge" {
  statement_id  = "AllowExecutionFromEventBridge"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.external_collector.function_name
  principal     = "events.amazonaws.com"
  source_arn    = aws_cloudwatch_event_rule.daily_collection.arn
}

# ============================================================================
# Outputs
# ============================================================================

output "surveillance_bucket" {
  description = "Surveillance data S3 bucket name"
  value       = aws_s3_bucket.surveillance_data.bucket
}

output "detections_table" {
  description = "Detections DynamoDB table name"
  value       = aws_dynamodb_table.surveillance_detections.name
}

output "external_updates_table" {
  description = "External updates DynamoDB table name"
  value       = aws_dynamodb_table.surveillance_external_updates.name
}

output "critical_alerts_topic" {
  description = "Critical alerts SNS topic ARN"
  value       = aws_sns_topic.critical_alerts.arn
}

output "high_alerts_topic" {
  description = "High alerts SNS topic ARN"
  value       = aws_sns_topic.high_alerts.arn
}

output "external_collector_function" {
  description = "External collector Lambda function name"
  value       = aws_lambda_function.external_collector.function_name
}
