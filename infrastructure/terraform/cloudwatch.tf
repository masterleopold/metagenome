# MinION Metagenomics Pipeline - CloudWatch Configuration
# Monitoring, logging, dashboards and alarms

# ===== CloudWatch Log Groups =====
resource "aws_cloudwatch_log_group" "lambda_logs" {
  for_each = {
    orchestrator     = "/aws/lambda/${var.project_name}-orchestrator"
    phase1_trigger   = "/aws/lambda/${var.project_name}-phase1-trigger"
    phase2_trigger   = "/aws/lambda/${var.project_name}-phase2-trigger"
    phase3_trigger   = "/aws/lambda/${var.project_name}-phase3-trigger"
    phase4_trigger   = "/aws/lambda/${var.project_name}-phase4-trigger"
    phase5_trigger   = "/aws/lambda/${var.project_name}-phase5-trigger"
    phase6_trigger   = "/aws/lambda/${var.project_name}-phase6-trigger"
    validate_input   = "/aws/lambda/${var.project_name}-validate-input"
    check_qc         = "/aws/lambda/${var.project_name}-check-qc"
    update_metadata  = "/aws/lambda/${var.project_name}-update-metadata"
    notify           = "/aws/lambda/${var.project_name}-notify"
  }

  name              = each.value
  retention_in_days = var.environment == "prod" ? 90 : 30

  tags = {
    Name        = each.value
    Environment = var.environment
    Function    = each.key
  }
}

resource "aws_cloudwatch_log_group" "ec2_logs" {
  for_each = {
    basecalling    = "/aws/ec2/${var.project_name}/basecalling"
    qc             = "/aws/ec2/${var.project_name}/qc"
    host_removal   = "/aws/ec2/${var.project_name}/host-removal"
    kraken2        = "/aws/ec2/${var.project_name}/kraken2"
    blast          = "/aws/ec2/${var.project_name}/blast"
    assembly       = "/aws/ec2/${var.project_name}/assembly"
    perv           = "/aws/ec2/${var.project_name}/perv"
    quantification = "/aws/ec2/${var.project_name}/quantification"
    reports        = "/aws/ec2/${var.project_name}/reports"
  }

  name              = each.value
  retention_in_days = var.environment == "prod" ? 90 : 30

  tags = {
    Name        = each.value
    Environment = var.environment
    Phase       = each.key
  }
}

resource "aws_cloudwatch_log_group" "application_logs" {
  name              = "/aws/application/${var.project_name}/main"
  retention_in_days = var.environment == "prod" ? 365 : 90

  tags = {
    Name        = "${var.project_name}-application-logs"
    Environment = var.environment
  }
}

# ===== CloudWatch Log Metric Filters =====
resource "aws_cloudwatch_log_metric_filter" "error_count" {
  name           = "${var.project_name}-error-count"
  pattern        = "[time, request_id, level = ERROR*, msg, ...]"
  log_group_name = aws_cloudwatch_log_group.application_logs.name

  metric_transformation {
    name      = "ErrorCount"
    namespace = "${var.project_name}/Pipeline"
    value     = "1"
    unit      = "Count"
  }
}

resource "aws_cloudwatch_log_metric_filter" "perv_detected" {
  name           = "${var.project_name}-perv-detected"
  pattern        = "[time, ..., msg = *PERV*DETECTED*]"
  log_group_name = aws_cloudwatch_log_group.ec2_logs["perv"].name

  metric_transformation {
    name      = "PERVDetected"
    namespace = "${var.project_name}/Pathogens"
    value     = "1"
    unit      = "Count"
  }
}

resource "aws_cloudwatch_log_metric_filter" "critical_pathogen" {
  name           = "${var.project_name}-critical-pathogen"
  pattern        = "[time, ..., msg = *BSL-3* || msg = *BSL-4*]"
  log_group_name = aws_cloudwatch_log_group.ec2_logs["kraken2"].name

  metric_transformation {
    name      = "CriticalPathogenDetected"
    namespace = "${var.project_name}/Pathogens"
    value     = "1"
    unit      = "Count"
  }
}

resource "aws_cloudwatch_log_metric_filter" "qc_failure" {
  name           = "${var.project_name}-qc-failure"
  pattern        = "[time, ..., msg = *QC*FAILED*]"
  log_group_name = aws_cloudwatch_log_group.ec2_logs["qc"].name

  metric_transformation {
    name      = "QCFailure"
    namespace = "${var.project_name}/Quality"
    value     = "1"
    unit      = "Count"
  }
}

# ===== Custom Metrics Namespace =====
locals {
  custom_namespaces = [
    "${var.project_name}/Pipeline",
    "${var.project_name}/Pathogens",
    "${var.project_name}/Quality",
    "${var.project_name}/Performance",
    "${var.project_name}/Cost"
  ]
}

# ===== CloudWatch Dashboard =====
resource "aws_cloudwatch_dashboard" "main" {
  dashboard_name = "${var.project_name}-main"

  dashboard_body = jsonencode({
    widgets = [
      # Pipeline Status Overview
      {
        type = "metric"
        properties = {
          title   = "Pipeline Status Overview"
          metrics = [
            ["${var.project_name}/Pipeline", "RunsInitiated", { stat = "Sum" }],
            [".", "RunsCompleted", { stat = "Sum" }],
            [".", "RunsFailed", { stat = "Sum" }],
            [".", "RunsInProgress", { stat = "Sum" }]
          ]
          view    = "timeSeries"
          stacked = false
          region  = var.aws_region
          period  = 300
        }
      },
      # Pathogen Detection Summary
      {
        type = "metric"
        properties = {
          title   = "Pathogen Detection Summary"
          metrics = [
            ["${var.project_name}/Pathogens", "TotalPathogensDetected", { stat = "Average" }],
            [".", "PERVDetected", { stat = "Sum" }],
            [".", "CriticalPathogenDetected", { stat = "Sum" }]
          ]
          view   = "singleValue"
          region = var.aws_region
          period = 3600
        }
      },
      # Quality Metrics
      {
        type = "metric"
        properties = {
          title   = "Quality Metrics"
          metrics = [
            ["${var.project_name}/Quality", "MeanQScore", { stat = "Average" }],
            [".", "MedianReadLength", { stat = "Average" }],
            [".", "TotalReads", { stat = "Average" }],
            [".", "HostRemovalRate", { stat = "Average" }]
          ]
          view   = "numberChart"
          region = var.aws_region
          period = 3600
        }
      },
      # EC2 Instance Metrics
      {
        type = "metric"
        properties = {
          title   = "EC2 Instance Utilization"
          metrics = [
            ["AWS/EC2", "CPUUtilization", { stat = "Average", label = "CPU %" }],
            [".", "NetworkIn", { stat = "Sum", label = "Network In" }],
            [".", "NetworkOut", { stat = "Sum", label = "Network Out" }],
            [".", "DiskReadBytes", { stat = "Sum", label = "Disk Read" }],
            [".", "DiskWriteBytes", { stat = "Sum", label = "Disk Write" }]
          ]
          view   = "timeSeries"
          region = var.aws_region
          period = 300
        }
      },
      # Lambda Function Performance
      {
        type = "metric"
        properties = {
          title   = "Lambda Function Performance"
          metrics = [
            ["AWS/Lambda", "Invocations", { stat = "Sum" }],
            [".", "Errors", { stat = "Sum" }],
            [".", "Duration", { stat = "Average" }],
            [".", "ConcurrentExecutions", { stat = "Maximum" }]
          ]
          view   = "timeSeries"
          region = var.aws_region
          period = 300
        }
      },
      # RDS Database Metrics
      {
        type = "metric"
        properties = {
          title   = "RDS Database Performance"
          metrics = [
            ["AWS/RDS", "CPUUtilization", { "dimensions": { "DBInstanceIdentifier": aws_db_instance.metadata.id }, stat = "Average" }],
            [".", "DatabaseConnections", { "dimensions": { "DBInstanceIdentifier": aws_db_instance.metadata.id }, stat = "Average" }],
            [".", "FreeableMemory", { "dimensions": { "DBInstanceIdentifier": aws_db_instance.metadata.id }, stat = "Average" }],
            [".", "ReadLatency", { "dimensions": { "DBInstanceIdentifier": aws_db_instance.metadata.id }, stat = "Average" }],
            [".", "WriteLatency", { "dimensions": { "DBInstanceIdentifier": aws_db_instance.metadata.id }, stat = "Average" }]
          ]
          view   = "timeSeries"
          region = var.aws_region
          period = 300
        }
      },
      # S3 Storage Metrics
      {
        type = "metric"
        properties = {
          title   = "S3 Storage Usage"
          metrics = [
            ["AWS/S3", "BucketSizeBytes", { "dimensions": { "BucketName": aws_s3_bucket.raw_data.id, "StorageType": "StandardStorage" }, stat = "Average", label = "Raw Data" }],
            [".", "BucketSizeBytes", { "dimensions": { "BucketName": aws_s3_bucket.analysis.id, "StorageType": "StandardStorage" }, stat = "Average", label = "Analysis" }],
            [".", "NumberOfObjects", { "dimensions": { "BucketName": aws_s3_bucket.raw_data.id, "StorageType": "AllStorageTypes" }, stat = "Average", label = "Raw Objects" }],
            [".", "NumberOfObjects", { "dimensions": { "BucketName": aws_s3_bucket.analysis.id, "StorageType": "AllStorageTypes" }, stat = "Average", label = "Analysis Objects" }]
          ]
          view   = "timeSeries"
          region = var.aws_region
          period = 86400
        }
      },
      # Error Logs
      {
        type = "log"
        properties = {
          title   = "Recent Errors"
          query   = <<-EOT
            SOURCE '${aws_cloudwatch_log_group.application_logs.name}'
            | fields @timestamp, @message
            | filter level = "ERROR"
            | sort @timestamp desc
            | limit 20
          EOT
          region  = var.aws_region
        }
      },
      # Cost Tracking
      {
        type = "metric"
        properties = {
          title   = "Estimated Daily Cost"
          metrics = [
            ["${var.project_name}/Cost", "EC2Cost", { stat = "Sum" }],
            [".", "S3Cost", { stat = "Sum" }],
            [".", "RDSCost", { stat = "Sum" }],
            [".", "LambdaCost", { stat = "Sum" }],
            [".", "DataTransferCost", { stat = "Sum" }]
          ]
          view   = "singleValue"
          region = var.aws_region
          period = 86400
        }
      }
    ]
  })
}

# ===== CloudWatch Alarms =====

# Critical pathogen detection alarm
resource "aws_cloudwatch_metric_alarm" "critical_pathogen_detected" {
  alarm_name          = "${var.project_name}-critical-pathogen-detected"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "1"
  metric_name         = "CriticalPathogenDetected"
  namespace           = "${var.project_name}/Pathogens"
  period              = "60"
  statistic           = "Sum"
  threshold           = "0"
  alarm_description   = "Critical pathogen (BSL-3/4) detected"
  alarm_actions       = [aws_sns_topic.critical_alerts.arn]

  tags = {
    Name        = "${var.project_name}-critical-pathogen-alarm"
    Environment = var.environment
    Priority    = "CRITICAL"
  }
}

# PERV detection alarm
resource "aws_cloudwatch_metric_alarm" "perv_detected" {
  alarm_name          = "${var.project_name}-perv-detected"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "1"
  metric_name         = "PERVDetected"
  namespace           = "${var.project_name}/Pathogens"
  period              = "60"
  statistic           = "Sum"
  threshold           = "0"
  alarm_description   = "PERV (Porcine Endogenous Retrovirus) detected"
  alarm_actions       = [aws_sns_topic.critical_alerts.arn]

  tags = {
    Name        = "${var.project_name}-perv-alarm"
    Environment = var.environment
    Priority    = "HIGH"
  }
}

# QC failure alarm
resource "aws_cloudwatch_metric_alarm" "qc_failure" {
  alarm_name          = "${var.project_name}-qc-failure"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "1"
  metric_name         = "QCFailure"
  namespace           = "${var.project_name}/Quality"
  period              = "300"
  statistic           = "Sum"
  threshold           = "0"
  alarm_description   = "Quality control check failed"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  tags = {
    Name        = "${var.project_name}-qc-alarm"
    Environment = var.environment
    Priority    = "MEDIUM"
  }
}

# Pipeline failure rate alarm
resource "aws_cloudwatch_metric_alarm" "high_failure_rate" {
  alarm_name          = "${var.project_name}-high-failure-rate"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  threshold           = "20"
  alarm_description   = "High pipeline failure rate detected"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  metric_query {
    id          = "failure_rate"
    expression  = "(m2/m1)*100"
    label       = "Failure Rate"
    return_data = true
  }

  metric_query {
    id = "m1"
    metric {
      metric_name = "RunsInitiated"
      namespace   = "${var.project_name}/Pipeline"
      period      = "3600"
      stat        = "Sum"
    }
  }

  metric_query {
    id = "m2"
    metric {
      metric_name = "RunsFailed"
      namespace   = "${var.project_name}/Pipeline"
      period      = "3600"
      stat        = "Sum"
    }
  }

  tags = {
    Name        = "${var.project_name}-failure-rate-alarm"
    Environment = var.environment
    Priority    = "HIGH"
  }
}

# Lambda error alarm
resource "aws_cloudwatch_metric_alarm" "lambda_errors" {
  for_each = toset([
    "orchestrator",
    "phase1-trigger",
    "phase2-trigger",
    "phase3-trigger",
    "phase4-trigger",
    "phase5-trigger",
    "phase6-trigger"
  ])

  alarm_name          = "${var.project_name}-lambda-${each.key}-errors"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "Errors"
  namespace           = "AWS/Lambda"
  period              = "300"
  statistic           = "Sum"
  threshold           = "5"
  alarm_description   = "Lambda function ${each.key} errors"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    FunctionName = "${var.project_name}-${each.key}"
  }

  tags = {
    Name        = "${var.project_name}-lambda-${each.key}-error-alarm"
    Environment = var.environment
  }
}

# EC2 high CPU alarm
resource "aws_cloudwatch_metric_alarm" "ec2_high_cpu" {
  alarm_name          = "${var.project_name}-ec2-high-cpu"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "3"
  metric_name         = "CPUUtilization"
  namespace           = "AWS/EC2"
  period              = "300"
  statistic           = "Average"
  threshold           = "90"
  alarm_description   = "EC2 instance CPU utilization is too high"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  tags = {
    Name        = "${var.project_name}-ec2-cpu-alarm"
    Environment = var.environment
  }
}

# RDS storage alarm (already defined in rds.tf, but adding metric filter here)
resource "aws_cloudwatch_metric_alarm" "rds_connection_count" {
  alarm_name          = "${var.project_name}-rds-high-connections"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "DatabaseConnections"
  namespace           = "AWS/RDS"
  period              = "300"
  statistic           = "Average"
  threshold           = "80"
  alarm_description   = "RDS connection count is high"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    DBInstanceIdentifier = aws_db_instance.metadata.id
  }

  tags = {
    Name        = "${var.project_name}-rds-connection-alarm"
    Environment = var.environment
  }
}

# ===== CloudWatch Events Rules (for scheduled tasks) =====
resource "aws_cloudwatch_event_rule" "daily_database_backup" {
  name                = "${var.project_name}-daily-db-backup"
  description         = "Trigger daily database backup"
  schedule_expression = "cron(0 2 * * ? *)"  # 2 AM UTC daily

  tags = {
    Name        = "${var.project_name}-db-backup-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_rule" "weekly_reference_db_update" {
  name                = "${var.project_name}-weekly-ref-db-update"
  description         = "Trigger weekly reference database update check"
  schedule_expression = "cron(0 3 ? * SUN *)"  # 3 AM UTC every Sunday

  tags = {
    Name        = "${var.project_name}-ref-db-update-rule"
    Environment = var.environment
  }
}

# ===== CloudWatch Logs Insights Queries (saved) =====
resource "aws_cloudwatch_query_definition" "pathogen_detection_summary" {
  name = "${var.project_name}-pathogen-detection-summary"

  log_group_names = [
    aws_cloudwatch_log_group.ec2_logs["kraken2"].name,
    aws_cloudwatch_log_group.ec2_logs["blast"].name
  ]

  query_string = <<-EOT
    fields @timestamp, pathogen_name, confidence_score, read_count
    | filter @message like /PATHOGEN_DETECTED/
    | stats count(*) as detection_count by pathogen_name
    | sort detection_count desc
  EOT
}

resource "aws_cloudwatch_query_definition" "qc_metrics_summary" {
  name = "${var.project_name}-qc-metrics-summary"

  log_group_names = [
    aws_cloudwatch_log_group.ec2_logs["qc"].name
  ]

  query_string = <<-EOT
    fields @timestamp, run_id, total_reads, mean_qscore, n50
    | filter @message like /QC_METRICS/
    | stats avg(mean_qscore) as avg_quality,
            min(mean_qscore) as min_quality,
            max(mean_qscore) as max_quality,
            avg(total_reads) as avg_reads
  EOT
}

resource "aws_cloudwatch_query_definition" "pipeline_performance" {
  name = "${var.project_name}-pipeline-performance"

  log_group_names = [
    aws_cloudwatch_log_group.lambda_logs["orchestrator"].name
  ]

  query_string = <<-EOT
    fields @timestamp, run_id, phase, duration_seconds
    | filter @message like /PHASE_COMPLETED/
    | stats avg(duration_seconds) as avg_duration,
            max(duration_seconds) as max_duration,
            min(duration_seconds) as min_duration
    by phase
  EOT
}

# ===== CloudWatch Composite Alarms =====
resource "aws_cloudwatch_composite_alarm" "pipeline_critical_failure" {
  alarm_name          = "${var.project_name}-pipeline-critical-failure"
  alarm_description   = "Multiple critical failures in pipeline"
  actions_enabled     = true
  alarm_actions       = [aws_sns_topic.critical_alerts.arn]

  alarm_rule = join(" OR ", [
    "ALARM(${aws_cloudwatch_metric_alarm.critical_pathogen_detected.alarm_name})",
    "ALARM(${aws_cloudwatch_metric_alarm.perv_detected.alarm_name})",
    "(ALARM(${aws_cloudwatch_metric_alarm.high_failure_rate.alarm_name}) AND ALARM(${aws_cloudwatch_metric_alarm.qc_failure.alarm_name}))"
  ])

  tags = {
    Name        = "${var.project_name}-critical-composite-alarm"
    Environment = var.environment
    Priority    = "CRITICAL"
  }
}

# ===== Outputs =====
output "cloudwatch_dashboard_url" {
  description = "URL to the CloudWatch dashboard"
  value       = "https://console.aws.amazon.com/cloudwatch/home?region=${var.aws_region}#dashboards:name=${aws_cloudwatch_dashboard.main.dashboard_name}"
}

output "log_group_names" {
  description = "Map of log group names"
  value = {
    lambda      = { for k, v in aws_cloudwatch_log_group.lambda_logs : k => v.name }
    ec2         = { for k, v in aws_cloudwatch_log_group.ec2_logs : k => v.name }
    application = aws_cloudwatch_log_group.application_logs.name
  }
}

output "alarm_arns" {
  description = "ARNs of CloudWatch alarms"
  value = {
    critical_pathogen = aws_cloudwatch_metric_alarm.critical_pathogen_detected.arn
    perv_detected     = aws_cloudwatch_metric_alarm.perv_detected.arn
    qc_failure        = aws_cloudwatch_metric_alarm.qc_failure.arn
    high_failure_rate = aws_cloudwatch_metric_alarm.high_failure_rate.arn
  }
}