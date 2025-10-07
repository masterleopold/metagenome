# MinION Metagenomics Pipeline - EventBridge Configuration
# Event-driven orchestration for pipeline phases

# ===== EventBridge Event Bus =====
resource "aws_cloudwatch_event_bus" "pipeline" {
  name = "${var.project_name}-pipeline-bus"

  tags = {
    Name        = "${var.project_name}-pipeline-bus"
    Environment = var.environment
    Purpose     = "Pipeline event orchestration"
  }
}

# ===== S3 Event Rules (Phase 1: Basecalling trigger) =====
resource "aws_cloudwatch_event_rule" "s3_raw_data_uploaded" {
  name           = "${var.project_name}-raw-data-uploaded"
  description    = "Trigger when new FAST5 files are uploaded to S3"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["aws.s3"]
    detail-type = ["Object Created"]
    detail = {
      bucket = {
        name = [aws_s3_bucket.raw_data.id]
      }
      object = {
        key = [
          { suffix = ".fast5" },
          { suffix = ".pod5" }
        ]
      }
    }
  })

  tags = {
    Name        = "${var.project_name}-s3-upload-rule"
    Environment = var.environment
    Phase       = "1"
  }
}

resource "aws_cloudwatch_event_target" "s3_to_orchestrator" {
  rule           = aws_cloudwatch_event_rule.s3_raw_data_uploaded.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "OrchestratorLambda"
  arn            = aws_lambda_function.orchestrator.arn

  input_transformer {
    input_paths = {
      bucket = "$.detail.bucket.name"
      key    = "$.detail.object.key"
      size   = "$.detail.object.size"
    }

    input_template = <<-EOT
      {
        "action": "start_pipeline",
        "source": "s3_upload",
        "bucket": "<bucket>",
        "key": "<key>",
        "size": <size>,
        "phase": 1
      }
    EOT
  }
}

# ===== Phase Completion Events =====
# Rule for phase completion events
resource "aws_cloudwatch_event_rule" "phase_completed" {
  name           = "${var.project_name}-phase-completed"
  description    = "Handle phase completion events"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.pipeline"]
    detail-type = ["Phase Completed"]
    detail = {
      status = ["SUCCESS"]
    }
  })

  tags = {
    Name        = "${var.project_name}-phase-completed-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "phase_to_orchestrator" {
  rule           = aws_cloudwatch_event_rule.phase_completed.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "PhaseOrchestratorLambda"
  arn            = aws_lambda_function.orchestrator.arn

  input_transformer {
    input_paths = {
      run_id        = "$.detail.run_id"
      phase         = "$.detail.phase"
      next_phase    = "$.detail.next_phase"
      metrics       = "$.detail.metrics"
    }

    input_template = <<-EOT
      {
        "action": "trigger_next_phase",
        "run_id": "<run_id>",
        "completed_phase": <phase>,
        "next_phase": <next_phase>,
        "metrics": <metrics>
      }
    EOT
  }
}

# ===== Phase Failure Events =====
resource "aws_cloudwatch_event_rule" "phase_failed" {
  name           = "${var.project_name}-phase-failed"
  description    = "Handle phase failure events"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.pipeline"]
    detail-type = ["Phase Failed"]
    detail = {
      status = ["FAILED"]
    }
  })

  tags = {
    Name        = "${var.project_name}-phase-failed-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "failure_to_sns" {
  rule           = aws_cloudwatch_event_rule.phase_failed.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "FailureSNS"
  arn            = aws_sns_topic.alerts.arn

  input_transformer {
    input_paths = {
      run_id     = "$.detail.run_id"
      phase      = "$.detail.phase"
      error      = "$.detail.error_message"
      timestamp  = "$.time"
    }

    input_template = <<-EOT
      "Pipeline Failure Alert:
      Run ID: <run_id>
      Phase: <phase>
      Error: <error>
      Time: <timestamp>

      Please check CloudWatch logs for details."
    EOT
  }
}

# ===== Critical Pathogen Detection Events =====
resource "aws_cloudwatch_event_rule" "critical_pathogen_detected" {
  name           = "${var.project_name}-critical-pathogen"
  description    = "Alert on critical pathogen detection"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.pathogens"]
    detail-type = ["Critical Pathogen Detected"]
    detail = {
      risk_level = ["BSL-3", "BSL-4"]
    }
  })

  tags = {
    Name        = "${var.project_name}-critical-pathogen-rule"
    Environment = var.environment
    Priority    = "CRITICAL"
  }
}

resource "aws_cloudwatch_event_target" "critical_to_sns" {
  rule           = aws_cloudwatch_event_rule.critical_pathogen_detected.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "CriticalSNS"
  arn            = aws_sns_topic.critical_alerts.arn

  input_transformer {
    input_paths = {
      run_id        = "$.detail.run_id"
      pathogen      = "$.detail.pathogen_name"
      risk_level    = "$.detail.risk_level"
      confidence    = "$.detail.confidence_score"
    }

    input_template = <<-EOT
      "CRITICAL ALERT - Pathogen Detected:
      Run ID: <run_id>
      Pathogen: <pathogen>
      Risk Level: <risk_level>
      Confidence: <confidence>

      IMMEDIATE ACTION REQUIRED"
    EOT
  }
}

# ===== PERV Detection Events =====
resource "aws_cloudwatch_event_rule" "perv_detected" {
  name           = "${var.project_name}-perv-detected"
  description    = "Alert on PERV detection"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.perv"]
    detail-type = ["PERV Detected"]
  })

  tags = {
    Name        = "${var.project_name}-perv-detected-rule"
    Environment = var.environment
    Priority    = "HIGH"
  }
}

resource "aws_cloudwatch_event_target" "perv_to_sns" {
  rule           = aws_cloudwatch_event_rule.perv_detected.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "PERVSNS"
  arn            = aws_sns_topic.critical_alerts.arn

  input_transformer {
    input_paths = {
      run_id     = "$.detail.run_id"
      perv_type  = "$.detail.perv_type"
      read_count = "$.detail.read_count"
      full_length = "$.detail.full_length_detected"
    }

    input_template = <<-EOT
      "PERV Detection Alert:
      Run ID: <run_id>
      PERV Type: <perv_type>
      Read Count: <read_count>
      Full-length Detected: <full_length>

      Review required for xenotransplantation safety"
    EOT
  }
}

# ===== QC Failure Events =====
resource "aws_cloudwatch_event_rule" "qc_failed" {
  name           = "${var.project_name}-qc-failed"
  description    = "Handle QC failure events"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.qc"]
    detail-type = ["QC Failed"]
  })

  tags = {
    Name        = "${var.project_name}-qc-failed-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "qc_to_lambda" {
  rule           = aws_cloudwatch_event_rule.qc_failed.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "QCFailureLambda"
  arn            = aws_lambda_function.check_qc.arn
}

# ===== Pipeline Complete Events =====
resource "aws_cloudwatch_event_rule" "pipeline_complete" {
  name           = "${var.project_name}-pipeline-complete"
  description    = "Pipeline completion notification"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source      = ["${var.project_name}.pipeline"]
    detail-type = ["Pipeline Complete"]
  })

  tags = {
    Name        = "${var.project_name}-pipeline-complete-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "complete_to_sns" {
  rule           = aws_cloudwatch_event_rule.pipeline_complete.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "CompleteSNS"
  arn            = aws_sns_topic.pipeline_status.arn

  input_transformer {
    input_paths = {
      run_id            = "$.detail.run_id"
      sample_id         = "$.detail.sample_id"
      total_duration    = "$.detail.total_duration_hours"
      pathogens_found   = "$.detail.pathogens_detected"
      report_url        = "$.detail.report_url"
    }

    input_template = <<-EOT
      "Pipeline Complete:
      Run ID: <run_id>
      Sample: <sample_id>
      Duration: <total_duration> hours
      Pathogens Detected: <pathogens_found>
      Report: <report_url>"
    EOT
  }
}

# ===== Scheduled Events =====
# Daily system check
resource "aws_cloudwatch_event_rule" "daily_system_check" {
  name                = "${var.project_name}-daily-system-check"
  description         = "Daily system health check"
  schedule_expression = "cron(0 6 * * ? *)"  # 6 AM UTC daily

  tags = {
    Name        = "${var.project_name}-daily-check-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "daily_check_lambda" {
  rule      = aws_cloudwatch_event_rule.daily_system_check.name
  target_id = "DailyCheckLambda"
  arn       = aws_lambda_function.orchestrator.arn

  input = jsonencode({
    action = "system_check"
    type   = "daily"
  })
}

# Weekly database maintenance
resource "aws_cloudwatch_event_rule" "weekly_maintenance" {
  name                = "${var.project_name}-weekly-maintenance"
  description         = "Weekly database maintenance"
  schedule_expression = "cron(0 2 ? * SUN *)"  # 2 AM UTC every Sunday

  tags = {
    Name        = "${var.project_name}-weekly-maintenance-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "maintenance_lambda" {
  rule      = aws_cloudwatch_event_rule.weekly_maintenance.name
  target_id = "MaintenanceLambda"
  arn       = aws_lambda_function.update_metadata.arn

  input = jsonencode({
    action = "database_maintenance"
    type   = "weekly"
  })
}

# Monthly reference database update check
resource "aws_cloudwatch_event_rule" "monthly_db_update" {
  name                = "${var.project_name}-monthly-db-update"
  description         = "Monthly reference database update check"
  schedule_expression = "cron(0 3 1 * ? *)"  # 3 AM UTC on the 1st of each month

  tags = {
    Name        = "${var.project_name}-monthly-db-update-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "db_update_lambda" {
  rule      = aws_cloudwatch_event_rule.monthly_db_update.name
  target_id = "DBUpdateLambda"
  arn       = aws_lambda_function.orchestrator.arn

  input = jsonencode({
    action = "check_database_updates"
    databases = ["kraken2", "blast", "rvdb", "pmda"]
  })
}

# ===== EC2 Instance State Change Events =====
resource "aws_cloudwatch_event_rule" "ec2_state_change" {
  name        = "${var.project_name}-ec2-state-change"
  description = "Monitor EC2 instance state changes"

  event_pattern = jsonencode({
    source      = ["aws.ec2"]
    detail-type = ["EC2 Instance State-change Notification"]
    detail = {
      state = ["running", "stopped", "terminated"]
    }
  })

  tags = {
    Name        = "${var.project_name}-ec2-state-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "ec2_state_to_lambda" {
  rule      = aws_cloudwatch_event_rule.ec2_state_change.name
  target_id = "EC2StateLambda"
  arn       = aws_lambda_function.update_metadata.arn

  input_transformer {
    input_paths = {
      instance_id = "$.detail.instance-id"
      state       = "$.detail.state"
      time        = "$.time"
    }

    input_template = <<-EOT
      {
        "action": "update_ec2_state",
        "instance_id": "<instance_id>",
        "state": "<state>",
        "timestamp": "<time>"
      }
    EOT
  }
}

# ===== Lambda Permissions for EventBridge =====
resource "aws_lambda_permission" "eventbridge_invoke_orchestrator" {
  for_each = toset([
    aws_cloudwatch_event_rule.s3_raw_data_uploaded.arn,
    aws_cloudwatch_event_rule.phase_completed.arn,
    aws_cloudwatch_event_rule.daily_system_check.arn,
    aws_cloudwatch_event_rule.monthly_db_update.arn
  ])

  statement_id  = "AllowEventBridgeInvoke-${sha256(each.key)}"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.orchestrator.function_name
  principal     = "events.amazonaws.com"
  source_arn    = each.key
}

resource "aws_lambda_permission" "eventbridge_invoke_check_qc" {
  statement_id  = "AllowEventBridgeInvokeCheckQC"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.check_qc.function_name
  principal     = "events.amazonaws.com"
  source_arn    = aws_cloudwatch_event_rule.qc_failed.arn
}

resource "aws_lambda_permission" "eventbridge_invoke_update_metadata" {
  for_each = toset([
    aws_cloudwatch_event_rule.weekly_maintenance.arn,
    aws_cloudwatch_event_rule.ec2_state_change.arn
  ])

  statement_id  = "AllowEventBridgeInvokeUpdateMetadata-${sha256(each.key)}"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.update_metadata.function_name
  principal     = "events.amazonaws.com"
  source_arn    = each.key
}

# ===== Event Archive =====
resource "aws_cloudwatch_event_archive" "pipeline_events" {
  name             = "${var.project_name}-pipeline-archive"
  description      = "Archive of all pipeline events for audit"
  event_source_arn = aws_cloudwatch_event_bus.pipeline.arn
  retention_days   = var.environment == "prod" ? 1827 : 90  # 5 years for prod (PMDA requirement)

  event_pattern = jsonencode({
    source = [
      "${var.project_name}.pipeline",
      "${var.project_name}.pathogens",
      "${var.project_name}.perv",
      "${var.project_name}.qc"
    ]
  })
}

# ===== Custom Event Patterns for Pipeline Monitoring =====
resource "aws_cloudwatch_event_rule" "pipeline_metrics" {
  name           = "${var.project_name}-pipeline-metrics"
  description    = "Capture pipeline metrics for CloudWatch"
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name

  event_pattern = jsonencode({
    source = ["${var.project_name}.pipeline"]
    detail-type = [
      "Phase Started",
      "Phase Completed",
      "Phase Failed",
      "Pipeline Complete"
    ]
  })

  tags = {
    Name        = "${var.project_name}-pipeline-metrics-rule"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_event_target" "metrics_to_cloudwatch" {
  rule           = aws_cloudwatch_event_rule.pipeline_metrics.name
  event_bus_name = aws_cloudwatch_event_bus.pipeline.name
  target_id      = "CloudWatchMetrics"
  arn            = aws_cloudwatch_log_group.application_logs.arn
  role_arn       = aws_iam_role.eventbridge_logs.arn
}

# IAM Role for EventBridge to write to CloudWatch Logs
resource "aws_iam_role" "eventbridge_logs" {
  name = "${var.project_name}-eventbridge-logs-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "events.amazonaws.com"
        }
      }
    ]
  })

  tags = {
    Name        = "${var.project_name}-eventbridge-logs-role"
    Environment = var.environment
  }
}

resource "aws_iam_role_policy" "eventbridge_logs" {
  name = "${var.project_name}-eventbridge-logs-policy"
  role = aws_iam_role.eventbridge_logs.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "logs:CreateLogStream",
          "logs:PutLogEvents"
        ]
        Resource = "${aws_cloudwatch_log_group.application_logs.arn}:*"
      }
    ]
  })
}

# ===== Outputs =====
output "event_bus_name" {
  description = "Name of the EventBridge event bus"
  value       = aws_cloudwatch_event_bus.pipeline.name
}

output "event_bus_arn" {
  description = "ARN of the EventBridge event bus"
  value       = aws_cloudwatch_event_bus.pipeline.arn
}

output "event_rules" {
  description = "Map of EventBridge rules"
  value = {
    s3_upload         = aws_cloudwatch_event_rule.s3_raw_data_uploaded.name
    phase_completed   = aws_cloudwatch_event_rule.phase_completed.name
    phase_failed      = aws_cloudwatch_event_rule.phase_failed.name
    critical_pathogen = aws_cloudwatch_event_rule.critical_pathogen_detected.name
    perv_detected     = aws_cloudwatch_event_rule.perv_detected.name
    qc_failed         = aws_cloudwatch_event_rule.qc_failed.name
    pipeline_complete = aws_cloudwatch_event_rule.pipeline_complete.name
  }
}

output "scheduled_rules" {
  description = "Map of scheduled EventBridge rules"
  value = {
    daily_check        = aws_cloudwatch_event_rule.daily_system_check.name
    weekly_maintenance = aws_cloudwatch_event_rule.weekly_maintenance.name
    monthly_db_update  = aws_cloudwatch_event_rule.monthly_db_update.name
  }
}

output "event_archive_name" {
  description = "Name of the event archive"
  value       = aws_cloudwatch_event_archive.pipeline_events.name
}