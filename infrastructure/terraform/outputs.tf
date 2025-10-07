# MinION Metagenomics Pipeline - Terraform Outputs
# Consolidated outputs from all modules for reference and integration

# ===== VPC and Networking Outputs =====
output "vpc_id" {
  description = "ID of the VPC"
  value       = module.vpc.vpc_id
}

output "vpc_cidr" {
  description = "CIDR block of the VPC"
  value       = module.vpc.vpc_cidr_block
}

output "public_subnets" {
  description = "List of IDs of public subnets"
  value       = module.vpc.public_subnets
}

output "private_subnets" {
  description = "List of IDs of private subnets"
  value       = module.vpc.private_subnets
}

output "database_subnets" {
  description = "List of IDs of database subnets"
  value       = module.vpc.database_subnets
}

output "nat_gateway_ids" {
  description = "List of NAT Gateway IDs"
  value       = module.vpc.natgw_ids
}

# ===== S3 Bucket Outputs =====
output "s3_buckets" {
  description = "S3 bucket information"
  value = {
    raw_data = {
      id          = aws_s3_bucket.raw_data.id
      arn         = aws_s3_bucket.raw_data.arn
      domain_name = aws_s3_bucket.raw_data.bucket_domain_name
      region      = aws_s3_bucket.raw_data.region
    }
    analysis = {
      id          = aws_s3_bucket.analysis.id
      arn         = aws_s3_bucket.analysis.arn
      domain_name = aws_s3_bucket.analysis.bucket_domain_name
      region      = aws_s3_bucket.analysis.region
    }
    scripts = {
      id          = aws_s3_bucket.scripts.id
      arn         = aws_s3_bucket.scripts.arn
      domain_name = aws_s3_bucket.scripts.bucket_domain_name
      region      = aws_s3_bucket.scripts.region
    }
  }
  sensitive = false
}

# ===== RDS Database Outputs =====
output "rds_database" {
  description = "RDS database connection information"
  value = {
    endpoint        = aws_db_instance.metadata.endpoint
    address         = aws_db_instance.metadata.address
    port            = aws_db_instance.metadata.port
    database_name   = aws_db_instance.metadata.db_name
    instance_id     = aws_db_instance.metadata.id
    arn             = aws_db_instance.metadata.arn
    secret_arn      = aws_secretsmanager_secret.rds_password.arn
  }
  sensitive = false
}

# ===== EFS File System Outputs =====
output "efs_file_system" {
  description = "EFS file system information"
  value = {
    id         = aws_efs_file_system.reference_db.id
    arn        = aws_efs_file_system.reference_db.arn
    dns_name   = aws_efs_file_system.reference_db.dns_name
    mount_targets = {
      for mt in aws_efs_mount_target.reference_db : mt.availability_zone => {
        id         = mt.id
        ip_address = mt.ip_address
        dns_name   = mt.mount_target_dns_name
      }
    }
    access_points = {
      pmda_db = {
        id  = aws_efs_access_point.pmda_db.id
        arn = aws_efs_access_point.pmda_db.arn
      }
      kraken2_db = {
        id  = aws_efs_access_point.kraken2_db.id
        arn = aws_efs_access_point.kraken2_db.arn
      }
      blast_db = {
        id  = aws_efs_access_point.blast_db.id
        arn = aws_efs_access_point.blast_db.arn
      }
      host_genome = {
        id  = aws_efs_access_point.host_genome.id
        arn = aws_efs_access_point.host_genome.arn
      }
      perv_db = {
        id  = aws_efs_access_point.perv_db.id
        arn = aws_efs_access_point.perv_db.arn
      }
      dorado_models = {
        id  = aws_efs_access_point.dorado_models.id
        arn = aws_efs_access_point.dorado_models.arn
      }
    }
  }
}

# ===== Lambda Functions Outputs =====
output "lambda_functions" {
  description = "Lambda function information"
  value = {
    orchestrator = {
      name         = aws_lambda_function.orchestrator.function_name
      arn          = aws_lambda_function.orchestrator.arn
      invoke_arn   = aws_lambda_function.orchestrator.invoke_arn
      version      = aws_lambda_function.orchestrator.version
    }
    phase_triggers = {
      phase1 = {
        name = aws_lambda_function.phase1_trigger.function_name
        arn  = aws_lambda_function.phase1_trigger.arn
      }
      phase2 = {
        name = aws_lambda_function.phase2_trigger.function_name
        arn  = aws_lambda_function.phase2_trigger.arn
      }
      phase3 = {
        name = aws_lambda_function.phase3_trigger.function_name
        arn  = aws_lambda_function.phase3_trigger.arn
      }
      phase4 = {
        name = aws_lambda_function.phase4_trigger.function_name
        arn  = aws_lambda_function.phase4_trigger.arn
      }
      phase5 = {
        name = aws_lambda_function.phase5_trigger.function_name
        arn  = aws_lambda_function.phase5_trigger.arn
      }
      phase6 = {
        name = aws_lambda_function.phase6_trigger.function_name
        arn  = aws_lambda_function.phase6_trigger.arn
      }
    }
    validators = {
      validate_input = {
        name = aws_lambda_function.validate_input.function_name
        arn  = aws_lambda_function.validate_input.arn
      }
      check_qc = {
        name = aws_lambda_function.check_qc.function_name
        arn  = aws_lambda_function.check_qc.arn
      }
    }
    metadata = {
      update_metadata = {
        name = aws_lambda_function.update_metadata.function_name
        arn  = aws_lambda_function.update_metadata.arn
      }
      notify = {
        name = aws_lambda_function.notify.function_name
        arn  = aws_lambda_function.notify.arn
      }
    }
  }
}

# ===== EC2 Launch Templates Outputs =====
output "ec2_launch_templates" {
  description = "EC2 launch template information"
  value = {
    basecalling = {
      id              = aws_launch_template.basecalling.id
      name            = aws_launch_template.basecalling.name
      latest_version  = aws_launch_template.basecalling.latest_version
      ami_id          = aws_launch_template.basecalling.image_id
    }
    analysis = {
      id              = aws_launch_template.analysis.id
      name            = aws_launch_template.analysis.name
      latest_version  = aws_launch_template.analysis.latest_version
      ami_id          = aws_launch_template.analysis.image_id
    }
  }
}

# ===== IAM Roles Outputs =====
output "iam_roles" {
  description = "IAM role ARNs"
  value = {
    lambda_execution_role = {
      name = aws_iam_role.lambda_execution.name
      arn  = aws_iam_role.lambda_execution.arn
    }
    ec2_instance_role = {
      name = aws_iam_role.ec2_instance.name
      arn  = aws_iam_role.ec2_instance.arn
    }
    ec2_instance_profile = {
      name = aws_iam_instance_profile.ec2_instance.name
      arn  = aws_iam_instance_profile.ec2_instance.arn
    }
  }
}

# ===== SNS Topics Outputs =====
output "sns_topics" {
  description = "SNS topic information"
  value = {
    alerts = {
      name = aws_sns_topic.alerts.name
      arn  = aws_sns_topic.alerts.arn
    }
    critical_alerts = {
      name = aws_sns_topic.critical_alerts.name
      arn  = aws_sns_topic.critical_alerts.arn
    }
    pipeline_status = {
      name = aws_sns_topic.pipeline_status.name
      arn  = aws_sns_topic.pipeline_status.arn
    }
  }
}

# ===== Security Groups Outputs =====
output "security_groups" {
  description = "Security group IDs"
  value = {
    lambda      = aws_security_group.lambda.id
    ec2         = aws_security_group.ec2_instances.id
    rds         = aws_security_group.rds.id
    efs         = aws_security_group.efs.id
  }
}

# ===== CloudWatch Resources Outputs =====
output "cloudwatch_resources" {
  description = "CloudWatch resource information"
  value = {
    dashboard_url = "https://console.aws.amazon.com/cloudwatch/home?region=${var.aws_region}#dashboards:name=${aws_cloudwatch_dashboard.main.dashboard_name}"
    log_groups = {
      lambda      = { for k, v in aws_cloudwatch_log_group.lambda_logs : k => v.name }
      ec2         = { for k, v in aws_cloudwatch_log_group.ec2_logs : k => v.name }
      application = aws_cloudwatch_log_group.application_logs.name
    }
    alarms = {
      critical_pathogen = aws_cloudwatch_metric_alarm.critical_pathogen_detected.arn
      perv_detected     = aws_cloudwatch_metric_alarm.perv_detected.arn
      qc_failure        = aws_cloudwatch_metric_alarm.qc_failure.arn
      high_failure_rate = aws_cloudwatch_metric_alarm.high_failure_rate.arn
    }
  }
}

# ===== EventBridge Resources Outputs =====
output "eventbridge_resources" {
  description = "EventBridge resource information"
  value = {
    event_bus = {
      name = aws_cloudwatch_event_bus.pipeline.name
      arn  = aws_cloudwatch_event_bus.pipeline.arn
    }
    rules = {
      s3_upload         = aws_cloudwatch_event_rule.s3_raw_data_uploaded.name
      phase_completed   = aws_cloudwatch_event_rule.phase_completed.name
      phase_failed      = aws_cloudwatch_event_rule.phase_failed.name
      critical_pathogen = aws_cloudwatch_event_rule.critical_pathogen_detected.name
      perv_detected     = aws_cloudwatch_event_rule.perv_detected.name
      pipeline_complete = aws_cloudwatch_event_rule.pipeline_complete.name
    }
    scheduled_rules = {
      daily_check        = aws_cloudwatch_event_rule.daily_system_check.name
      weekly_maintenance = aws_cloudwatch_event_rule.weekly_maintenance.name
      monthly_db_update  = aws_cloudwatch_event_rule.monthly_db_update.name
    }
  }
}

# ===== Environment and Deployment Information =====
output "deployment_info" {
  description = "Deployment configuration information"
  value = {
    environment    = var.environment
    project_name   = var.project_name
    aws_region     = var.aws_region
    aws_account_id = data.aws_caller_identity.current.account_id
    deployment_timestamp = timestamp()
  }
}

# ===== Pipeline Configuration URLs (for easy access) =====
output "console_urls" {
  description = "AWS Console URLs for quick access"
  value = {
    s3_console       = "https://s3.console.aws.amazon.com/s3/buckets/${var.project_name}"
    rds_console      = "https://console.aws.amazon.com/rds/home?region=${var.aws_region}#database:id=${aws_db_instance.metadata.id}"
    lambda_console   = "https://console.aws.amazon.com/lambda/home?region=${var.aws_region}#/functions"
    cloudwatch_logs  = "https://console.aws.amazon.com/cloudwatch/home?region=${var.aws_region}#logsV2:log-groups"
    ec2_console      = "https://console.aws.amazon.com/ec2/v2/home?region=${var.aws_region}#Instances:"
    efs_console      = "https://console.aws.amazon.com/efs/home?region=${var.aws_region}#/file-systems/${aws_efs_file_system.reference_db.id}"
    sns_console      = "https://console.aws.amazon.com/sns/v3/home?region=${var.aws_region}#/topics"
    eventbridge_console = "https://console.aws.amazon.com/events/home?region=${var.aws_region}#/eventbuses/${aws_cloudwatch_event_bus.pipeline.name}"
  }
}

# ===== Connection Strings and Commands =====
output "connection_commands" {
  description = "Useful connection commands"
  value = {
    rds_connection = "postgresql://${aws_db_instance.metadata.username}:<password>@${aws_db_instance.metadata.endpoint}/${aws_db_instance.metadata.db_name}"
    efs_mount = "sudo mount -t efs -o tls ${aws_efs_file_system.reference_db.id}:/ /mnt/efs"
    s3_sync_raw = "aws s3 sync ./data s3://${aws_s3_bucket.raw_data.id}/"
    s3_sync_results = "aws s3 sync s3://${aws_s3_bucket.analysis.id}/reports/ ./reports/"
  }
  sensitive = true
}

# ===== Cost Tracking Tags =====
output "resource_tags" {
  description = "Common resource tags applied"
  value = {
    Environment = var.environment
    Project     = var.project_name
    ManagedBy   = "Terraform"
    Purpose     = "MinION Metagenomics Pipeline"
    Compliance  = "PMDA"
  }
}

# ===== Pipeline Entry Points =====
output "pipeline_entry_points" {
  description = "Entry points for pipeline operations"
  value = {
    upload_bucket     = "s3://${aws_s3_bucket.raw_data.id}/fast5/"
    orchestrator_arn  = aws_lambda_function.orchestrator.arn
    api_endpoint      = "Not yet configured - will be added with API Gateway"
    reports_bucket    = "s3://${aws_s3_bucket.analysis.id}/reports/"
  }
}

# ===== Monitoring and Alerting =====
output "monitoring_endpoints" {
  description = "Monitoring and alerting endpoints"
  value = {
    alerts_topic         = aws_sns_topic.alerts.arn
    critical_alerts_topic = aws_sns_topic.critical_alerts.arn
    dashboard_name       = aws_cloudwatch_dashboard.main.dashboard_name
    log_insights_group   = aws_cloudwatch_log_group.application_logs.name
  }
}

# ===== Database Access Information =====
output "database_access" {
  description = "Database access information for applications"
  value = {
    host        = aws_db_instance.metadata.address
    port        = aws_db_instance.metadata.port
    database    = aws_db_instance.metadata.db_name
    username    = aws_db_instance.metadata.username
    secret_name = aws_secretsmanager_secret.rds_password.name
    connection_string = format(
      "postgresql://%s:<get_from_secrets_manager>@%s:%s/%s",
      aws_db_instance.metadata.username,
      aws_db_instance.metadata.address,
      aws_db_instance.metadata.port,
      aws_db_instance.metadata.db_name
    )
  }
  sensitive = true
}

# ===== Reference Database Mount Points =====
output "reference_databases" {
  description = "EFS mount points for reference databases"
  value = {
    pmda_db       = "/mnt/efs/pmda"
    kraken2_db    = "/mnt/efs/kraken2"
    blast_db      = "/mnt/efs/blast"
    host_genome   = "/mnt/efs/host_genome"
    perv_db       = "/mnt/efs/perv"
    dorado_models = "/mnt/efs/dorado"
  }
}

# ===== Summary Output for Documentation =====
output "infrastructure_summary" {
  description = "Summary of deployed infrastructure"
  value = <<-EOT
    MinION Metagenomics Pipeline Infrastructure
    ============================================
    Environment: ${var.environment}
    Region: ${var.aws_region}

    Core Components:
    - VPC with ${length(module.vpc.private_subnets)} private subnets
    - RDS PostgreSQL database (${var.rds_instance_class})
    - EFS file system with 6 access points
    - 11 Lambda functions for orchestration
    - 3 S3 buckets for data storage
    - EC2 launch templates for GPU and CPU instances

    Monitoring:
    - CloudWatch dashboard: ${aws_cloudwatch_dashboard.main.dashboard_name}
    - SNS alert topics: 3 configured
    - EventBridge rules: ${length(aws_cloudwatch_event_rule.s3_raw_data_uploaded.name) + 6} active rules

    Security:
    - All data encrypted at rest
    - VPC isolation for compute resources
    - IAM roles with least privilege
    - Secrets Manager for credentials

    Compliance:
    - PMDA 5-year data retention configured
    - Audit logging enabled
    - Data integrity checksums implemented
  EOT
}