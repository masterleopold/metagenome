# MinION Metagenomics Pipeline - RDS PostgreSQL Configuration
# Database for metadata, pathogen quantification results, and audit logs

# DB Subnet Group
resource "aws_db_subnet_group" "metadata" {
  name       = "${var.project_name}-db-subnet-group"
  subnet_ids = module.vpc.database_subnets

  tags = {
    Name        = "${var.project_name}-db-subnet-group"
    Environment = var.environment
  }
}

# RDS PostgreSQL Instance
resource "aws_db_instance" "metadata" {
  identifier     = "${var.project_name}-metadata-db"
  engine         = "postgres"
  engine_version = "15.4"

  instance_class        = var.rds_instance_class
  allocated_storage     = var.rds_allocated_storage
  max_allocated_storage = var.rds_max_allocated_storage
  storage_type          = "gp3"
  storage_encrypted     = true

  db_name  = var.rds_database_name
  username = var.rds_username
  password = random_password.rds_password.result

  vpc_security_group_ids = [aws_security_group.rds.id]
  db_subnet_group_name   = aws_db_subnet_group.metadata.name

  # Backup configuration for PMDA compliance
  backup_retention_period = 30
  backup_window          = "03:00-04:00"
  maintenance_window     = "Mon:04:00-Mon:05:00"

  # High availability
  multi_az               = true

  # Performance and monitoring
  enabled_cloudwatch_logs_exports = ["postgresql"]
  performance_insights_enabled    = true
  performance_insights_retention_period = 7
  monitoring_interval            = 60
  monitoring_role_arn           = aws_iam_role.rds_monitoring.arn

  # Security
  deletion_protection = true
  skip_final_snapshot = false
  final_snapshot_identifier = "${var.project_name}-metadata-db-final-snapshot-${formatdate("YYYY-MM-DD-hhmm", timestamp())}"

  # Auto minor version upgrade for security patches
  auto_minor_version_upgrade = true
  apply_immediately         = false

  tags = {
    Name        = "${var.project_name}-metadata-db"
    Environment = var.environment
    Purpose     = "MinION analysis metadata and results"
    Compliance  = "PMDA"
  }
}

# Parameter Group for PostgreSQL tuning
resource "aws_db_parameter_group" "postgres15" {
  name   = "${var.project_name}-postgres15-params"
  family = "postgres15"

  parameter {
    name  = "shared_preload_libraries"
    value = "pg_stat_statements,pgaudit"
  }

  parameter {
    name  = "log_statement"
    value = "all"
  }

  parameter {
    name  = "log_min_duration_statement"
    value = "1000"  # Log queries longer than 1 second
  }

  parameter {
    name  = "pgaudit.log"
    value = "ALL"
  }

  parameter {
    name  = "max_connections"
    value = "200"
  }

  parameter {
    name  = "shared_buffers"
    value = "{DBInstanceClassMemory/4}"
  }

  parameter {
    name  = "effective_cache_size"
    value = "{DBInstanceClassMemory*3/4}"
  }

  parameter {
    name  = "work_mem"
    value = "16384"  # 16MB
  }

  parameter {
    name  = "maintenance_work_mem"
    value = "524288"  # 512MB
  }

  tags = {
    Name        = "${var.project_name}-postgres15-params"
    Environment = var.environment
  }
}

# Apply parameter group to RDS instance
resource "aws_db_instance" "metadata_params" {
  identifier = aws_db_instance.metadata.id

  parameter_group_name = aws_db_parameter_group.postgres15.name

  depends_on = [aws_db_instance.metadata]

  lifecycle {
    ignore_changes = [identifier, allocated_storage, instance_class]
  }
}

# IAM Role for RDS Enhanced Monitoring
resource "aws_iam_role" "rds_monitoring" {
  name = "${var.project_name}-rds-monitoring-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "monitoring.rds.amazonaws.com"
        }
      }
    ]
  })

  tags = {
    Name        = "${var.project_name}-rds-monitoring-role"
    Environment = var.environment
  }
}

# Attach AWS managed policy for RDS Enhanced Monitoring
resource "aws_iam_role_policy_attachment" "rds_monitoring" {
  role       = aws_iam_role.rds_monitoring.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonRDSEnhancedMonitoringRole"
}

# CloudWatch Alarms for RDS
resource "aws_cloudwatch_metric_alarm" "rds_cpu" {
  alarm_name          = "${var.project_name}-rds-high-cpu"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "CPUUtilization"
  namespace           = "AWS/RDS"
  period              = "300"
  statistic           = "Average"
  threshold           = "80"
  alarm_description   = "This metric monitors RDS CPU utilization"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    DBInstanceIdentifier = aws_db_instance.metadata.id
  }

  tags = {
    Name        = "${var.project_name}-rds-cpu-alarm"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_metric_alarm" "rds_storage" {
  alarm_name          = "${var.project_name}-rds-low-storage"
  comparison_operator = "LessThanThreshold"
  evaluation_periods  = "1"
  metric_name         = "FreeStorageSpace"
  namespace           = "AWS/RDS"
  period              = "300"
  statistic           = "Average"
  threshold           = "10737418240"  # 10 GB in bytes
  alarm_description   = "Alert when RDS free storage is below 10GB"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    DBInstanceIdentifier = aws_db_instance.metadata.id
  }

  tags = {
    Name        = "${var.project_name}-rds-storage-alarm"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_metric_alarm" "rds_connections" {
  alarm_name          = "${var.project_name}-rds-high-connections"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "DatabaseConnections"
  namespace           = "AWS/RDS"
  period              = "300"
  statistic           = "Average"
  threshold           = "180"  # 90% of max_connections (200)
  alarm_description   = "Alert when database connections exceed 90% of maximum"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    DBInstanceIdentifier = aws_db_instance.metadata.id
  }

  tags = {
    Name        = "${var.project_name}-rds-connections-alarm"
    Environment = var.environment
  }
}

# RDS Proxy (Optional - for better connection management)
resource "aws_db_proxy" "metadata" {
  count = var.environment == "prod" ? 1 : 0  # Only create in production

  name                   = "${var.project_name}-rds-proxy"
  engine_family          = "POSTGRESQL"
  auth {
    auth_scheme = "SECRETS"
    secret_arn  = aws_secretsmanager_secret.rds_password.arn
  }

  role_arn               = aws_iam_role.rds_proxy[0].arn
  vpc_subnet_ids         = module.vpc.database_subnets

  max_connections_percent        = 100
  max_idle_connections_percent   = 50
  connection_borrow_timeout      = 120

  require_tls = true

  target {
    db_instance_identifier = aws_db_instance.metadata.id
  }

  tags = {
    Name        = "${var.project_name}-rds-proxy"
    Environment = var.environment
  }
}

# IAM Role for RDS Proxy
resource "aws_iam_role" "rds_proxy" {
  count = var.environment == "prod" ? 1 : 0

  name = "${var.project_name}-rds-proxy-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "rds.amazonaws.com"
        }
      }
    ]
  })
}

# IAM Policy for RDS Proxy to access Secrets Manager
resource "aws_iam_role_policy" "rds_proxy" {
  count = var.environment == "prod" ? 1 : 0

  name = "${var.project_name}-rds-proxy-policy"
  role = aws_iam_role.rds_proxy[0].id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "secretsmanager:GetSecretValue",
          "secretsmanager:DescribeSecret"
        ]
        Resource = aws_secretsmanager_secret.rds_password.arn
      },
      {
        Effect = "Allow"
        Action = [
          "kms:Decrypt"
        ]
        Resource = "*"
        Condition = {
          StringEquals = {
            "kms:ViaService" = "secretsmanager.${var.aws_region}.amazonaws.com"
          }
        }
      }
    ]
  })
}

# Outputs for RDS
output "rds_endpoint" {
  description = "RDS instance endpoint"
  value       = aws_db_instance.metadata.endpoint
  sensitive   = true
}

output "rds_port" {
  description = "RDS instance port"
  value       = aws_db_instance.metadata.port
}

output "rds_database_name" {
  description = "RDS database name"
  value       = aws_db_instance.metadata.db_name
}

output "rds_proxy_endpoint" {
  description = "RDS Proxy endpoint (if created)"
  value       = try(aws_db_proxy.metadata[0].endpoint, "")
  sensitive   = true
}