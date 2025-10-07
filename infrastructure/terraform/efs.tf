# MinION Metagenomics Pipeline - EFS Configuration
# Elastic File System for reference databases and shared data

# EFS File System
resource "aws_efs_file_system" "reference_db" {
  creation_token = "${var.project_name}-reference-db"
  encrypted      = true

  performance_mode = var.efs_performance_mode
  throughput_mode  = var.efs_throughput_mode

  # Enable automatic backups
  lifecycle_policy {
    transition_to_ia = "AFTER_30_DAYS"
  }

  lifecycle_policy {
    transition_to_primary_storage_class = "AFTER_1_ACCESS"
  }

  tags = {
    Name        = "${var.project_name}-reference-db"
    Environment = var.environment
    Purpose     = "Reference genomes and PMDA pathogen databases"
  }
}

# EFS Mount Targets - one per availability zone
resource "aws_efs_mount_target" "reference_db" {
  count = length(var.availability_zones)

  file_system_id  = aws_efs_file_system.reference_db.id
  subnet_id       = module.vpc.private_subnets[count.index]
  security_groups = [aws_security_group.efs.id]
}

# EFS Access Points for different database types
resource "aws_efs_access_point" "pmda_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/pmda_91_pathogens"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-pmda-db-access-point"
    Environment = var.environment
    Database    = "PMDA 91 pathogens"
  }
}

resource "aws_efs_access_point" "kraken2_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/kraken2_db"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-kraken2-db-access-point"
    Environment = var.environment
    Database    = "Kraken2 reference database"
  }
}

resource "aws_efs_access_point" "blast_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/blast_db"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-blast-db-access-point"
    Environment = var.environment
    Database    = "BLAST nucleotide database"
  }
}

resource "aws_efs_access_point" "host_genome" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/host_genome"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-host-genome-access-point"
    Environment = var.environment
    Database    = "Sus scrofa genome"
  }
}

resource "aws_efs_access_point" "perv_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/perv_references"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-perv-db-access-point"
    Environment = var.environment
    Database    = "PERV reference sequences"
  }
}

resource "aws_efs_access_point" "dorado_models" {
  file_system_id = aws_efs_file_system.reference_db.id

  root_directory {
    path = "/dorado_models"
    creation_info {
      owner_uid   = 1000
      owner_gid   = 1000
      permissions = "755"
    }
  }

  posix_user {
    uid = 1000
    gid = 1000
  }

  tags = {
    Name        = "${var.project_name}-dorado-models-access-point"
    Environment = var.environment
    Database    = "Dorado basecalling models"
  }
}

# EFS Backup Policy
resource "aws_efs_backup_policy" "reference_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  backup_policy {
    status = "ENABLED"
  }
}

# CloudWatch Alarms for EFS
resource "aws_cloudwatch_metric_alarm" "efs_client_connections" {
  alarm_name          = "${var.project_name}-efs-high-connections"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "ClientConnections"
  namespace           = "AWS/EFS"
  period              = "300"
  statistic           = "Sum"
  threshold           = "100"
  alarm_description   = "Alert when EFS client connections are high"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    FileSystemId = aws_efs_file_system.reference_db.id
  }

  tags = {
    Name        = "${var.project_name}-efs-connections-alarm"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_metric_alarm" "efs_burst_credits" {
  alarm_name          = "${var.project_name}-efs-low-burst-credits"
  comparison_operator = "LessThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "BurstCreditBalance"
  namespace           = "AWS/EFS"
  period              = "300"
  statistic           = "Average"
  threshold           = "1000000000000"  # 1 TiB worth of burst credits
  alarm_description   = "Alert when EFS burst credits are low"
  alarm_actions       = [aws_sns_topic.alerts.arn]

  dimensions = {
    FileSystemId = aws_efs_file_system.reference_db.id
  }

  tags = {
    Name        = "${var.project_name}-efs-burst-credits-alarm"
    Environment = var.environment
  }
}

# IAM Policy for EFS access
resource "aws_iam_policy" "efs_access" {
  name        = "${var.project_name}-efs-access-policy"
  description = "Policy for EC2 instances to access EFS"

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "elasticfilesystem:ClientMount",
          "elasticfilesystem:ClientWrite",
          "elasticfilesystem:ClientRootAccess"
        ]
        Resource = aws_efs_file_system.reference_db.arn
      },
      {
        Effect = "Allow"
        Action = [
          "elasticfilesystem:DescribeAccessPoints",
          "elasticfilesystem:DescribeFileSystems",
          "elasticfilesystem:DescribeMountTargets"
        ]
        Resource = "*"
      }
    ]
  })

  tags = {
    Name        = "${var.project_name}-efs-access-policy"
    Environment = var.environment
  }
}

# EFS File System Policy
resource "aws_efs_file_system_policy" "reference_db" {
  file_system_id = aws_efs_file_system.reference_db.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "EnforceInTransitEncryption"
        Effect = "Deny"
        Principal = {
          AWS = "*"
        }
        Action = "*"
        Resource = aws_efs_file_system.reference_db.arn
        Condition = {
          Bool = {
            "aws:SecureTransport" = "false"
          }
        }
      },
      {
        Sid    = "AllowEC2Access"
        Effect = "Allow"
        Principal = {
          AWS = aws_iam_role.ec2_instance.arn
        }
        Action = [
          "elasticfilesystem:ClientMount",
          "elasticfilesystem:ClientWrite"
        ]
        Resource = aws_efs_file_system.reference_db.arn
      }
    ]
  })
}

# Outputs for EFS
output "efs_file_system_id" {
  description = "ID of the EFS file system"
  value       = aws_efs_file_system.reference_db.id
}

output "efs_file_system_dns_name" {
  description = "DNS name of the EFS file system"
  value       = aws_efs_file_system.reference_db.dns_name
}

output "efs_mount_target_ids" {
  description = "IDs of the EFS mount targets"
  value       = aws_efs_mount_target.reference_db[*].id
}

output "efs_access_point_pmda" {
  description = "Access point ID for PMDA database"
  value       = aws_efs_access_point.pmda_db.id
}

output "efs_access_point_kraken2" {
  description = "Access point ID for Kraken2 database"
  value       = aws_efs_access_point.kraken2_db.id
}

output "efs_access_point_blast" {
  description = "Access point ID for BLAST database"
  value       = aws_efs_access_point.blast_db.id
}

output "efs_access_point_host_genome" {
  description = "Access point ID for host genome"
  value       = aws_efs_access_point.host_genome.id
}

output "efs_access_point_perv" {
  description = "Access point ID for PERV references"
  value       = aws_efs_access_point.perv_db.id
}

output "efs_access_point_dorado" {
  description = "Access point ID for Dorado models"
  value       = aws_efs_access_point.dorado_models.id
}