# MinION Metagenomics Pipeline - Main Terraform Configuration
# PMDA 91 Pathogen Screening System

terraform {
  required_version = ">= 1.5.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
    random = {
      source  = "hashicorp/random"
      version = "~> 3.5"
    }
  }

  # Backend configuration for remote state (optional)
  # backend "s3" {
  #   bucket         = "minion-terraform-state"
  #   key            = "metagenomics/terraform.tfstate"
  #   region         = "ap-northeast-1"
  #   dynamodb_table = "terraform-locks"
  #   encrypt        = true
  # }
}

provider "aws" {
  region = var.aws_region

  default_tags {
    tags = var.tags
  }
}

# Data sources
data "aws_caller_identity" "current" {}

data "aws_region" "current" {}

# Random password for RDS
resource "random_password" "rds_password" {
  length  = 32
  special = true
  override_special = "!#$%&*()-_=+[]{}<>:?"
}

# VPC Module
module "vpc" {
  source  = "terraform-aws-modules/vpc/aws"
  version = "~> 5.0"

  name = "${var.project_name}-vpc"
  cidr = var.vpc_cidr

  azs             = var.availability_zones
  private_subnets = [
    cidrsubnet(var.vpc_cidr, 8, 1),  # 10.0.1.0/24
    cidrsubnet(var.vpc_cidr, 8, 2),  # 10.0.2.0/24
  ]
  public_subnets = [
    cidrsubnet(var.vpc_cidr, 8, 101),  # 10.0.101.0/24
    cidrsubnet(var.vpc_cidr, 8, 102),  # 10.0.102.0/24
  ]
  database_subnets = [
    cidrsubnet(var.vpc_cidr, 8, 201),  # 10.0.201.0/24
    cidrsubnet(var.vpc_cidr, 8, 202),  # 10.0.202.0/24
  ]

  enable_nat_gateway   = true
  single_nat_gateway   = false  # High availability
  enable_dns_hostnames = true
  enable_dns_support   = true

  # VPC Flow Logs for audit
  enable_flow_log                      = true
  create_flow_log_cloudwatch_iam_role  = true
  create_flow_log_cloudwatch_log_group = true

  tags = {
    Name = "${var.project_name}-vpc"
  }
}

# Security Groups
resource "aws_security_group" "ec2_instances" {
  name_description = "Security group for MinION EC2 instances"
  vpc_id      = module.vpc.vpc_id

  # Outbound: Allow all (for S3, RDS, EFS access)
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound traffic"
  }

  tags = {
    Name = "${var.project_name}-ec2-sg"
  }
}

resource "aws_security_group" "rds" {
  name_description = "Security group for RDS PostgreSQL"
  vpc_id      = module.vpc.vpc_id

  # Inbound: PostgreSQL from EC2 instances and Lambda
  ingress {
    from_port       = 5432
    to_port         = 5432
    protocol        = "tcp"
    security_groups = [aws_security_group.ec2_instances.id, aws_security_group.lambda.id]
    description     = "PostgreSQL from EC2 and Lambda"
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound traffic"
  }

  tags = {
    Name = "${var.project_name}-rds-sg"
  }
}

resource "aws_security_group" "efs" {
  name_description = "Security group for EFS"
  vpc_id      = module.vpc.vpc_id

  # Inbound: NFS from EC2 instances
  ingress {
    from_port       = 2049
    to_port         = 2049
    protocol        = "tcp"
    security_groups = [aws_security_group.ec2_instances.id]
    description     = "NFS from EC2"
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound traffic"
  }

  tags = {
    Name = "${var.project_name}-efs-sg"
  }
}

resource "aws_security_group" "lambda" {
  name_description = "Security group for Lambda functions"
  vpc_id      = module.vpc.vpc_id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound traffic"
  }

  tags = {
    Name = "${var.project_name}-lambda-sg"
  }
}

# VPC Endpoints for S3 (cost optimization, avoid NAT charges)
resource "aws_vpc_endpoint" "s3" {
  vpc_id       = module.vpc.vpc_id
  service_name = "com.amazonaws.${var.aws_region}.s3"

  route_table_ids = concat(
    module.vpc.private_route_table_ids,
    module.vpc.public_route_table_ids
  )

  tags = {
    Name = "${var.project_name}-s3-endpoint"
  }
}

# Secrets Manager for RDS password
resource "aws_secretsmanager_secret" "rds_password" {
  name        = "${var.project_name}-rds-password"
  description = "RDS PostgreSQL master password"

  recovery_window_in_days = 7
}

resource "aws_secretsmanager_secret_version" "rds_password" {
  secret_id     = aws_secretsmanager_secret.rds_password.id
  secret_string = jsonencode({
    username = var.rds_username
    password = random_password.rds_password.result
    engine   = "postgres"
    host     = aws_db_instance.metadata.address
    port     = 5432
    dbname   = var.rds_database_name
  })
}

# CloudWatch Log Group for VPC Flow Logs
resource "aws_cloudwatch_log_group" "lambda_logs" {
  name              = "/aws/lambda/${var.project_name}"
  retention_in_days = 30

  tags = {
    Name = "${var.project_name}-lambda-logs"
  }
}

resource "aws_cloudwatch_log_group" "ec2_logs" {
  name              = "/aws/ec2/${var.project_name}"
  retention_in_days = 90  # PMDA compliance: retain logs for audit

  tags = {
    Name = "${var.project_name}-ec2-logs"
  }
}

# Outputs
output "vpc_id" {
  description = "VPC ID"
  value       = module.vpc.vpc_id
}

output "private_subnet_ids" {
  description = "Private subnet IDs"
  value       = module.vpc.private_subnets
}

output "public_subnet_ids" {
  description = "Public subnet IDs"
  value       = module.vpc.public_subnets
}

output "database_subnet_ids" {
  description = "Database subnet IDs"
  value       = module.vpc.database_subnets
}

output "ec2_security_group_id" {
  description = "EC2 instances security group ID"
  value       = aws_security_group.ec2_instances.id
}

output "rds_security_group_id" {
  description = "RDS security group ID"
  value       = aws_security_group.rds.id
}

output "efs_security_group_id" {
  description = "EFS security group ID"
  value       = aws_security_group.efs.id
}

output "lambda_security_group_id" {
  description = "Lambda security group ID"
  value       = aws_security_group.lambda.id
}

output "rds_password_secret_arn" {
  description = "ARN of the RDS password secret"
  value       = aws_secretsmanager_secret.rds_password.arn
  sensitive   = true
}
