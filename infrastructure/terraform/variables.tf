# MinION Metagenomics Pipeline - Terraform Variables
# PMDA 91 Pathogen Screening System

variable "aws_region" {
  description = "AWS region for deployment"
  type        = string
  default     = "ap-northeast-1"  # Tokyo
}

variable "environment" {
  description = "Environment name (dev, staging, prod)"
  type        = string
  default     = "prod"
}

variable "project_name" {
  description = "Project name prefix"
  type        = string
  default     = "minion-metagenomics"
}

# VPC Configuration
variable "vpc_cidr" {
  description = "CIDR block for VPC"
  type        = string
  default     = "10.0.0.0/16"
}

variable "availability_zones" {
  description = "Availability zones"
  type        = list(string)
  default     = ["ap-northeast-1a", "ap-northeast-1c"]
}

# EC2 Configuration
variable "basecalling_ami_id" {
  description = "AMI ID for basecalling (GPU instance)"
  type        = string
  default     = ""  # To be filled after AMI creation
}

variable "analysis_ami_id" {
  description = "AMI ID for analysis (CPU instance)"
  type        = string
  default     = ""  # To be filled after AMI creation
}

variable "basecalling_instance_type" {
  description = "Instance type for basecalling (GPU required)"
  type        = string
  default     = "g5.xlarge"  # NVIDIA A10G GPU, 4 vCPU, 16GB RAM
}

variable "qc_instance_type" {
  description = "Instance type for QC analysis"
  type        = string
  default     = "c6i.xlarge"  # 4 vCPU, 8GB RAM
}

variable "host_removal_instance_type" {
  description = "Instance type for host genome removal"
  type        = string
  default     = "c6i.4xlarge"  # 16 vCPU, 32GB RAM
}

variable "kraken2_instance_type" {
  description = "Instance type for Kraken2 classification"
  type        = string
  default     = "r6i.xlarge"  # 4 vCPU, 32GB RAM (memory-optimized)
}

variable "blast_instance_type" {
  description = "Instance type for BLAST"
  type        = string
  default     = "c6i.8xlarge"  # 32 vCPU, 64GB RAM
}

variable "assembly_instance_type" {
  description = "Instance type for de novo assembly"
  type        = string
  default     = "c6i.8xlarge"  # 32 vCPU, 64GB RAM
}

variable "perv_instance_type" {
  description = "Instance type for PERV analysis"
  type        = string
  default     = "c6i.xlarge"  # 4 vCPU, 8GB RAM
}

variable "quantification_instance_type" {
  description = "Instance type for quantification"
  type        = string
  default     = "c6i.xlarge"  # 4 vCPU, 8GB RAM
}

variable "reports_instance_type" {
  description = "Instance type for report generation"
  type        = string
  default     = "c6i.xlarge"  # 4 vCPU, 8GB RAM
}

# RDS Configuration
variable "rds_instance_class" {
  description = "RDS instance class"
  type        = string
  default     = "db.t3.medium"  # 2 vCPU, 4GB RAM
}

variable "rds_allocated_storage" {
  description = "RDS allocated storage (GB)"
  type        = number
  default     = 100
}

variable "rds_max_allocated_storage" {
  description = "RDS max allocated storage for autoscaling (GB)"
  type        = number
  default     = 500
}

variable "rds_username" {
  description = "RDS master username"
  type        = string
  default     = "admin"
  sensitive   = true
}

variable "rds_database_name" {
  description = "RDS database name"
  type        = string
  default     = "minion_metagenomics"
}

# EFS Configuration
variable "efs_performance_mode" {
  description = "EFS performance mode"
  type        = string
  default     = "generalPurpose"
}

variable "efs_throughput_mode" {
  description = "EFS throughput mode"
  type        = string
  default     = "bursting"
}

# S3 Configuration
variable "s3_lifecycle_glacier_days" {
  description = "Days before transitioning to Glacier Deep Archive"
  type        = number
  default     = 90
}

variable "s3_expiration_days" {
  description = "Days before expiration (5 years for PMDA compliance)"
  type        = number
  default     = 1825  # 5 years
}

# Lambda Configuration
variable "lambda_runtime" {
  description = "Lambda runtime version"
  type        = string
  default     = "python3.11"
}

variable "lambda_timeout" {
  description = "Lambda timeout in seconds"
  type        = number
  default     = 900  # 15 minutes
}

variable "lambda_memory_size" {
  description = "Lambda memory size in MB"
  type        = number
  default     = 512
}

# SNS Configuration
variable "notification_email" {
  description = "Email address for SNS notifications"
  type        = string
  default     = ""  # To be filled by user
}

# Tags
variable "tags" {
  description = "Common tags for all resources"
  type        = map(string)
  default = {
    Project     = "MinION-Metagenomics"
    Environment = "Production"
    Compliance  = "PMDA"
    ManagedBy   = "Terraform"
  }
}
