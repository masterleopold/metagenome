# MinION Metagenomics Pipeline - S3 Configuration
# S3 buckets for raw data, analysis results, and scripts

# Raw data bucket - stores FAST5 files from MinION
resource "aws_s3_bucket" "raw_data" {
  bucket = "${var.project_name}-raw-data-${var.environment}"

  tags = {
    Name        = "${var.project_name}-raw-data"
    Environment = var.environment
    Purpose     = "MinION FAST5 raw data storage"
  }
}

# Analysis bucket - stores intermediate and final results
resource "aws_s3_bucket" "analysis" {
  bucket = "${var.project_name}-analysis-${var.environment}"

  tags = {
    Name        = "${var.project_name}-analysis"
    Environment = var.environment
    Purpose     = "Analysis results and intermediate files"
  }
}

# Scripts bucket - stores analysis scripts and configurations
resource "aws_s3_bucket" "scripts" {
  bucket = "${var.project_name}-scripts-${var.environment}"

  tags = {
    Name        = "${var.project_name}-scripts"
    Environment = var.environment
    Purpose     = "Analysis scripts and pipeline code"
  }
}

# Enable versioning for all buckets
resource "aws_s3_bucket_versioning" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_versioning" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_versioning" "scripts" {
  bucket = aws_s3_bucket.scripts.id

  versioning_configuration {
    status = "Enabled"
  }
}

# Server-side encryption for all buckets
resource "aws_s3_bucket_server_side_encryption_configuration" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "scripts" {
  bucket = aws_s3_bucket.scripts.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

# Block public access for all buckets
resource "aws_s3_bucket_public_access_block" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_public_access_block" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_public_access_block" "scripts" {
  bucket = aws_s3_bucket.scripts.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# Lifecycle policies for cost optimization
resource "aws_s3_bucket_lifecycle_configuration" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  rule {
    id     = "archive_old_data"
    status = "Enabled"

    transition {
      days          = var.s3_lifecycle_glacier_days
      storage_class = "GLACIER_IR"
    }

    transition {
      days          = 365
      storage_class = "DEEP_ARCHIVE"
    }

    expiration {
      days = var.s3_expiration_days  # 5 years for PMDA compliance
    }

    noncurrent_version_expiration {
      noncurrent_days = 30
    }
  }

  rule {
    id     = "delete_incomplete_multipart_uploads"
    status = "Enabled"

    abort_incomplete_multipart_upload {
      days_after_initiation = 7
    }
  }
}

resource "aws_s3_bucket_lifecycle_configuration" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  rule {
    id     = "archive_old_results"
    status = "Enabled"

    transition {
      days          = var.s3_lifecycle_glacier_days
      storage_class = "GLACIER_IR"
    }

    transition {
      days          = 365
      storage_class = "DEEP_ARCHIVE"
    }

    expiration {
      days = var.s3_expiration_days  # 5 years for PMDA compliance
    }

    noncurrent_version_expiration {
      noncurrent_days = 30
    }
  }

  rule {
    id     = "delete_temp_files"
    status = "Enabled"

    filter {
      prefix = "temp/"
    }

    expiration {
      days = 7
    }
  }
}

# S3 Event Notifications for triggering Lambda
resource "aws_s3_bucket_notification" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.orchestrator.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "fast5/"
    filter_suffix       = ".fast5"
  }

  depends_on = [aws_lambda_permission.allow_s3_invoke_orchestrator]
}

resource "aws_s3_bucket_notification" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  # Phase 1 completion triggers Phase 2
  lambda_function {
    id                  = "phase1_complete"
    lambda_function_arn = aws_lambda_function.phase2_trigger.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "phase1/"
    filter_suffix       = "basecalling_summary.json"
  }

  # Phase 2 completion triggers Phase 3
  lambda_function {
    id                  = "phase2_complete"
    lambda_function_arn = aws_lambda_function.phase3_trigger.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "phase2/"
    filter_suffix       = "qc_passed.json"
  }

  # Phase 3 completion triggers Phase 4
  lambda_function {
    id                  = "phase3_complete"
    lambda_function_arn = aws_lambda_function.phase4_trigger.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "phase3/"
    filter_suffix       = "host_removal_complete.json"
  }

  # Phase 4 completion triggers Phase 5
  lambda_function {
    id                  = "phase4_complete"
    lambda_function_arn = aws_lambda_function.phase5_trigger.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "phase4/"
    filter_suffix       = "pathogen_detection_complete.json"
  }

  # Phase 5 completion triggers Phase 6
  lambda_function {
    id                  = "phase5_complete"
    lambda_function_arn = aws_lambda_function.phase6_trigger.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "phase5/"
    filter_suffix       = "quantification_complete.json"
  }

  depends_on = [
    aws_lambda_permission.allow_s3_invoke_phase2,
    aws_lambda_permission.allow_s3_invoke_phase3,
    aws_lambda_permission.allow_s3_invoke_phase4,
    aws_lambda_permission.allow_s3_invoke_phase5,
    aws_lambda_permission.allow_s3_invoke_phase6
  ]
}

# CORS configuration for analysis bucket (for web dashboard access)
resource "aws_s3_bucket_cors_configuration" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  cors_rule {
    allowed_headers = ["*"]
    allowed_methods = ["GET", "HEAD"]
    allowed_origins = ["https://quicksight.aws.amazon.com"]
    expose_headers  = ["ETag"]
    max_age_seconds = 3000
  }
}

# S3 bucket policies for controlled access
resource "aws_s3_bucket_policy" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "EnforceSSLRequestsOnly"
        Effect = "Deny"
        Principal = "*"
        Action = "s3:*"
        Resource = [
          aws_s3_bucket.raw_data.arn,
          "${aws_s3_bucket.raw_data.arn}/*"
        ]
        Condition = {
          Bool = {
            "aws:SecureTransport" = "false"
          }
        }
      }
    ]
  })
}

resource "aws_s3_bucket_policy" "analysis" {
  bucket = aws_s3_bucket.analysis.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "EnforceSSLRequestsOnly"
        Effect = "Deny"
        Principal = "*"
        Action = "s3:*"
        Resource = [
          aws_s3_bucket.analysis.arn,
          "${aws_s3_bucket.analysis.arn}/*"
        ]
        Condition = {
          Bool = {
            "aws:SecureTransport" = "false"
          }
        }
      }
    ]
  })
}

resource "aws_s3_bucket_policy" "scripts" {
  bucket = aws_s3_bucket.scripts.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "EnforceSSLRequestsOnly"
        Effect = "Deny"
        Principal = "*"
        Action = "s3:*"
        Resource = [
          aws_s3_bucket.scripts.arn,
          "${aws_s3_bucket.scripts.arn}/*"
        ]
        Condition = {
          Bool = {
            "aws:SecureTransport" = "false"
          }
        }
      }
    ]
  })
}

# Intelligent Tiering configuration for cost optimization
resource "aws_s3_bucket_intelligent_tiering_configuration" "raw_data" {
  bucket = aws_s3_bucket.raw_data.id
  name   = "EntireBucket"

  tiering {
    access_tier = "ARCHIVE_ACCESS"
    days        = 90
  }

  tiering {
    access_tier = "DEEP_ARCHIVE_ACCESS"
    days        = 180
  }
}

resource "aws_s3_bucket_intelligent_tiering_configuration" "analysis" {
  bucket = aws_s3_bucket.analysis.id
  name   = "EntireBucket"

  tiering {
    access_tier = "ARCHIVE_ACCESS"
    days        = 90
  }

  tiering {
    access_tier = "DEEP_ARCHIVE_ACCESS"
    days        = 180
  }
}

# Outputs for S3 buckets
output "s3_raw_data_bucket_name" {
  description = "Name of the raw data S3 bucket"
  value       = aws_s3_bucket.raw_data.id
}

output "s3_raw_data_bucket_arn" {
  description = "ARN of the raw data S3 bucket"
  value       = aws_s3_bucket.raw_data.arn
}

output "s3_analysis_bucket_name" {
  description = "Name of the analysis S3 bucket"
  value       = aws_s3_bucket.analysis.id
}

output "s3_analysis_bucket_arn" {
  description = "ARN of the analysis S3 bucket"
  value       = aws_s3_bucket.analysis.arn
}

output "s3_scripts_bucket_name" {
  description = "Name of the scripts S3 bucket"
  value       = aws_s3_bucket.scripts.id
}

output "s3_scripts_bucket_arn" {
  description = "ARN of the scripts S3 bucket"
  value       = aws_s3_bucket.scripts.arn
}