# MinION Metagenomics Pipeline - Lambda Functions Configuration
# Lambda functions for pipeline orchestration and phase triggers

# ===== Lambda Function: Orchestrator =====
resource "aws_lambda_function" "orchestrator" {
  function_name = "${var.project_name}-orchestrator"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_orchestrator.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  # Deployment package (to be created)
  filename         = "${path.module}/../../lambda/orchestrator/deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/orchestrator/deployment.zip")

  environment {
    variables = {
      PROJECT_NAME          = var.project_name
      ENVIRONMENT           = var.environment
      S3_RAW_BUCKET         = aws_s3_bucket.raw_data.id
      S3_ANALYSIS_BUCKET    = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET     = aws_s3_bucket.scripts.id
      RDS_ENDPOINT          = aws_db_instance.metadata.endpoint
      RDS_DATABASE          = aws_db_instance.metadata.db_name
      RDS_SECRET_ARN        = aws_secretsmanager_secret.rds_password.arn
      SNS_TOPIC_ARN         = aws_sns_topic.alerts.arn
      EC2_INSTANCE_PROFILE  = aws_iam_instance_profile.ec2_instance.name
      BASECALLING_AMI       = var.basecalling_ami_id
      BASECALLING_INSTANCE  = var.basecalling_instance_type
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-orchestrator"
    Environment = var.environment
    Purpose     = "Main workflow orchestrator"
  }
}

# ===== Lambda Functions: Phase Triggers =====

# Phase 1 Trigger - Basecalling
resource "aws_lambda_function" "phase1_trigger" {
  function_name = "${var.project_name}-phase1-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase1_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  filename         = "${path.module}/../../lambda/phase_triggers/phase1_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase1_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME         = var.project_name
      S3_ANALYSIS_BUCKET   = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET    = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE = aws_iam_instance_profile.ec2_instance.name
      BASECALLING_AMI      = var.basecalling_ami_id
      INSTANCE_TYPE        = var.basecalling_instance_type
      SNS_TOPIC_ARN        = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase1-trigger"
    Environment = var.environment
    Phase       = "1-Basecalling"
  }
}

# Phase 2 Trigger - QC
resource "aws_lambda_function" "phase2_trigger" {
  function_name = "${var.project_name}-phase2-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase2_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  filename         = "${path.module}/../../lambda/phase_triggers/phase2_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase2_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME         = var.project_name
      S3_ANALYSIS_BUCKET   = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET    = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE = aws_iam_instance_profile.ec2_instance.name
      ANALYSIS_AMI         = var.analysis_ami_id
      INSTANCE_TYPE        = var.qc_instance_type
      SNS_TOPIC_ARN        = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase2-trigger"
    Environment = var.environment
    Phase       = "2-QC"
  }
}

# Phase 3 Trigger - Host Removal
resource "aws_lambda_function" "phase3_trigger" {
  function_name = "${var.project_name}-phase3-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase3_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  filename         = "${path.module}/../../lambda/phase_triggers/phase3_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase3_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME         = var.project_name
      S3_ANALYSIS_BUCKET   = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET    = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE = aws_iam_instance_profile.ec2_instance.name
      ANALYSIS_AMI         = var.analysis_ami_id
      INSTANCE_TYPE        = var.host_removal_instance_type
      EFS_HOST_GENOME_AP   = aws_efs_access_point.host_genome.id
      SNS_TOPIC_ARN        = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase3-trigger"
    Environment = var.environment
    Phase       = "3-HostRemoval"
  }
}

# Phase 4 Trigger - Pathogen Detection (launches 4 parallel EC2)
resource "aws_lambda_function" "phase4_trigger" {
  function_name = "${var.project_name}-phase4-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase4_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = 1024  # Increased memory for parallel EC2 launches

  filename         = "${path.module}/../../lambda/phase_triggers/phase4_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase4_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME           = var.project_name
      S3_ANALYSIS_BUCKET     = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET      = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE   = aws_iam_instance_profile.ec2_instance.name
      ANALYSIS_AMI           = var.analysis_ami_id
      KRAKEN2_INSTANCE_TYPE  = var.kraken2_instance_type
      BLAST_INSTANCE_TYPE    = var.blast_instance_type
      ASSEMBLY_INSTANCE_TYPE = var.assembly_instance_type
      PERV_INSTANCE_TYPE     = var.perv_instance_type
      EFS_KRAKEN2_AP         = aws_efs_access_point.kraken2_db.id
      EFS_BLAST_AP           = aws_efs_access_point.blast_db.id
      EFS_PERV_AP            = aws_efs_access_point.perv_db.id
      SNS_TOPIC_ARN          = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase4-trigger"
    Environment = var.environment
    Phase       = "4-PathogenDetection"
  }
}

# Phase 5 Trigger - Quantification
resource "aws_lambda_function" "phase5_trigger" {
  function_name = "${var.project_name}-phase5-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase5_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  filename         = "${path.module}/../../lambda/phase_triggers/phase5_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase5_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME         = var.project_name
      S3_ANALYSIS_BUCKET   = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET    = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE = aws_iam_instance_profile.ec2_instance.name
      ANALYSIS_AMI         = var.analysis_ami_id
      INSTANCE_TYPE        = var.quantification_instance_type
      RDS_ENDPOINT         = aws_db_instance.metadata.endpoint
      RDS_DATABASE         = aws_db_instance.metadata.db_name
      RDS_SECRET_ARN       = aws_secretsmanager_secret.rds_password.arn
      SNS_TOPIC_ARN        = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase5-trigger"
    Environment = var.environment
    Phase       = "5-Quantification"
  }
}

# Phase 6 Trigger - Report Generation
resource "aws_lambda_function" "phase6_trigger" {
  function_name = "${var.project_name}-phase6-trigger"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_phase6_trigger.handler"
  runtime       = var.lambda_runtime
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  filename         = "${path.module}/../../lambda/phase_triggers/phase6_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/phase_triggers/phase6_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME         = var.project_name
      S3_ANALYSIS_BUCKET   = aws_s3_bucket.analysis.id
      S3_SCRIPTS_BUCKET    = aws_s3_bucket.scripts.id
      EC2_INSTANCE_PROFILE = aws_iam_instance_profile.ec2_instance.name
      ANALYSIS_AMI         = var.analysis_ami_id
      INSTANCE_TYPE        = var.reports_instance_type
      RDS_ENDPOINT         = aws_db_instance.metadata.endpoint
      RDS_DATABASE         = aws_db_instance.metadata.db_name
      RDS_SECRET_ARN       = aws_secretsmanager_secret.rds_password.arn
      SNS_TOPIC_ARN        = aws_sns_topic.alerts.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-phase6-trigger"
    Environment = var.environment
    Phase       = "6-ReportGeneration"
  }
}

# ===== Validator Lambda Functions =====

# Input Validator
resource "aws_lambda_function" "validate_input" {
  function_name = "${var.project_name}-validate-input"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_validate_input.handler"
  runtime       = var.lambda_runtime
  timeout       = 60
  memory_size   = 256

  filename         = "${path.module}/../../lambda/validators/validate_input_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/validators/validate_input_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME       = var.project_name
      S3_RAW_BUCKET      = aws_s3_bucket.raw_data.id
      S3_ANALYSIS_BUCKET = aws_s3_bucket.analysis.id
    }
  }

  tags = {
    Name        = "${var.project_name}-validate-input"
    Environment = var.environment
    Purpose     = "Validate FAST5 input files"
  }
}

# QC Checker
resource "aws_lambda_function" "check_qc" {
  function_name = "${var.project_name}-check-qc"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_check_qc.handler"
  runtime       = var.lambda_runtime
  timeout       = 60
  memory_size   = 256

  filename         = "${path.module}/../../lambda/validators/check_qc_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/validators/check_qc_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME       = var.project_name
      S3_ANALYSIS_BUCKET = aws_s3_bucket.analysis.id
      MIN_TOTAL_BASES    = "25000000000"  # 25GB
      MIN_MEAN_QSCORE    = "18"
      MIN_N50            = "8000"
      MIN_READS          = "1000000"
    }
  }

  tags = {
    Name        = "${var.project_name}-check-qc"
    Environment = var.environment
    Purpose     = "Check basecalling QC metrics"
  }
}

# ===== Metadata Lambda Functions =====

# Metadata Updater
resource "aws_lambda_function" "update_metadata" {
  function_name = "${var.project_name}-update-metadata"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_update_metadata.handler"
  runtime       = var.lambda_runtime
  timeout       = 60
  memory_size   = 256

  filename         = "${path.module}/../../lambda/metadata/update_metadata_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/metadata/update_metadata_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME   = var.project_name
      RDS_ENDPOINT   = aws_db_instance.metadata.endpoint
      RDS_DATABASE   = aws_db_instance.metadata.db_name
      RDS_SECRET_ARN = aws_secretsmanager_secret.rds_password.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-update-metadata"
    Environment = var.environment
    Purpose     = "Update RDS metadata"
  }
}

# Notifier
resource "aws_lambda_function" "notify" {
  function_name = "${var.project_name}-notify"
  role          = aws_iam_role.lambda_execution.arn
  handler       = "lambda_notify.handler"
  runtime       = var.lambda_runtime
  timeout       = 60
  memory_size   = 256

  filename         = "${path.module}/../../lambda/metadata/notify_deployment.zip"
  source_code_hash = filebase64sha256("${path.module}/../../lambda/metadata/notify_deployment.zip")

  environment {
    variables = {
      PROJECT_NAME   = var.project_name
      SNS_TOPIC_ARN  = aws_sns_topic.alerts.arn
      RDS_ENDPOINT   = aws_db_instance.metadata.endpoint
      RDS_DATABASE   = aws_db_instance.metadata.db_name
      RDS_SECRET_ARN = aws_secretsmanager_secret.rds_password.arn
    }
  }

  vpc_config {
    subnet_ids         = module.vpc.private_subnets
    security_group_ids = [aws_security_group.lambda.id]
  }

  tags = {
    Name        = "${var.project_name}-notify"
    Environment = var.environment
    Purpose     = "Send notifications via SNS"
  }
}

# ===== Lambda Aliases for Production =====
resource "aws_lambda_alias" "orchestrator_prod" {
  name             = "PROD"
  description      = "Production alias for orchestrator"
  function_name    = aws_lambda_function.orchestrator.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase1_prod" {
  name             = "PROD"
  description      = "Production alias for phase1 trigger"
  function_name    = aws_lambda_function.phase1_trigger.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase2_prod" {
  name             = "PROD"
  description      = "Production alias for phase2 trigger"
  function_name    = aws_lambda_function.phase2_trigger.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase3_prod" {
  name             = "PROD"
  description      = "Production alias for phase3 trigger"
  function_name    = aws_lambda_function.phase3_trigger.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase4_prod" {
  name             = "PROD"
  description      = "Production alias for phase4 trigger"
  function_name    = aws_lambda_function.phase4_trigger.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase5_prod" {
  name             = "PROD"
  description      = "Production alias for phase5 trigger"
  function_name    = aws_lambda_function.phase5_trigger.function_name
  function_version = "$LATEST"
}

resource "aws_lambda_alias" "phase6_prod" {
  name             = "PROD"
  description      = "Production alias for phase6 trigger"
  function_name    = aws_lambda_function.phase6_trigger.function_name
  function_version = "$LATEST"
}

# ===== CloudWatch Log Groups for Lambda =====
resource "aws_cloudwatch_log_group" "lambda_orchestrator" {
  name              = "/aws/lambda/${aws_lambda_function.orchestrator.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-orchestrator-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase1" {
  name              = "/aws/lambda/${aws_lambda_function.phase1_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase1-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase2" {
  name              = "/aws/lambda/${aws_lambda_function.phase2_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase2-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase3" {
  name              = "/aws/lambda/${aws_lambda_function.phase3_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase3-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase4" {
  name              = "/aws/lambda/${aws_lambda_function.phase4_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase4-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase5" {
  name              = "/aws/lambda/${aws_lambda_function.phase5_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase5-logs"
    Environment = var.environment
  }
}

resource "aws_cloudwatch_log_group" "lambda_phase6" {
  name              = "/aws/lambda/${aws_lambda_function.phase6_trigger.function_name}"
  retention_in_days = 30

  tags = {
    Name        = "${var.project_name}-phase6-logs"
    Environment = var.environment
  }
}

# Outputs for Lambda functions
output "lambda_orchestrator_arn" {
  description = "ARN of the orchestrator Lambda function"
  value       = aws_lambda_function.orchestrator.arn
}

output "lambda_orchestrator_name" {
  description = "Name of the orchestrator Lambda function"
  value       = aws_lambda_function.orchestrator.function_name
}

output "lambda_phase_triggers" {
  description = "Map of phase trigger Lambda functions"
  value = {
    phase1 = aws_lambda_function.phase1_trigger.function_name
    phase2 = aws_lambda_function.phase2_trigger.function_name
    phase3 = aws_lambda_function.phase3_trigger.function_name
    phase4 = aws_lambda_function.phase4_trigger.function_name
    phase5 = aws_lambda_function.phase5_trigger.function_name
    phase6 = aws_lambda_function.phase6_trigger.function_name
  }
}