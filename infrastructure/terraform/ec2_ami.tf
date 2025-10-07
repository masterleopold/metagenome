# MinION Metagenomics Pipeline - EC2 AMI Configuration
# AMI definitions for basecalling and analysis EC2 instances

# ===== Data source for Amazon Linux 2023 AMI (base) =====
data "aws_ami" "amazon_linux_2023" {
  most_recent = true
  owners      = ["amazon"]

  filter {
    name   = "name"
    values = ["al2023-ami-*-x86_64"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

  filter {
    name   = "root-device-type"
    values = ["ebs"]
  }

  filter {
    name   = "state"
    values = ["available"]
  }
}

# ===== Data source for Deep Learning AMI (for GPU/basecalling) =====
data "aws_ami" "deep_learning" {
  most_recent = true
  owners      = ["amazon"]

  filter {
    name   = "name"
    values = ["Deep Learning AMI GPU PyTorch * (Amazon Linux 2023)*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

  filter {
    name   = "root-device-type"
    values = ["ebs"]
  }

  filter {
    name   = "state"
    values = ["available"]
  }
}

# ===== Custom AMI for Basecalling (if already created) =====
data "aws_ami" "basecalling_custom" {
  count = var.basecalling_ami_id != "" ? 0 : 1

  most_recent = true
  owners      = ["self"]

  filter {
    name   = "name"
    values = ["${var.project_name}-basecalling-ami-*"]
  }

  filter {
    name   = "state"
    values = ["available"]
  }

  filter {
    name   = "tag:Environment"
    values = [var.environment]
  }
}

# ===== Custom AMI for Analysis (if already created) =====
data "aws_ami" "analysis_custom" {
  count = var.analysis_ami_id != "" ? 0 : 1

  most_recent = true
  owners      = ["self"]

  filter {
    name   = "name"
    values = ["${var.project_name}-analysis-ami-*"]
  }

  filter {
    name   = "state"
    values = ["available"]
  }

  filter {
    name   = "tag:Environment"
    values = [var.environment]
  }
}

# ===== Launch Template for Basecalling Instances =====
resource "aws_launch_template" "basecalling" {
  name_prefix   = "${var.project_name}-basecalling-"
  description   = "Launch template for MinION basecalling instances (GPU)"

  # Use custom AMI if available, otherwise use Deep Learning AMI
  image_id = var.basecalling_ami_id != "" ? var.basecalling_ami_id : (
    length(data.aws_ami.basecalling_custom) > 0 ?
    data.aws_ami.basecalling_custom[0].id :
    data.aws_ami.deep_learning.id
  )

  instance_type = var.basecalling_instance_type

  # Instance profile for IAM role
  iam_instance_profile {
    name = aws_iam_instance_profile.ec2_instance.name
  }

  # Key pair (optional, for SSH access)
  # key_name = var.key_pair_name

  # VPC security groups
  vpc_security_group_ids = [aws_security_group.ec2_instances.id]

  # Block device mappings
  block_device_mappings {
    device_name = "/dev/xvda"

    ebs {
      volume_size           = 200  # GB for OS and software
      volume_type           = "gp3"
      iops                  = 3000
      throughput            = 125
      delete_on_termination = true
      encrypted             = true
    }
  }

  # Additional storage for FAST5/FASTQ data
  block_device_mappings {
    device_name = "/dev/xvdf"

    ebs {
      volume_size           = 500  # GB for data
      volume_type           = "gp3"
      iops                  = 6000
      throughput            = 250
      delete_on_termination = true
      encrypted             = true
    }
  }

  # User data for instance initialization
  user_data = base64encode(templatefile("${path.module}/../../ec2_setup/userdata_basecalling.sh", {
    project_name   = var.project_name
    environment    = var.environment
    s3_scripts     = aws_s3_bucket.scripts.id
    efs_dns        = aws_efs_file_system.reference_db.dns_name
    dorado_ap_id   = aws_efs_access_point.dorado_models.id
  }))

  # Instance metadata options
  metadata_options {
    http_tokens                 = "required"  # IMDSv2
    http_put_response_hop_limit = 1
    http_endpoint               = "enabled"
  }

  # Monitoring
  monitoring {
    enabled = true
  }

  # Tags
  tag_specifications {
    resource_type = "instance"

    tags = {
      Name        = "${var.project_name}-basecalling"
      Environment = var.environment
      Purpose     = "MinION basecalling"
      Phase       = "1"
    }
  }

  tag_specifications {
    resource_type = "volume"

    tags = {
      Name        = "${var.project_name}-basecalling-volume"
      Environment = var.environment
    }
  }

  tags = {
    Name        = "${var.project_name}-basecalling-lt"
    Environment = var.environment
  }
}

# ===== Launch Template for Analysis Instances =====
resource "aws_launch_template" "analysis" {
  name_prefix   = "${var.project_name}-analysis-"
  description   = "Launch template for MinION analysis instances (CPU)"

  # Use custom AMI if available, otherwise use Amazon Linux 2023
  image_id = var.analysis_ami_id != "" ? var.analysis_ami_id : (
    length(data.aws_ami.analysis_custom) > 0 ?
    data.aws_ami.analysis_custom[0].id :
    data.aws_ami.amazon_linux_2023.id
  )

  # Note: instance_type will be overridden per phase
  instance_type = var.qc_instance_type

  # Instance profile for IAM role
  iam_instance_profile {
    name = aws_iam_instance_profile.ec2_instance.name
  }

  # VPC security groups
  vpc_security_group_ids = [aws_security_group.ec2_instances.id]

  # Block device mappings
  block_device_mappings {
    device_name = "/dev/xvda"

    ebs {
      volume_size           = 100  # GB for OS and software
      volume_type           = "gp3"
      iops                  = 3000
      throughput            = 125
      delete_on_termination = true
      encrypted             = true
    }
  }

  # Additional storage for analysis data
  block_device_mappings {
    device_name = "/dev/xvdf"

    ebs {
      volume_size           = 200  # GB for data
      volume_type           = "gp3"
      iops                  = 3000
      throughput            = 125
      delete_on_termination = true
      encrypted             = true
    }
  }

  # User data for instance initialization
  user_data = base64encode(templatefile("${path.module}/../../ec2_setup/userdata_analysis.sh", {
    project_name     = var.project_name
    environment      = var.environment
    s3_scripts       = aws_s3_bucket.scripts.id
    s3_analysis      = aws_s3_bucket.analysis.id
    efs_dns          = aws_efs_file_system.reference_db.dns_name
    kraken2_ap_id    = aws_efs_access_point.kraken2_db.id
    blast_ap_id      = aws_efs_access_point.blast_db.id
    host_genome_ap_id = aws_efs_access_point.host_genome.id
    perv_ap_id       = aws_efs_access_point.perv_db.id
  }))

  # Instance metadata options
  metadata_options {
    http_tokens                 = "required"  # IMDSv2
    http_put_response_hop_limit = 1
    http_endpoint               = "enabled"
  }

  # Monitoring
  monitoring {
    enabled = true
  }

  # Tags
  tag_specifications {
    resource_type = "instance"

    tags = {
      Name        = "${var.project_name}-analysis"
      Environment = var.environment
      Purpose     = "MinION analysis"
    }
  }

  tag_specifications {
    resource_type = "volume"

    tags = {
      Name        = "${var.project_name}-analysis-volume"
      Environment = var.environment
    }
  }

  tags = {
    Name        = "${var.project_name}-analysis-lt"
    Environment = var.environment
  }
}

# ===== Spot Fleet Request Configuration (Optional for cost optimization) =====
resource "aws_spot_fleet_request" "analysis_spot" {
  count = var.environment == "dev" ? 1 : 0  # Only use spot in dev environment

  iam_fleet_role                      = aws_iam_role.spot_fleet[0].arn
  allocation_strategy                  = "lowestPrice"
  target_capacity                      = 0  # Start with 0, scale as needed
  valid_until                         = timeadd(timestamp(), "8760h")  # 1 year
  terminate_instances_with_expiration = true
  instance_interruption_behaviour      = "terminate"

  launch_template_config {
    launch_template_specification {
      id      = aws_launch_template.analysis.id
      version = "$Latest"
    }

    overrides {
      instance_type = var.qc_instance_type
      spot_price    = "0.10"  # Max price per hour
    }

    overrides {
      instance_type = var.host_removal_instance_type
      spot_price    = "0.30"
    }

    overrides {
      instance_type = var.blast_instance_type
      spot_price    = "0.50"
    }
  }

  tags = {
    Name        = "${var.project_name}-spot-fleet"
    Environment = var.environment
  }
}

# IAM Role for Spot Fleet (if using)
resource "aws_iam_role" "spot_fleet" {
  count = var.environment == "dev" ? 1 : 0

  name = "${var.project_name}-spot-fleet-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "spotfleet.amazonaws.com"
        }
      }
    ]
  })

  tags = {
    Name        = "${var.project_name}-spot-fleet-role"
    Environment = var.environment
  }
}

resource "aws_iam_role_policy_attachment" "spot_fleet" {
  count = var.environment == "dev" ? 1 : 0

  role       = aws_iam_role.spot_fleet[0].name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
}

# ===== EC2 Instance Connect Endpoint (for secure access) =====
resource "aws_ec2_instance_connect_endpoint" "main" {
  subnet_id          = module.vpc.private_subnets[0]
  security_group_ids = [aws_security_group.ec2_instances.id]

  tags = {
    Name        = "${var.project_name}-eice"
    Environment = var.environment
  }
}

# ===== Outputs =====
output "basecalling_ami_id" {
  description = "AMI ID used for basecalling instances"
  value = var.basecalling_ami_id != "" ? var.basecalling_ami_id : (
    length(data.aws_ami.basecalling_custom) > 0 ?
    data.aws_ami.basecalling_custom[0].id :
    data.aws_ami.deep_learning.id
  )
}

output "analysis_ami_id" {
  description = "AMI ID used for analysis instances"
  value = var.analysis_ami_id != "" ? var.analysis_ami_id : (
    length(data.aws_ami.analysis_custom) > 0 ?
    data.aws_ami.analysis_custom[0].id :
    data.aws_ami.amazon_linux_2023.id
  )
}

output "basecalling_launch_template_id" {
  description = "Launch template ID for basecalling instances"
  value       = aws_launch_template.basecalling.id
}

output "analysis_launch_template_id" {
  description = "Launch template ID for analysis instances"
  value       = aws_launch_template.analysis.id
}

output "ec2_instance_connect_endpoint_id" {
  description = "EC2 Instance Connect Endpoint ID"
  value       = aws_ec2_instance_connect_endpoint.main.id
}