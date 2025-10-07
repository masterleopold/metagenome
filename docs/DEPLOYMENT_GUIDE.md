# MinION Pipeline Deployment Guide

## Prerequisites

### AWS Account Setup

1. **AWS Account**: Active AWS account with appropriate permissions
2. **IAM User**: Administrator or PowerUser access
3. **AWS CLI**: Version 2.0+ configured with credentials
4. **Service Quotas**: Ensure adequate quotas for:
   - EC2 instances (especially GPU instances)
   - Lambda concurrent executions
   - S3 buckets
   - EFS file systems

### Local Environment

```bash
# Required tools
aws --version        # AWS CLI 2.0+
terraform --version  # Terraform 1.0+
python3 --version    # Python 3.9+
docker --version     # Docker 20.10+ (optional)

# Python dependencies
pip3 install -r requirements.txt
```

## Step-by-Step Deployment

### Step 1: Environment Configuration

```bash
# Set environment variables
export AWS_REGION=ap-northeast-1
export ENVIRONMENT=production
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Create deployment configuration
cat > deployment.env << EOF
AWS_REGION=$AWS_REGION
ENVIRONMENT=$ENVIRONMENT
AWS_ACCOUNT_ID=$AWS_ACCOUNT_ID
ALERT_EMAIL=admin@your-domain.com
DOMAIN_NAME=api.your-domain.com  # Optional
EOF

source deployment.env
```

### Step 2: Infrastructure Deployment

#### 2.1 Initialize Terraform

```bash
cd infrastructure/terraform

# Initialize backend
terraform init \
  -backend-config="bucket=terraform-state-$AWS_ACCOUNT_ID" \
  -backend-config="key=minion-pipeline/$ENVIRONMENT/terraform.tfstate" \
  -backend-config="region=$AWS_REGION"
```

#### 2.2 Plan Infrastructure

```bash
terraform plan \
  -var="environment=$ENVIRONMENT" \
  -var="region=$AWS_REGION" \
  -out=tfplan

# Review the plan
terraform show tfplan
```

#### 2.3 Apply Infrastructure

```bash
terraform apply tfplan

# Save outputs
terraform output -json > ../outputs.json
```

### Step 3: Build AMIs

#### 3.1 Basecalling AMI (GPU)

```bash
cd ec2_setup

./build_basecalling_ami.sh

# Note the AMI ID
export BASECALLING_AMI_ID=ami-xxxxxxxxx
```

#### 3.2 Analysis AMI

```bash
./build_analysis_ami.sh

# Note the AMI ID
export ANALYSIS_AMI_ID=ami-xxxxxxxxx
```

### Step 4: Database Setup

#### 4.1 Mount EFS

Launch a temporary EC2 instance to set up databases:

```bash
# Launch setup instance
aws ec2 run-instances \
  --image-id $ANALYSIS_AMI_ID \
  --instance-type t3.xlarge \
  --subnet-id $(terraform output -raw private_subnet_id) \
  --security-group-ids $(terraform output -raw security_group_id) \
  --iam-instance-profile Name=MinIONEC2Role \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=database-setup}]'

# Get instance ID
INSTANCE_ID=i-xxxxxxxxx

# Connect via SSM
aws ssm start-session --target $INSTANCE_ID
```

#### 4.2 Install Databases

On the EC2 instance:

```bash
# Mount EFS
sudo mkdir -p /mnt/efs
sudo mount -t nfs4 $(terraform output -raw efs_dns)/ /mnt/efs

# Run database setup
cd /opt/minion
./tools/database_setup.sh --all

# Verify installation
./tools/database_setup.sh --check
```

### Step 5: Deploy Lambda Functions

#### 5.1 Package Functions

```bash
cd lambda

./package_functions.sh

# Creates ZIP files for each function
ls -la *.zip
```

#### 5.2 Deploy Functions

```bash
# Deploy orchestration functions
aws lambda create-function \
  --function-name minion-pipeline-orchestrator-$ENVIRONMENT \
  --runtime python3.9 \
  --role arn:aws:iam::$AWS_ACCOUNT_ID:role/MinIONLambdaRole \
  --handler pipeline_orchestrator.lambda_handler \
  --zip-file fileb://orchestration.zip \
  --timeout 30 \
  --memory-size 256 \
  --environment Variables="{
    ENVIRONMENT=$ENVIRONMENT,
    STATE_MACHINE_ARN=arn:aws:states:$AWS_REGION:$AWS_ACCOUNT_ID:stateMachine:minion-pipeline-$ENVIRONMENT,
    SNS_TOPIC_ARN=$(terraform output -raw sns_topic_arn)
  }"

# Repeat for other functions...
```

### Step 6: Deploy Step Functions

```bash
cd templates/stepfunctions

# Replace variables
sed -e "s/\${AWS_REGION}/$AWS_REGION/g" \
    -e "s/\${AWS_ACCOUNT}/$AWS_ACCOUNT_ID/g" \
    workflow_definition.json > /tmp/workflow.json

# Create state machine
aws stepfunctions create-state-machine \
  --name minion-pipeline-$ENVIRONMENT \
  --definition file:///tmp/workflow.json \
  --role-arn arn:aws:iam::$AWS_ACCOUNT_ID:role/StepFunctionsRole-$ENVIRONMENT
```

### Step 7: Setup Monitoring

```bash
cd templates/cloudformation

# Deploy monitoring stack
aws cloudformation deploy \
  --template-file monitoring_stack.yaml \
  --stack-name minion-monitoring-$ENVIRONMENT \
  --parameter-overrides \
    Environment=$ENVIRONMENT \
    AlertEmail=$ALERT_EMAIL \
  --capabilities CAPABILITY_IAM
```

### Step 8: Deploy API Gateway

```bash
# Deploy API stack
aws cloudformation deploy \
  --template-file api_gateway.yaml \
  --stack-name minion-api-$ENVIRONMENT \
  --parameter-overrides \
    Environment=$ENVIRONMENT \
    DomainName=$DOMAIN_NAME \
  --capabilities CAPABILITY_IAM CAPABILITY_AUTO_EXPAND

# Get API endpoint
API_ENDPOINT=$(aws cloudformation describe-stacks \
  --stack-name minion-api-$ENVIRONMENT \
  --query 'Stacks[0].Outputs[?OutputKey==`ApiEndpoint`].OutputValue' \
  --output text)

echo "API Endpoint: $API_ENDPOINT"
```

### Step 9: Configure S3 Event Triggers

```bash
# Create event notification configuration
cat > s3-notification.json << EOF
{
  "LambdaFunctionConfigurations": [
    {
      "LambdaFunctionArn": "arn:aws:lambda:$AWS_REGION:$AWS_ACCOUNT_ID:function:minion-pipeline-orchestrator-$ENVIRONMENT",
      "Events": ["s3:ObjectCreated:*"],
      "Filter": {
        "Key": {
          "FilterRules": [
            {
              "Name": "prefix",
              "Value": "runs/"
            },
            {
              "Name": "suffix",
              "Value": ".fast5"
            }
          ]
        }
      }
    }
  ]
}
EOF

# Apply to S3 bucket
aws s3api put-bucket-notification-configuration \
  --bucket minion-data-$ENVIRONMENT \
  --notification-configuration file://s3-notification.json
```

### Step 10: Validation

```bash
# Run deployment validation
./tools/deployment_script.sh validate

# Check all services
./tools/deployment_script.sh status
```

## Post-Deployment Configuration

### 1. API Keys

Generate and secure API keys:

```bash
# Generate API key
API_KEY=$(openssl rand -hex 32)

# Store in Secrets Manager
aws secretsmanager create-secret \
  --name minion-api-key-$ENVIRONMENT \
  --secret-string "$API_KEY"
```

### 2. DNS Configuration

If using custom domain:

```bash
# Get API Gateway domain
aws apigateway get-domain-names

# Create Route53 record
aws route53 change-resource-record-sets \
  --hosted-zone-id ZXXXXXXXXXXXXX \
  --change-batch file://route53-change.json
```

### 3. Cost Alerts

Set up billing alerts:

```bash
aws cloudwatch put-metric-alarm \
  --alarm-name minion-cost-alarm-$ENVIRONMENT \
  --alarm-description "MinION pipeline cost alert" \
  --metric-name EstimatedCharges \
  --namespace AWS/Billing \
  --statistic Maximum \
  --period 86400 \
  --evaluation-periods 1 \
  --threshold 500 \
  --comparison-operator GreaterThanThreshold
```

## Verification Tests

### 1. Test Workflow Execution

```bash
# Create test data
echo "test" > test.fast5
aws s3 cp test.fast5 s3://minion-data-$ENVIRONMENT/runs/TEST-001/fast5/

# Start test workflow
./tools/workflow_cli.py start \
  --run-id TEST-001 \
  --bucket minion-data-$ENVIRONMENT \
  --input-prefix runs/TEST-001/fast5/

# Monitor progress
./tools/workflow_cli.py status --run-id TEST-001 --watch
```

### 2. Test API Endpoints

```bash
# Test API
curl -X GET $API_ENDPOINT/workflows \
  -H "x-api-key: $API_KEY"
```

### 3. Test Monitoring

```bash
# Check CloudWatch dashboard
echo "Dashboard URL: https://console.aws.amazon.com/cloudwatch/home?region=$AWS_REGION#dashboards:name=minion-pipeline-$ENVIRONMENT"

# Test SNS notifications
aws sns publish \
  --topic-arn $(terraform output -raw sns_topic_arn) \
  --subject "Test Alert" \
  --message "Testing SNS notifications"
```

## Rollback Procedure

If deployment fails:

### 1. Terraform Rollback

```bash
cd infrastructure/terraform

# Revert to previous state
terraform state pull > current_state.json
terraform state push previous_state.json
```

### 2. CloudFormation Rollback

```bash
# Rollback stack
aws cloudformation cancel-update-stack \
  --stack-name minion-monitoring-$ENVIRONMENT

# Or delete and recreate
aws cloudformation delete-stack \
  --stack-name minion-monitoring-$ENVIRONMENT
```

### 3. Lambda Rollback

```bash
# Revert to previous version
aws lambda update-function-code \
  --function-name minion-pipeline-orchestrator-$ENVIRONMENT \
  --s3-bucket lambda-deployments-$AWS_ACCOUNT_ID \
  --s3-key previous/orchestration.zip
```

## Troubleshooting

### Common Issues

#### 1. Terraform State Lock

```bash
# Force unlock
terraform force-unlock LOCK_ID
```

#### 2. Lambda Timeout

```bash
# Increase timeout
aws lambda update-function-configuration \
  --function-name FUNCTION_NAME \
  --timeout 300
```

#### 3. EFS Mount Failed

```bash
# Check security group rules
aws ec2 describe-security-groups --group-ids sg-xxxxxxxxx

# Add NFS rule
aws ec2 authorize-security-group-ingress \
  --group-id sg-xxxxxxxxx \
  --protocol tcp \
  --port 2049 \
  --source-group sg-xxxxxxxxx
```

#### 4. Step Functions Failed

```bash
# Check execution history
aws stepfunctions get-execution-history \
  --execution-arn arn:aws:states:... \
  --max-results 10
```

## Maintenance

### Daily Tasks

- Check CloudWatch dashboard
- Review error logs
- Monitor costs

### Weekly Tasks

- Review pathogen detection reports
- Check database updates
- Validate backups

### Monthly Tasks

- Update AMIs with security patches
- Review and optimize costs
- Update reference databases

## Support

For deployment support:
- Documentation: `/docs`
- Issues: GitHub Issues
- Email: devops@your-org.com