import { Sidebar } from "@/components/layout/Sidebar";
import { Badge } from "@/components/ui/Badge";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/Card";
import { CodeBlock } from "@/components/ui/CodeBlock";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/Alert";
import {
  CheckCircleIcon,
  AlertTriangleIcon,
  CloudIcon,
  ServerIcon,
  DatabaseIcon,
  SettingsIcon,
} from "lucide-react";

export default function DeploymentPage() {
  return (
    <div className="flex-1">
      <div className="container mx-auto flex gap-6 px-4">
        <Sidebar />
        <main className="flex-1 py-6 px-4 md:px-6 lg:px-8">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <Badge variant="secondary" className="mb-4">Development</Badge>
              <h1 className="text-4xl font-bold mb-4">Deployment Guide</h1>
              <p className="text-lg text-muted-foreground">
                Step-by-step guide to deploy the MinION pipeline to AWS.
              </p>
            </div>

            <Alert variant="warning" className="mb-8">
              <AlertTriangleIcon className="h-4 w-4" />
              <AlertTitle>Deployment Time & Cost</AlertTitle>
              <AlertDescription>
                Full deployment takes approximately 30-45 minutes. AWS infrastructure will incur charges (estimated $100-300/month for development environment).
              </AlertDescription>
            </Alert>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Prerequisites</h2>

              <div className="grid md:grid-cols-2 gap-6">
                <Card>
                  <CardHeader>
                    <CloudIcon className="h-8 w-8 text-primary mb-2" />
                    <CardTitle>AWS Account</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Active AWS account with admin/PowerUser access</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>AWS CLI 2.0+ installed and configured</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Region: ap-northeast-1 (Tokyo) recommended</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Service quotas verified for GPU instances</span>
                      </li>
                    </ul>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <ServerIcon className="h-8 w-8 text-secondary mb-2" />
                    <CardTitle>Local Environment</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <ul className="text-sm space-y-2">
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Terraform 1.0+ installed</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Python 3.9+ with pip</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Git for repository access</span>
                      </li>
                      <li className="flex items-start gap-2">
                        <CheckCircleIcon className="h-5 w-5 text-primary mt-0.5 shrink-0" />
                        <span>Docker (optional, for local testing)</span>
                      </li>
                    </ul>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 1: Environment Configuration</h2>

              <Card>
                <CardHeader>
                  <CardTitle>Set Environment Variables</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <CodeBlock
                    code={`# Set AWS region and environment
export AWS_REGION=ap-northeast-1
export ENVIRONMENT=production
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Create deployment configuration file
cat > deployment.env << EOF
AWS_REGION=$AWS_REGION
ENVIRONMENT=$ENVIRONMENT
AWS_ACCOUNT_ID=$AWS_ACCOUNT_ID
ALERT_EMAIL=admin@your-domain.com
DOMAIN_NAME=api.your-domain.com  # Optional
EOF

# Load configuration
source deployment.env`}
                    language="bash"
                  />
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 2: Deploy Infrastructure with Terraform</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Initialize Terraform</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`cd infrastructure/terraform

# Initialize Terraform backend
terraform init \\
  -backend-config="bucket=terraform-state-$AWS_ACCOUNT_ID" \\
  -backend-config="key=minion-pipeline/$ENVIRONMENT/terraform.tfstate" \\
  -backend-config="region=$AWS_REGION"`}
                      language="bash"
                    />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Plan Infrastructure</CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    <CodeBlock
                      code={`# Generate execution plan
terraform plan \\
  -var="environment=$ENVIRONMENT" \\
  -var="region=$AWS_REGION" \\
  -out=tfplan

# Review the plan
terraform show tfplan`}
                      language="bash"
                    />
                    <p className="text-sm text-muted-foreground">
                      Review the planned changes carefully. This will create ~30 AWS resources including VPC, S3, RDS, Lambda, etc.
                    </p>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Apply Infrastructure</CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    <CodeBlock
                      code={`# Apply the infrastructure
terraform apply tfplan

# Save outputs for later use
terraform output -json > outputs.json

# View key outputs
terraform output`}
                      language="bash"
                    />
                    <Alert>
                      <AlertDescription>
                        This step takes 10-15 minutes to complete. Terraform will create all AWS resources defined in the configuration.
                      </AlertDescription>
                    </Alert>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 3: Build AMIs</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Basecalling AMI (GPU-enabled)</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`cd ec2_setup

# Build GPU-optimized AMI with Dorado
./build_basecalling_ami.sh

# Note the AMI ID from output
export BASECALLING_AMI_ID=ami-xxxxxxxxx`}
                      language="bash"
                    />
                    <p className="text-sm text-muted-foreground mt-3">
                      This AMI includes CUDA drivers, NVIDIA Docker, and Dorado basecaller. Build time: ~20 minutes.
                    </p>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Analysis AMI</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`# Build general analysis AMI
./build_analysis_ami.sh

# Note the AMI ID from output
export ANALYSIS_AMI_ID=ami-yyyyyyyyy`}
                      language="bash"
                    />
                    <p className="text-sm text-muted-foreground mt-3">
                      This AMI includes Kraken2, BLAST, Minimap2, SAMtools, and all Python dependencies. Build time: ~15 minutes.
                    </p>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 4: Setup Reference Databases</h2>

              <Card>
                <CardHeader>
                  <DatabaseIcon className="h-8 w-8 text-primary mb-2" />
                  <CardTitle>Install Databases on EFS</CardTitle>
                  <CardDescription>Kraken2, RVDB, BLAST, PMDA databases</CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <Alert variant="warning">
                    <AlertTriangleIcon className="h-4 w-4" />
                    <AlertDescription>
                      Database download is ~150GB and takes 2-4 hours. Ensure stable internet connection.
                    </AlertDescription>
                  </Alert>

                  <CodeBlock
                    code={`# Launch temporary EC2 instance for database setup
aws ec2 run-instances \\
  --image-id $ANALYSIS_AMI_ID \\
  --instance-type t3.xlarge \\
  --subnet-id $(terraform output -raw private_subnet_id) \\
  --security-group-ids $(terraform output -raw security_group_id) \\
  --iam-instance-profile Name=MinIONEC2Role \\
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=database-setup}]'

# Get instance ID and connect via SSM
INSTANCE_ID=i-xxxxxxxxx
aws ssm start-session --target $INSTANCE_ID

# On the EC2 instance:
sudo mkdir -p /mnt/efs
sudo mount -t nfs4 $(terraform output -raw efs_dns)/ /mnt/efs

cd /opt/minion
./tools/database_setup.sh --all

# Verify installation
./tools/database_setup.sh --check`}
                    language="bash"
                  />

                  <div>
                    <h4 className="font-semibold mb-2">Databases Installed</h4>
                    <ul className="text-sm space-y-1">
                      <li>✓ Kraken2 Standard Database (~50GB)</li>
                      <li>✓ RVDB v30.0 Viral Database (~10GB)</li>
                      <li>✓ BLAST PMDA Pathogen Database (~5GB)</li>
                      <li>✓ Sus scrofa Reference Genome (~3GB)</li>
                      <li>✓ PERV Reference Sequences (~1MB)</li>
                    </ul>
                  </div>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 5: Deploy Lambda Functions</h2>

              <Card>
                <CardHeader>
                  <CardTitle>Package and Deploy Functions</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <CodeBlock
                    code={`cd lambda

# Package all Lambda functions
./package_functions.sh

# Deploy orchestration functions
aws lambda create-function \\
  --function-name minion-pipeline-orchestrator-$ENVIRONMENT \\
  --runtime python3.9 \\
  --role arn:aws:iam::$AWS_ACCOUNT_ID:role/MinIONLambdaRole \\
  --handler pipeline_orchestrator.lambda_handler \\
  --zip-file fileb://orchestration.zip \\
  --timeout 30 \\
  --memory-size 256 \\
  --environment Variables="{\\"ENVIRONMENT\\":\\"$ENVIRONMENT\\",\\"STATE_MACHINE_ARN\\":\\"arn:aws:states:$AWS_REGION:$AWS_ACCOUNT_ID:stateMachine:minion-pipeline-$ENVIRONMENT\\"}"

# Deploy remaining functions (repeat for each)
# - EC2 management functions (3)
# - Data processing functions (3)
# - Monitoring functions (4)
# - Reporting functions (3)`}
                    language="bash"
                  />
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 6: Configure S3 Event Triggers</h2>

              <Card>
                <CardHeader>
                  <CardTitle>S3 to Lambda Integration</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <CodeBlock
                    code={`# Create S3 event notification configuration
cat > s3-notification.json << EOF
{
  "LambdaFunctionConfigurations": [
    {
      "LambdaFunctionArn": "arn:aws:lambda:$AWS_REGION:$AWS_ACCOUNT_ID:function:minion-pipeline-orchestrator-$ENVIRONMENT",
      "Events": ["s3:ObjectCreated:*"],
      "Filter": {
        "Key": {
          "FilterRules": [
            {"Name": "prefix", "Value": "runs/"},
            {"Name": "suffix", "Value": ".fast5"}
          ]
        }
      }
    }
  ]
}
EOF

# Apply to S3 bucket
aws s3api put-bucket-notification-configuration \\
  --bucket minion-data-$ENVIRONMENT \\
  --notification-configuration file://s3-notification.json`}
                    language="bash"
                  />
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Step 7: Validation</h2>

              <Card>
                <CardHeader>
                  <CardTitle>Verify Deployment</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <CodeBlock
                    code={`# Run deployment validation
./tools/deployment_script.sh validate

# Check all services status
./tools/deployment_script.sh status

# Test with sample data
echo "test" > test.fast5
aws s3 cp test.fast5 s3://minion-data-$ENVIRONMENT/runs/TEST-001/fast5/

# Start test workflow
./tools/workflow_cli.py start \\
  --run-id TEST-001 \\
  --bucket minion-data-$ENVIRONMENT \\
  --input-prefix runs/TEST-001/fast5/

# Monitor progress
./tools/workflow_cli.py status --run-id TEST-001 --watch`}
                    language="bash"
                  />

                  <Alert variant="success">
                    <CheckCircleIcon className="h-4 w-4" />
                    <AlertTitle>Deployment Complete!</AlertTitle>
                    <AlertDescription>
                      If all validation checks pass, your MinION pipeline is ready for production use.
                    </AlertDescription>
                  </Alert>
                </CardContent>
              </Card>
            </section>

            <section className="mb-12">
              <h2 className="text-2xl font-semibold mb-6">Post-Deployment</h2>

              <div className="space-y-6">
                <Card>
                  <CardHeader>
                    <CardTitle>Generate API Keys</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`# Generate secure API key
API_KEY=$(openssl rand -hex 32)

# Store in AWS Secrets Manager
aws secretsmanager create-secret \\
  --name minion-api-key-$ENVIRONMENT \\
  --secret-string "$API_KEY"

# Retrieve when needed
aws secretsmanager get-secret-value \\
  --secret-id minion-api-key-$ENVIRONMENT \\
  --query SecretString --output text`}
                      language="bash"
                    />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Setup Cost Alerts</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`# Create billing alarm
aws cloudwatch put-metric-alarm \\
  --alarm-name minion-cost-alarm-$ENVIRONMENT \\
  --alarm-description "MinION pipeline cost alert" \\
  --metric-name EstimatedCharges \\
  --namespace AWS/Billing \\
  --statistic Maximum \\
  --period 86400 \\
  --evaluation-periods 1 \\
  --threshold 500 \\
  --comparison-operator GreaterThanThreshold`}
                      language="bash"
                    />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle>Configure Monitoring Dashboard</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground mb-3">
                      Access CloudWatch dashboard at:
                    </p>
                    <code className="bg-muted px-3 py-2 rounded block text-sm">
                      {`https://console.aws.amazon.com/cloudwatch/home?region=$AWS_REGION#dashboards:name=minion-pipeline-$ENVIRONMENT`}
                    </code>
                  </CardContent>
                </Card>
              </div>
            </section>

            <section>
              <h2 className="text-2xl font-semibold mb-6">Troubleshooting</h2>

              <div className="space-y-4">
                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Terraform State Lock</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock code={`terraform force-unlock LOCK_ID`} language="bash" />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">EFS Mount Failed</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`# Check security group rules
aws ec2 describe-security-groups --group-ids sg-xxxxxxxxx

# Add NFS rule (port 2049)
aws ec2 authorize-security-group-ingress \\
  --group-id sg-xxxxxxxxx \\
  --protocol tcp \\
  --port 2049 \\
  --source-group sg-xxxxxxxxx`}
                      language="bash"
                    />
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader>
                    <CardTitle className="text-base">Lambda Timeout</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <CodeBlock
                      code={`# Increase Lambda timeout
aws lambda update-function-configuration \\
  --function-name FUNCTION_NAME \\
  --timeout 300`}
                      language="bash"
                    />
                  </CardContent>
                </Card>
              </div>
            </section>
          </div>
        </main>
      </div>
    </div>
  );
}
