#!/usr/bin/env python3
"""
Trigger report generation phase on EC2
"""

import json
import boto3
import os
from typing import Dict, Any, List

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')
sns = boto3.client('sns')

ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('REPORT_INSTANCE_TYPE', 'm5.xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']
SNS_TOPIC = os.environ['SNS_TOPIC_ARN']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for report generation.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    formats = event.get('config', {}).get('formats', ['pdf', 'json'])

    # Get all analysis outputs
    bucket = f'minion-analysis-{os.environ["ENVIRONMENT"]}'
    base_prefix = f'runs/{run_id}/analysis'

    # Launch instance
    instance_id = launch_report_instance(
        run_id, bucket, base_prefix
    )

    # Wait for instance
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute report generation
    response = execute_report_generation(
        instance_id, run_id, bucket, base_prefix, formats
    )

    return {
        'phase': 'reporting',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': response['Command']['CommandId'],
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{base_prefix}/reports/'
    }

def launch_report_instance(run_id: str, bucket: str, base_prefix: str) -> str:
    """Launch EC2 instance for report generation."""

    user_data = f"""#!/bin/bash
export RUN_ID='{run_id}'
export S3_BASE='s3://{bucket}/{base_prefix}'

# Install LaTeX for PDF generation
apt-get update && apt-get install -y texlive-full pandoc

echo "Report generation instance started for run $RUN_ID" | logger

# Auto-terminate after 2 hours
echo "sudo shutdown -h +120" | at now + 2 hours
"""

    response = ec2.run_instances(
        ImageId=ANALYSIS_AMI,
        InstanceType=INSTANCE_TYPE,
        MinCount=1,
        MaxCount=1,
        SubnetId=SUBNET_ID,
        SecurityGroupIds=[SECURITY_GROUP],
        IamInstanceProfile={'Name': 'MinIONEC2Role'},
        UserData=user_data,
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/xvda',
                'Ebs': {
                    'VolumeSize': 50,
                    'VolumeType': 'gp3',
                    'DeleteOnTermination': True
                }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'reporting-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': 'reporting'}
                ]
            }
        ]
    )

    return response['Instances'][0]['InstanceId']

def execute_report_generation(instance_id: str, run_id: str, bucket: str,
                             base_prefix: str, formats: List[str]) -> Dict:
    """Execute report generation on EC2 instance."""

    # Build format-specific commands
    format_commands = []

    if 'json' in formats:
        format_commands.append("""
# Generate JSON report
/opt/minion/scripts/phase6_reports/generate_pmda_report.py \\
    --input-dir analysis/ \\
    --output-dir reports/ \\
    --format json \\
    --run-id "$RUN_ID"
""")

    if 'pdf' in formats:
        format_commands.append("""
# Generate PDF report
/opt/minion/scripts/phase6_reports/generate_pdf_report.py \\
    --input-dir analysis/ \\
    --output "reports/${RUN_ID}_report.pdf" \\
    --run-id "$RUN_ID"
""")

    if 'html' in formats:
        format_commands.append("""
# Generate HTML report
/opt/minion/scripts/phase6_reports/generate_html_report.py \\
    --input-dir analysis/ \\
    --output "reports/${RUN_ID}_report.html" \\
    --run-id "$RUN_ID"
""")

    command = f"""
#!/bin/bash
set -euo pipefail

export RUN_ID='{run_id}'
WORK_DIR="/mnt/analysis/$RUN_ID/reporting"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Download all analysis results
echo "Downloading analysis results..."
mkdir -p analysis
aws s3 sync "s3://{bucket}/{base_prefix}/" analysis/ \\
    --exclude "*.fastq*" \\
    --exclude "*.fast5" \\
    --exclude "*.pod5"

# Create reports directory
mkdir -p reports

# Generate PMDA compliance checklist (always)
/opt/minion/scripts/phase6_reports/generate_pmda_checklist.py \\
    --input-dir analysis/ \\
    --output reports/pmda_checklist.json \\
    --run-id "$RUN_ID"

{''.join(format_commands)}

# Create final summary
cat > reports/analysis_complete.json << EOF
{{
    "run_id": "$RUN_ID",
    "status": "COMPLETED",
    "timestamp": "$(date -Iseconds)",
    "reports_generated": {formats},
    "pmda_compliance": true
}}
EOF

# Upload reports
aws s3 sync reports/ "s3://{bucket}/{base_prefix}/reports/"

# Generate presigned URL for PDF report (if exists)
if [ -f "reports/${{RUN_ID}}_report.pdf" ]; then
    PDF_URL=$(aws s3 presign "s3://{bucket}/{base_prefix}/reports/${{RUN_ID}}_report.pdf" --expires-in 604800)

    # Send completion notification
    aws sns publish \\
        --topic-arn {SNS_TOPIC} \\
        --subject "Analysis Complete - Run $RUN_ID" \\
        --message "Pathogen screening analysis completed for run $RUN_ID.\\n\\nReport available at: $PDF_URL\\n\\nPMDA compliance checklist generated."
fi

# Clean up
rm -rf "$WORK_DIR"

echo "Report generation completed for run $RUN_ID"
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command.format(formats=json.dumps(formats))],
            'executionTimeout': ['7200']  # 2 hours
        }
    )

    return response