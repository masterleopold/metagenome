#!/usr/bin/env python3
"""
Handle critical alerts and notifications
"""

import json
import boto3
import os
from typing import Dict, Any
from datetime import datetime

sns = boto3.client('sns')
ses = boto3.client('ses')
rds = boto3.client('rds-data')

SNS_TOPIC = os.environ['SNS_TOPIC_ARN']
ALERT_EMAIL = os.environ.get('ALERT_EMAIL', '')
CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']

# Critical pathogen codes requiring immediate alert
CRITICAL_PATHOGENS = [
    'PERV-A', 'PERV-B', 'PERV-C',  # PERV variants
    'ASFV',  # African swine fever virus
    'CSFV',  # Classical swine fever virus
    'FMDV',  # Foot-and-mouth disease virus
    'JEV',   # Japanese encephalitis virus
    'HEV'    # Hepatitis E virus
]

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Handle alerts from various sources.

    Triggered by:
    1. PERV detection
    2. Critical pathogen detection
    3. QC failure
    4. Workflow errors
    """

    alert_type = event.get('alert_type')
    run_id = event.get('run_id')
    workflow_id = event.get('workflow_id')
    details = event.get('details', {})

    # Process alert based on type
    if alert_type == 'PERV_DETECTION':
        handle_perv_alert(run_id, workflow_id, details)

    elif alert_type == 'CRITICAL_PATHOGEN':
        handle_critical_pathogen(run_id, workflow_id, details)

    elif alert_type == 'QC_FAILURE':
        handle_qc_failure(run_id, workflow_id, details)

    elif alert_type == 'WORKFLOW_ERROR':
        handle_workflow_error(run_id, workflow_id, details)

    elif alert_type == 'HIGH_CONTAMINATION':
        handle_high_contamination(run_id, workflow_id, details)

    else:
        # Generic alert
        send_generic_alert(run_id, workflow_id, alert_type, details)

    # Log alert to database
    log_alert(workflow_id, alert_type, details)

    return {
        'status': 'ALERT_SENT',
        'alert_type': alert_type,
        'run_id': run_id,
        'timestamp': datetime.utcnow().isoformat()
    }

def handle_perv_alert(run_id: str, workflow_id: str, details: Dict):
    """Handle PERV detection alert - highest priority."""

    subject = f'[CRITICAL] PERV Detection - Run {run_id}'

    message = f"""
CRITICAL ALERT: Porcine Endogenous Retrovirus (PERV) Detected

Run ID: {run_id}
Workflow ID: {workflow_id}
Detection Time: {datetime.utcnow().isoformat()}

PERV Details:
- Subtype: {details.get('perv_subtype', 'Unknown')}
- Read Count: {details.get('read_count', 0)}
- Confidence: {details.get('confidence', 'Unknown')}
- Copy Number: {details.get('copy_number', 0)} copies/mL

IMMEDIATE ACTIONS REQUIRED:
1. Quarantine the donor animal
2. Halt any scheduled procedures
3. Perform confirmatory PCR testing
4. Review contamination protocols
5. Notify regulatory authorities

Analysis Reports:
{details.get('report_url', 'Not available')}

This is an automated alert from the MinION Pathogen Screening System.
PMDA compliance requires immediate action on PERV detection.
"""

    # Send high-priority notifications
    send_sns_alert(subject, message)
    if ALERT_EMAIL:
        send_email_alert(subject, message, priority='High')

    # Send to multiple endpoints if configured
    if 'EMERGENCY_PHONE' in os.environ:
        send_sms_alert(os.environ['EMERGENCY_PHONE'], f'PERV detected in run {run_id}')

def handle_critical_pathogen(run_id: str, workflow_id: str, details: Dict):
    """Handle critical pathogen detection."""

    pathogen = details.get('pathogen_code', 'Unknown')
    pathogen_name = details.get('pathogen_name', pathogen)

    subject = f'[WARNING] Critical Pathogen Detected: {pathogen} - Run {run_id}'

    message = f"""
CRITICAL PATHOGEN DETECTION ALERT

Pathogen: {pathogen_name} ({pathogen})
Run ID: {run_id}
Workflow ID: {workflow_id}
Detection Time: {datetime.utcnow().isoformat()}

Detection Details:
- Read Count: {details.get('read_count', 0)}
- Confidence Level: {details.get('confidence', 'Unknown')}
- Detection Method: {details.get('method', 'Unknown')}
- Copy Number: {details.get('copy_number', 0)} copies/mL

Required Actions:
1. Verify detection with orthogonal method
2. Review sample handling procedures
3. Check for potential contamination
4. Document in compliance reports
5. Notify designated personnel

Full Report: {details.get('report_url', 'Not available')}
"""

    send_sns_alert(subject, message)
    if ALERT_EMAIL:
        send_email_alert(subject, message)

def handle_qc_failure(run_id: str, workflow_id: str, details: Dict):
    """Handle QC failure alert."""

    subject = f'QC Failure - Run {run_id}'

    message = f"""
QUALITY CONTROL FAILURE

Run ID: {run_id}
Workflow ID: {workflow_id}
Failure Time: {datetime.utcnow().isoformat()}

QC Metrics:
- Mean Quality Score: {details.get('mean_quality', 'N/A')}
- Reads Passed: {details.get('reads_passed', 0)}
- Reads Failed: {details.get('reads_failed', 0)}
- Failure Reason: {details.get('reason', 'Quality below threshold')}

Recommended Actions:
1. Review sequencing parameters
2. Check sample quality
3. Verify library preparation
4. Consider re-sequencing

The workflow will continue but results may be unreliable.
"""

    send_sns_alert(subject, message)

def handle_workflow_error(run_id: str, workflow_id: str, details: Dict):
    """Handle workflow execution error."""

    phase = details.get('phase', 'Unknown')
    error_message = details.get('error', 'Unknown error')

    subject = f'Workflow Error - {phase} - Run {run_id}'

    message = f"""
WORKFLOW EXECUTION ERROR

Run ID: {run_id}
Workflow ID: {workflow_id}
Failed Phase: {phase}
Error Time: {datetime.utcnow().isoformat()}

Error Details:
{error_message}

Instance ID: {details.get('instance_id', 'N/A')}
Command ID: {details.get('command_id', 'N/A')}

The workflow has been stopped. Manual intervention required.
"""

    send_sns_alert(subject, message)

def handle_high_contamination(run_id: str, workflow_id: str, details: Dict):
    """Handle high host contamination alert."""

    contamination_rate = details.get('contamination_rate', 0)

    subject = f'High Host Contamination - Run {run_id}'

    message = f"""
HIGH HOST CONTAMINATION DETECTED

Run ID: {run_id}
Workflow ID: {workflow_id}
Detection Time: {datetime.utcnow().isoformat()}

Contamination Details:
- Host DNA: {contamination_rate:.1f}%
- Reads Before Depletion: {details.get('reads_before', 0)}
- Reads After Depletion: {details.get('reads_after', 0)}
- Depletion Efficiency: {details.get('depletion_rate', 0):.1f}%

This high contamination level may affect pathogen detection sensitivity.

Recommended Actions:
1. Review host depletion protocol
2. Check sample collection procedure
3. Consider re-extraction with improved depletion
4. Adjust detection thresholds
"""

    if contamination_rate > 50:  # Very high contamination
        send_sns_alert(subject, message)
        if ALERT_EMAIL:
            send_email_alert(subject, message)
    else:
        # Just log for moderate contamination
        send_sns_alert(subject, message)

def send_generic_alert(run_id: str, workflow_id: str, alert_type: str, details: Dict):
    """Send generic alert."""

    subject = f'Alert: {alert_type} - Run {run_id}'
    message = f"""
Alert Type: {alert_type}
Run ID: {run_id}
Workflow ID: {workflow_id}
Time: {datetime.utcnow().isoformat()}

Details:
{json.dumps(details, indent=2, default=str)}
"""

    send_sns_alert(subject, message)

def send_sns_alert(subject: str, message: str):
    """Send alert via SNS."""

    sns.publish(
        TopicArn=SNS_TOPIC,
        Subject=subject,
        Message=message
    )

def send_email_alert(subject: str, message: str, priority: str = 'Normal'):
    """Send alert via SES email."""

    if not ALERT_EMAIL:
        return

    try:
        ses.send_email(
            Source=f'MinION Alerts <alerts@{os.environ["DOMAIN"]}>',
            Destination={'ToAddresses': [ALERT_EMAIL]},
            Message={
                'Subject': {'Data': subject},
                'Body': {'Text': {'Data': message}}
            },
            ReplyToAddresses=['noreply@' + os.environ.get('DOMAIN', 'example.com')],
            Tags=[
                {'Name': 'Priority', 'Value': priority},
                {'Name': 'System', 'Value': 'MinION'}
            ]
        )
    except Exception as e:
        print(f'Error sending email alert: {e}')

def send_sms_alert(phone_number: str, message: str):
    """Send SMS alert for critical issues."""

    try:
        sns.publish(
            PhoneNumber=phone_number,
            Message=message[:160]  # SMS limit
        )
    except Exception as e:
        print(f'Error sending SMS alert: {e}')

def log_alert(workflow_id: str, alert_type: str, details: Dict):
    """Log alert to database."""

    sql = """
        INSERT INTO alerts
        (workflow_id, alert_type, details, created_at)
        VALUES (:workflow_id, :alert_type, :details, NOW())
    """

    # Convert workflow_id to int, handling both string and int types
    try:
        wf_id = int(workflow_id) if workflow_id else 0
    except (ValueError, TypeError):
        wf_id = 0
        print(f'Warning: Invalid workflow_id "{workflow_id}", using 0')

    try:
        rds.execute_statement(
            resourceArn=CLUSTER_ARN,
            secretArn=SECRET_ARN,
            database=DATABASE,
            sql=sql,
            parameters=[
                {'name': 'workflow_id', 'value': {'longValue': wf_id}},
                {'name': 'alert_type', 'value': {'stringValue': alert_type}},
                {'name': 'details', 'value': {'stringValue': json.dumps(details, default=str)}}
            ]
        )
    except Exception as e:
        print(f'Error logging alert to database: {e}')