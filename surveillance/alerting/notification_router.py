"""
Notification Router for 4-Virus Surveillance System

Routes alerts to appropriate channels (SNS, SES, SMS, Slack) based on severity level.
Handles deduplication, rate limiting, and delivery tracking.
"""

import os
import json
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
from pathlib import Path
import boto3
from botocore.exceptions import ClientError

# Import Slack client
from .slack_client import SlackClient


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class NotificationRouter:
    """
    Router for distributing virus detection alerts

    Manages multi-channel notifications based on severity levels
    and handles deduplication to prevent alert fatigue.
    """

    def __init__(self, region: str = "ap-northeast-1"):
        """
        Initialize notification router

        Args:
            region: AWS region
        """
        self.region = region

        # Initialize AWS clients
        self.sns_client = boto3.client('sns', region_name=region)
        self.ses_client = boto3.client('ses', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)

        # Initialize Slack client
        self.slack_client = SlackClient()

        # Notification history for deduplication
        self.notification_history: List[Dict[str, Any]] = []

        logger.info("Initialized Notification Router (Slack: {})".format(
            "enabled" if self.slack_client.enabled else "disabled"
        ))

    def route_alert(self, detection: Dict[str, Any]) -> Dict[str, Any]:
        """
        Route alert based on detection severity

        Args:
            detection: Detection record with severity and notification_config

        Returns:
            Delivery status for each channel
        """
        severity = detection.get('severity', 'low')
        notification_config = detection.get('notification_config', {})

        logger.info(f"Routing {severity} alert for {detection.get('virus_type')}")

        delivery_status = {
            'timestamp': datetime.now().isoformat(),
            'detection_id': detection.get('detection_id'),
            'severity': severity,
            'channels': {}
        }

        # Check if should deduplicate
        if detection.get('alert_sent', False):
            logger.info(f"Alert already sent for detection {detection.get('detection_id')}")
            delivery_status['deduplicated'] = True
            return delivery_status

        # Send to SNS topics
        sns_topics = notification_config.get('sns_topics', [])
        for topic_arn in sns_topics:
            try:
                result = self._send_sns(detection, topic_arn)
                delivery_status['channels'][f"sns_{topic_arn.split(':')[-1]}"] = result
            except Exception as e:
                logger.error(f"Failed to send SNS to {topic_arn}: {e}")
                delivery_status['channels'][f"sns_{topic_arn.split(':')[-1]}"] = {'error': str(e)}

        # Send emails via SES
        email_recipients = notification_config.get('email', [])
        if email_recipients:
            try:
                result = self._send_email(detection, email_recipients)
                delivery_status['channels']['email'] = result
            except Exception as e:
                logger.error(f"Failed to send email: {e}")
                delivery_status['channels']['email'] = {'error': str(e)}

        # Send SMS (critical only)
        if severity == 'critical':
            sms_recipients = notification_config.get('sms', [])
            for phone_number in sms_recipients:
                try:
                    result = self._send_sms(detection, phone_number)
                    delivery_status['channels'][f"sms_{phone_number}"] = result
                except Exception as e:
                    logger.error(f"Failed to send SMS to {phone_number}: {e}")
                    delivery_status['channels'][f"sms_{phone_number}"] = {'error': str(e)}

        # Send to Slack (if configured)
        slack_channels = notification_config.get('slack', [])
        for channel in slack_channels:
            try:
                result = self._send_slack(detection, channel)
                delivery_status['channels'][f"slack_{channel}"] = result
            except Exception as e:
                logger.error(f"Failed to send Slack to {channel}: {e}")
                delivery_status['channels'][f"slack_{channel}"] = {'error': str(e)}

        # Mark detection as alerted
        detection['alert_sent'] = True
        detection['alert_timestamp'] = datetime.now().isoformat()

        # Save notification record
        self._save_notification_record(detection, delivery_status)

        return delivery_status

    def _send_sns(self, detection: Dict[str, Any], topic_arn: str) -> Dict[str, str]:
        """
        Send SNS notification

        Args:
            detection: Detection record
            topic_arn: SNS topic ARN

        Returns:
            Delivery status
        """
        message = self._format_sns_message(detection)
        subject = self._format_subject(detection)

        try:
            response = self.sns_client.publish(
                TopicArn=topic_arn,
                Subject=subject,
                Message=message,
                MessageAttributes={
                    'severity': {
                        'DataType': 'String',
                        'StringValue': detection.get('severity', 'unknown')
                    },
                    'virus_type': {
                        'DataType': 'String',
                        'StringValue': detection.get('virus_type', 'unknown')
                    }
                }
            )

            logger.info(f"SNS sent: MessageId={response['MessageId']}")
            return {
                'status': 'sent',
                'message_id': response['MessageId'],
                'timestamp': datetime.now().isoformat()
            }

        except ClientError as e:
            logger.error(f"SNS publish failed: {e}")
            raise

    def _send_email(self, detection: Dict[str, Any], recipients: List[str]) -> Dict[str, Any]:
        """
        Send email via SES

        Args:
            detection: Detection record
            recipients: List of email addresses

        Returns:
            Delivery status
        """
        subject = self._format_subject(detection)
        body_html = self._format_email_html(detection)
        body_text = self._format_email_text(detection)

        try:
            response = self.ses_client.send_email(
                Source='surveillance@example.com',  # Configure sender
                Destination={
                    'ToAddresses': recipients
                },
                Message={
                    'Subject': {
                        'Data': subject,
                        'Charset': 'UTF-8'
                    },
                    'Body': {
                        'Text': {
                            'Data': body_text,
                            'Charset': 'UTF-8'
                        },
                        'Html': {
                            'Data': body_html,
                            'Charset': 'UTF-8'
                        }
                    }
                }
            )

            logger.info(f"Email sent: MessageId={response['MessageId']}")
            return {
                'status': 'sent',
                'message_id': response['MessageId'],
                'recipients': recipients,
                'timestamp': datetime.now().isoformat()
            }

        except ClientError as e:
            logger.error(f"SES send failed: {e}")
            raise

    def _send_sms(self, detection: Dict[str, Any], phone_number: str) -> Dict[str, str]:
        """
        Send SMS via SNS

        Args:
            detection: Detection record
            phone_number: Phone number (E.164 format)

        Returns:
            Delivery status
        """
        message = self._format_sms_message(detection)

        try:
            response = self.sns_client.publish(
                PhoneNumber=phone_number,
                Message=message,
                MessageAttributes={
                    'AWS.SNS.SMS.SenderID': {
                        'DataType': 'String',
                        'StringValue': 'MinION'
                    },
                    'AWS.SNS.SMS.SMSType': {
                        'DataType': 'String',
                        'StringValue': 'Transactional'
                    }
                }
            )

            logger.info(f"SMS sent to {phone_number}: MessageId={response['MessageId']}")
            return {
                'status': 'sent',
                'message_id': response['MessageId'],
                'timestamp': datetime.now().isoformat()
            }

        except ClientError as e:
            logger.error(f"SMS send failed: {e}")
            raise

    def _send_slack(self, detection: Dict[str, Any], channel: str) -> Dict[str, str]:
        """
        Send Slack notification using Slack client

        Args:
            detection: Detection record
            channel: Slack channel name (e.g., #critical-alerts)

        Returns:
            Delivery status
        """
        if not self.slack_client.enabled:
            logger.warning(f"Slack client disabled - skipping notification to {channel}")
            return {
                'status': 'disabled',
                'channel': channel,
                'timestamp': datetime.now().isoformat()
            }

        try:
            result = self.slack_client.send_alert(detection, channel)
            logger.info(f"Slack notification sent to {channel}: {result.get('status')}")
            return result

        except Exception as e:
            logger.error(f"Failed to send Slack notification to {channel}: {e}")
            raise

    def _format_subject(self, detection: Dict[str, Any]) -> str:
        """
        Format alert subject line

        Args:
            detection: Detection record

        Returns:
            Subject string
        """
        severity = detection.get('severity', 'UNKNOWN').upper()
        virus = detection.get('virus_type', 'UNKNOWN').upper()
        sample = detection.get('sample_id', detection.get('run_id', 'UNKNOWN'))

        # Add emoji for visual scanning
        emoji = {
            'CRITICAL': '🚨',
            'HIGH': '⚠️',
            'MEDIUM': 'ℹ️',
            'LOW': '📝'
        }.get(severity, '')

        return f"{emoji} [{severity}] {virus} Detected - {sample}"

    def _format_sns_message(self, detection: Dict[str, Any]) -> str:
        """Format message for SNS notification"""
        return json.dumps({
            'alert_type': 'virus_detection',
            'severity': detection.get('severity'),
            'virus_type': detection.get('virus_type'),
            'sample_id': detection.get('sample_id'),
            'run_id': detection.get('run_id'),
            'source': detection.get('source'),
            'copies_per_ml': detection.get('copies_per_ml'),
            'timestamp': detection.get('timestamp'),
            'reason': detection.get('reason'),
            'detection_id': detection.get('detection_id')
        }, indent=2)

    def _format_sms_message(self, detection: Dict[str, Any]) -> str:
        """Format concise SMS message (160 char limit)"""
        severity = detection.get('severity', 'UNK').upper()
        virus = detection.get('virus_type', 'UNKNOWN')
        copies = detection.get('copies_per_ml', 0)

        return f"ALERT [{severity}]: {virus} detected - {copies:.0f} copies/mL. Check dashboard immediately."

    def _format_email_text(self, detection: Dict[str, Any]) -> str:
        """Format plain text email body"""
        return f"""
4-Virus Surveillance System Alert
=================================

Severity: {detection.get('severity', 'unknown').upper()}
Virus Type: {detection.get('virus_type', 'unknown')}
Sample ID: {detection.get('sample_id', 'N/A')}
Run ID: {detection.get('run_id', 'N/A')}

Detection Details:
------------------
Source: {detection.get('source', 'unknown')}
Estimated Viral Load: {detection.get('copies_per_ml', 0):.2f} copies/mL
Detection Time: {detection.get('timestamp', 'unknown')}
Reason: {detection.get('reason', 'N/A')}

Metadata:
---------
{json.dumps(detection.get('metadata', {}), indent=2)}

Validation Required:
--------------------
{chr(10).join('- ' + req for req in detection.get('validation_required', ['None']))}

Next Steps:
-----------
1. Review detection in real-time dashboard
2. Verify sequence data quality
3. Check for external confirmation
4. Follow validation protocol for {detection.get('severity')} severity

--
MinION 4-Virus Surveillance System
Generated: {datetime.now().isoformat()}
        """.strip()

    def _format_email_html(self, detection: Dict[str, Any]) -> str:
        """Format HTML email body"""
        severity = detection.get('severity', 'unknown')
        severity_colors = {
            'critical': '#DC3545',
            'high': '#FFC107',
            'medium': '#17A2B8',
            'low': '#6C757D'
        }
        color = severity_colors.get(severity, '#6C757D')

        return f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; }}
        .header {{ background-color: {color}; color: white; padding: 20px; }}
        .content {{ padding: 20px; }}
        .severity {{ font-size: 24px; font-weight: bold; }}
        .details {{ background-color: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid {color}; }}
        .footer {{ padding: 20px; background-color: #f1f1f1; font-size: 12px; color: #666; }}
    </style>
</head>
<body>
    <div class="header">
        <div class="severity">{severity.upper()} ALERT</div>
        <p>{detection.get('virus_type', 'UNKNOWN').upper()} Detected</p>
    </div>

    <div class="content">
        <h2>Detection Details</h2>
        <div class="details">
            <p><strong>Sample ID:</strong> {detection.get('sample_id', 'N/A')}</p>
            <p><strong>Run ID:</strong> {detection.get('run_id', 'N/A')}</p>
            <p><strong>Source:</strong> {detection.get('source', 'unknown')}</p>
            <p><strong>Viral Load:</strong> {detection.get('copies_per_ml', 0):.2f} copies/mL</p>
            <p><strong>Time:</strong> {detection.get('timestamp', 'unknown')}</p>
            <p><strong>Reason:</strong> {detection.get('reason', 'N/A')}</p>
        </div>

        <h3>Required Validation Steps</h3>
        <ul>
            {''.join(f"<li>{req}</li>" for req in detection.get('validation_required', ['None']))}
        </ul>

        <h3>Next Actions</h3>
        <ol>
            <li>Review detection in real-time dashboard</li>
            <li>Verify sequence data quality</li>
            <li>Check for external confirmation</li>
            <li>Follow validation protocol for {severity} severity</li>
        </ol>
    </div>

    <div class="footer">
        <p>MinION 4-Virus Surveillance System</p>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Detection ID: {detection.get('detection_id', 'unknown')}</p>
    </div>
</body>
</html>
        """.strip()

    def _save_notification_record(
        self,
        detection: Dict[str, Any],
        delivery_status: Dict[str, Any]
    ) -> None:
        """
        Save notification record to DynamoDB for tracking

        Args:
            detection: Detection record
            delivery_status: Delivery status for each channel
        """
        try:
            table = self.dynamodb.Table('surveillance-notifications')

            table.put_item(Item={
                'notification_id': detection.get('detection_id'),
                'timestamp': datetime.now().isoformat(),
                'virus_type': detection.get('virus_type'),
                'severity': detection.get('severity'),
                'delivery_status': delivery_status,
                'ttl': int((datetime.now() + timedelta(days=90)).timestamp())  # 90-day retention
            })

            logger.info(f"Saved notification record: {detection.get('detection_id')}")

        except ClientError as e:
            logger.error(f"Failed to save notification record: {e}")

    def send_daily_summary(self, summary_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Send daily summary via email and Slack

        Args:
            summary_data: Summary of daily surveillance activities

        Returns:
            Delivery status for each channel
        """
        logger.info("Sending daily summary")

        delivery_status = {
            'timestamp': datetime.now().isoformat(),
            'channels': {}
        }

        # Send email summary
        recipients = ['surveillance-team@example.com']  # Configure
        subject = f"Daily 4-Virus Surveillance Summary - {datetime.now().strftime('%Y-%m-%d')}"

        body_text = self._format_daily_summary_text(summary_data)
        body_html = self._format_daily_summary_html(summary_data)

        try:
            response = self.ses_client.send_email(
                Source='surveillance@example.com',
                Destination={'ToAddresses': recipients},
                Message={
                    'Subject': {'Data': subject, 'Charset': 'UTF-8'},
                    'Body': {
                        'Text': {'Data': body_text, 'Charset': 'UTF-8'},
                        'Html': {'Data': body_html, 'Charset': 'UTF-8'}
                    }
                }
            )

            delivery_status['channels']['email'] = {
                'status': 'sent',
                'message_id': response['MessageId']
            }

        except ClientError as e:
            logger.error(f"Failed to send email summary: {e}")
            delivery_status['channels']['email'] = {'error': str(e)}

        # Send Slack summary
        if self.slack_client.enabled:
            try:
                slack_result = self.slack_client.send_daily_summary(
                    summary_data,
                    channel='#pathogen-monitoring'
                )
                delivery_status['channels']['slack'] = slack_result

            except Exception as e:
                logger.error(f"Failed to send Slack summary: {e}")
                delivery_status['channels']['slack'] = {'error': str(e)}

        return delivery_status

    def _format_daily_summary_text(self, summary: Dict[str, Any]) -> str:
        """Format daily summary as plain text"""
        return f"""
Daily 4-Virus Surveillance Summary
Date: {datetime.now().strftime('%Y-%m-%d')}
====================================

External Sources Checked:
-------------------------
MAFF Reports: {summary.get('maff_reports_checked', 0)} checked, {summary.get('maff_new_reports', 0)} new
E-Stat Statistics: {summary.get('estat_tables_checked', 0)} tables queried
Academic Publications: {summary.get('academic_new_papers', 0)} new papers

Internal Pipeline Detections:
------------------------------
Total Samples Processed: {summary.get('samples_processed', 0)}
Virus Detections: {summary.get('total_detections', 0)}

Detection Breakdown:
{json.dumps(summary.get('detections_by_virus', {}), indent=2)}

Alert Summary:
--------------
Critical: {summary.get('alerts_critical', 0)}
High: {summary.get('alerts_high', 0)}
Medium: {summary.get('alerts_medium', 0)}
Low: {summary.get('alerts_low', 0)}

--
MinION 4-Virus Surveillance System
        """.strip()

    def _format_daily_summary_html(self, summary: Dict[str, Any]) -> str:
        """Format daily summary as HTML"""
        # Simplified HTML format
        return f"""
<!DOCTYPE html>
<html>
<body style="font-family: Arial, sans-serif;">
    <h1>Daily 4-Virus Surveillance Summary</h1>
    <p><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d')}</p>

    <h2>External Sources</h2>
    <ul>
        <li>MAFF: {summary.get('maff_new_reports', 0)} new reports</li>
        <li>E-Stat: {summary.get('estat_tables_checked', 0)} tables</li>
        <li>Academic: {summary.get('academic_new_papers', 0)} papers</li>
    </ul>

    <h2>Pipeline Activity</h2>
    <p>Samples Processed: {summary.get('samples_processed', 0)}</p>
    <p>Total Detections: {summary.get('total_detections', 0)}</p>

    <h2>Alerts</h2>
    <ul>
        <li>Critical: {summary.get('alerts_critical', 0)}</li>
        <li>High: {summary.get('alerts_high', 0)}</li>
        <li>Medium: {summary.get('alerts_medium', 0)}</li>
        <li>Low: {summary.get('alerts_low', 0)}</li>
    </ul>
</body>
</html>
        """.strip()


if __name__ == "__main__":
    # Example usage
    router = NotificationRouter()

    test_detection = {
        'detection_id': 'test-123',
        'virus_type': 'spumavirus',
        'severity': 'high',
        'source': 'internal_pipeline',
        'copies_per_ml': 250,
        'sample_id': 'SAMPLE-001',
        'run_id': 'RUN-001',
        'timestamp': datetime.now().isoformat(),
        'reason': 'High viral load detected',
        'validation_required': ['Sequence confirmation', 'Repeat testing'],
        'notification_config': {
            'email': ['test@example.com']
        }
    }

    # status = router.route_alert(test_detection)
    # print(json.dumps(status, indent=2))
