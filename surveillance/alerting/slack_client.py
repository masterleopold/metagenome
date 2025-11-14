"""
Slack Notification Client for 4-Virus Surveillance System

Supports multiple notification methods:
1. Bot Token API (recommended) - Posts to specific channels with rich formatting
2. Incoming Webhooks - Simple fallback method

Environment Variables Required:
- SLACK_BOT_TOKEN: Bot User OAuth Token (xoxb-...)
- SLACK_WEBHOOK_URL: Optional fallback webhook URL
- SLACK_APP_ID: Application ID (for verification)
- SLACK_SIGNING_SECRET: For webhook signature verification
"""

import os
import json
import logging
import requests
from typing import Dict, List, Optional, Any
from datetime import datetime


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SlackClient:
    """
    Client for sending notifications to Slack

    Supports both Bot API and Webhook methods with automatic fallback.
    """

    def __init__(
        self,
        bot_token: Optional[str] = None,
        webhook_url: Optional[str] = None,
        app_id: Optional[str] = None,
        signing_secret: Optional[str] = None
    ):
        """
        Initialize Slack client

        Args:
            bot_token: Bot User OAuth Token (from env if not provided)
            webhook_url: Incoming Webhook URL (from env if not provided)
            app_id: Slack App ID
            signing_secret: Signing secret for verification
        """
        self.bot_token = bot_token or os.environ.get('SLACK_BOT_TOKEN')
        self.webhook_url = webhook_url or os.environ.get('SLACK_WEBHOOK_URL')
        self.app_id = app_id or os.environ.get('SLACK_APP_ID')
        self.signing_secret = signing_secret or os.environ.get('SLACK_SIGNING_SECRET')

        # Slack API endpoint
        self.api_base = "https://slack.com/api"

        # Validate configuration
        if not self.bot_token and not self.webhook_url:
            logger.warning("No Slack credentials configured - notifications will not be sent")
            self.enabled = False
        else:
            self.enabled = True
            logger.info("Initialized Slack client (Bot API: {}, Webhook: {})".format(
                bool(self.bot_token), bool(self.webhook_url)
            ))

    def send_alert(
        self,
        detection: Dict[str, Any],
        channel: str
    ) -> Dict[str, Any]:
        """
        Send virus detection alert to Slack channel

        Args:
            detection: Detection record with severity, virus_type, etc.
            channel: Slack channel name (e.g., #critical-alerts)

        Returns:
            Delivery status
        """
        if not self.enabled:
            logger.warning("Slack notifications disabled - no credentials configured")
            return {
                'status': 'disabled',
                'timestamp': datetime.now().isoformat()
            }

        # Build message payload
        message = self._build_alert_message(detection)

        # Try Bot API first (preferred method)
        if self.bot_token:
            try:
                return self._send_via_bot_api(message, channel)
            except Exception as e:
                logger.error(f"Bot API failed: {e}, falling back to webhook")

        # Fallback to webhook if bot API not available or failed
        if self.webhook_url:
            try:
                return self._send_via_webhook(message)
            except Exception as e:
                logger.error(f"Webhook failed: {e}")
                raise

        raise RuntimeError("No valid Slack notification method available")

    def _send_via_bot_api(
        self,
        message: Dict[str, Any],
        channel: str
    ) -> Dict[str, Any]:
        """
        Send message using Bot Token API (chat.postMessage)

        Args:
            message: Message payload with blocks
            channel: Channel name or ID

        Returns:
            API response
        """
        url = f"{self.api_base}/chat.postMessage"

        headers = {
            "Authorization": f"Bearer {self.bot_token}",
            "Content-Type": "application/json"
        }

        # Ensure channel starts with # or is a channel ID
        if not channel.startswith('#') and not channel.startswith('C'):
            channel = f"#{channel}"

        payload = {
            "channel": channel,
            "text": message.get('text', 'Virus Detection Alert'),  # Fallback text
            "blocks": message.get('blocks', []),
            "unfurl_links": False,
            "unfurl_media": False
        }

        response = requests.post(url, headers=headers, json=payload, timeout=10)
        response.raise_for_status()

        result = response.json()

        if not result.get('ok'):
            error = result.get('error', 'Unknown error')
            logger.error(f"Slack API error: {error}")
            raise RuntimeError(f"Slack API error: {error}")

        logger.info(f"Slack message sent via Bot API: ts={result.get('ts')}")

        return {
            'status': 'sent',
            'method': 'bot_api',
            'channel': channel,
            'timestamp': result.get('ts'),
            'sent_at': datetime.now().isoformat()
        }

    def _send_via_webhook(self, message: Dict[str, Any]) -> Dict[str, Any]:
        """
        Send message using Incoming Webhook

        Args:
            message: Message payload with blocks

        Returns:
            Delivery status
        """
        if not self.webhook_url:
            raise ValueError("No webhook URL configured")

        payload = {
            "text": message.get('text', 'Virus Detection Alert'),
            "blocks": message.get('blocks', [])
        }

        response = requests.post(self.webhook_url, json=payload, timeout=10)
        response.raise_for_status()

        if response.text != 'ok':
            logger.error(f"Webhook error: {response.text}")
            raise RuntimeError(f"Slack webhook error: {response.text}")

        logger.info("Slack message sent via Webhook")

        return {
            'status': 'sent',
            'method': 'webhook',
            'sent_at': datetime.now().isoformat()
        }

    def _build_alert_message(self, detection: Dict[str, Any]) -> Dict[str, Any]:
        """
        Build Slack message with Block Kit formatting

        Args:
            detection: Detection record

        Returns:
            Message payload with blocks
        """
        severity = detection.get('severity', 'unknown')
        virus_type = detection.get('virus_type', 'UNKNOWN').upper()
        copies_per_ml = detection.get('copies_per_ml', 0)
        sample_id = detection.get('sample_id', 'N/A')
        run_id = detection.get('run_id', 'N/A')
        source = detection.get('source', 'unknown')
        reason = detection.get('reason', 'N/A')
        timestamp = detection.get('timestamp', datetime.now().isoformat())

        # Severity-specific formatting
        severity_config = {
            'critical': {
                'emoji': ':rotating_light:',
                'color': '#DC3545',
                'header': 'CRITICAL ALERT'
            },
            'high': {
                'emoji': ':warning:',
                'color': '#FFC107',
                'header': 'HIGH PRIORITY ALERT'
            },
            'medium': {
                'emoji': ':information_source:',
                'color': '#17A2B8',
                'header': 'MEDIUM PRIORITY ALERT'
            },
            'low': {
                'emoji': ':memo:',
                'color': '#6C757D',
                'header': 'INFORMATIONAL ALERT'
            }
        }

        config = severity_config.get(severity, severity_config['low'])

        # Build message blocks
        blocks = [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": f"{config['emoji']} {config['header']}",
                    "emoji": True
                }
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"*{virus_type}* detected in MinION surveillance system"
                }
            },
            {
                "type": "divider"
            },
            {
                "type": "section",
                "fields": [
                    {
                        "type": "mrkdwn",
                        "text": f"*Severity:*\n{severity.upper()}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Viral Load:*\n{copies_per_ml:.2f} copies/mL"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Sample ID:*\n{sample_id}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Run ID:*\n{run_id}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Source:*\n{source}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Timestamp:*\n{timestamp[:19].replace('T', ' ')}"
                    }
                ]
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"*Reason:*\n{reason}"
                }
            }
        ]

        # Add validation requirements for critical/high severity
        if severity in ['critical', 'high']:
            validation_required = detection.get('validation_required', [])
            if validation_required:
                validation_text = "\n".join([f"• {req}" for req in validation_required])
                blocks.append({
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": f"*Validation Required:*\n{validation_text}"
                    }
                })

        # Add action buttons for critical alerts
        if severity == 'critical':
            blocks.extend([
                {
                    "type": "divider"
                },
                {
                    "type": "actions",
                    "elements": [
                        {
                            "type": "button",
                            "text": {
                                "type": "plain_text",
                                "text": "View Dashboard",
                                "emoji": True
                            },
                            "style": "primary",
                            "url": "http://localhost:8501",  # Update with actual dashboard URL
                            "action_id": "view_dashboard"
                        },
                        {
                            "type": "button",
                            "text": {
                                "type": "plain_text",
                                "text": "Acknowledge",
                                "emoji": True
                            },
                            "style": "danger",
                            "action_id": "acknowledge_alert",
                            "value": detection.get('detection_id', 'unknown')
                        }
                    ]
                }
            ])

        # Add footer
        blocks.append({
            "type": "context",
            "elements": [
                {
                    "type": "mrkdwn",
                    "text": f"MinION 4-Virus Surveillance System | Detection ID: {detection.get('detection_id', 'unknown')}"
                }
            ]
        })

        # Fallback text for notifications
        fallback_text = (
            f"{config['emoji']} [{severity.upper()}] {virus_type} detected - "
            f"{copies_per_ml:.0f} copies/mL in {sample_id}"
        )

        return {
            'text': fallback_text,
            'blocks': blocks,
            'attachments': [
                {
                    'color': config['color'],
                    'fallback': fallback_text
                }
            ]
        }

    def send_daily_summary(
        self,
        summary_data: Dict[str, Any],
        channel: str = "#pathogen-monitoring"
    ) -> Dict[str, Any]:
        """
        Send daily summary to Slack

        Args:
            summary_data: Summary statistics
            channel: Target channel

        Returns:
            Delivery status
        """
        if not self.enabled:
            logger.warning("Slack notifications disabled")
            return {'status': 'disabled'}

        message = self._build_summary_message(summary_data)

        if self.bot_token:
            return self._send_via_bot_api(message, channel)
        elif self.webhook_url:
            return self._send_via_webhook(message)

        raise RuntimeError("No Slack notification method available")

    def _build_summary_message(self, summary: Dict[str, Any]) -> Dict[str, Any]:
        """Build daily summary message"""
        date_str = datetime.now().strftime('%Y-%m-%d')

        blocks = [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": f":bar_chart: Daily Surveillance Summary - {date_str}",
                    "emoji": True
                }
            },
            {
                "type": "section",
                "fields": [
                    {
                        "type": "mrkdwn",
                        "text": f"*Samples Processed:*\n{summary.get('samples_processed', 0)}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Total Detections:*\n{summary.get('total_detections', 0)}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*MAFF Reports:*\n{summary.get('maff_new_reports', 0)} new"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f"*Academic Papers:*\n{summary.get('academic_new_papers', 0)} new"
                    }
                ]
            },
            {
                "type": "divider"
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": "*Alert Summary:*"
                }
            },
            {
                "type": "section",
                "fields": [
                    {
                        "type": "mrkdwn",
                        "text": f":rotating_light: Critical: {summary.get('alerts_critical', 0)}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f":warning: High: {summary.get('alerts_high', 0)}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f":information_source: Medium: {summary.get('alerts_medium', 0)}"
                    },
                    {
                        "type": "mrkdwn",
                        "text": f":memo: Low: {summary.get('alerts_low', 0)}"
                    }
                ]
            },
            {
                "type": "context",
                "elements": [
                    {
                        "type": "mrkdwn",
                        "text": "MinION 4-Virus Surveillance System"
                    }
                ]
            }
        ]

        return {
            'text': f'Daily Surveillance Summary - {date_str}',
            'blocks': blocks
        }

    def test_connection(self) -> bool:
        """
        Test Slack connection

        Returns:
            True if connection successful
        """
        if not self.enabled:
            logger.error("Slack not configured")
            return False

        if self.bot_token:
            try:
                url = f"{self.api_base}/auth.test"
                headers = {"Authorization": f"Bearer {self.bot_token}"}
                response = requests.post(url, headers=headers, timeout=10)
                result = response.json()

                if result.get('ok'):
                    logger.info(f"Slack connection successful: team={result.get('team')}, "
                              f"user={result.get('user')}")
                    return True
                else:
                    logger.error(f"Slack auth failed: {result.get('error')}")
                    return False

            except Exception as e:
                logger.error(f"Slack connection test failed: {e}")
                return False

        # For webhook, just verify URL is configured
        return bool(self.webhook_url)


if __name__ == "__main__":
    # Example usage and testing
    import sys

    # Set test credentials (normally from environment)
    os.environ['SLACK_BOT_TOKEN'] = 'xoxb-your-token-here'  # Replace with actual token
    os.environ['SLACK_APP_ID'] = 'A09TVLTGDSL'

    client = SlackClient()

    # Test connection
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        success = client.test_connection()
        print(f"Connection test: {'SUCCESS' if success else 'FAILED'}")
        sys.exit(0 if success else 1)

    # Test alert
    test_detection = {
        'detection_id': 'test-123',
        'virus_type': 'spumavirus',
        'severity': 'high',
        'source': 'internal_pipeline_kraken2',
        'copies_per_ml': 250.5,
        'sample_id': 'SAMPLE-001',
        'run_id': 'RUN-001',
        'timestamp': datetime.now().isoformat(),
        'reason': 'High viral load detected',
        'validation_required': ['Sequence confirmation', 'Repeat testing']
    }

    if client.enabled:
        try:
            result = client.send_alert(test_detection, '#test-alerts')
            print(f"Alert sent: {json.dumps(result, indent=2)}")
        except Exception as e:
            print(f"Failed to send alert: {e}")
    else:
        print("Slack not configured - set SLACK_BOT_TOKEN or SLACK_WEBHOOK_URL")
