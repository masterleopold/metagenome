# Slack Integration Setup Guide

Complete guide for setting up Slack notifications in the 4-Virus Surveillance System.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Slack Workspace Setup](#slack-workspace-setup)
- [Configuration](#configuration)
- [Testing](#testing)
- [Deployment](#deployment)
- [Troubleshooting](#troubleshooting)

---

## Overview

The Slack integration provides real-time notifications for virus detections with severity-based routing:

| Severity | Channel | Notifications |
|----------|---------|---------------|
| **CRITICAL** | `#critical-alerts` | Immediate alerts with action buttons |
| **HIGH** | `#pathogen-alerts` | Priority notifications |
| **MEDIUM** | `#pathogen-monitoring` | Standard monitoring alerts |
| **LOW** | Dashboard only | No Slack notification |

### Features

- **Rich Formatting**: Block Kit messages with detection details
- **Severity-Based Routing**: Automatic channel selection
- **Action Buttons**: Quick links to dashboard and acknowledgment
- **Daily Summaries**: Automated daily reports
- **Dual Methods**: Bot API (preferred) with Webhook fallback

---

## Prerequisites

### 1. Slack Workspace

- Admin access to Slack workspace
- Permission to create apps and channels

### 2. Required Channels

Create these channels in your Slack workspace:

```bash
#critical-alerts          # Critical severity (PERV-level)
#pathogen-alerts          # High severity
#pathogen-monitoring      # Medium severity + daily summaries
```

### 3. Python Dependencies

```bash
pip install requests boto3 pyyaml
```

---

## Slack Workspace Setup

### Step 1: Create Slack App

1. Go to https://api.slack.com/apps
2. Click **"Create New App"** → **"From scratch"**
3. App Name: `4-Virus Surveillance`
4. Workspace: Select your workspace
5. Click **"Create App"**

Your App ID: `A09TVLTGDSL` (already created)

### Step 2: Configure Bot User

1. In your app settings, go to **"OAuth & Permissions"**
2. Scroll to **"Scopes"** → **"Bot Token Scopes"**
3. Add these scopes:
   - `chat:write` - Post messages
   - `chat:write.public` - Post to channels without joining
   - `channels:read` - List public channels

### Step 3: Install App to Workspace

1. In app settings, go to **"Install App"**
2. Click **"Install to Workspace"**
3. Review permissions and click **"Allow"**
4. Copy the **Bot User OAuth Token** (starts with `xoxb-`)

### Step 4: Get Credentials

You should now have:

```
App ID: A09TVLTGDSL
Client ID: 580133244533.9947707557904
Client Secret: 78bbc76e4c6e2cb7d591d6b4a359d370
Signing Secret: 5ebeb36dadc1df3a68714befc819c615
Bot Token: xoxb-XXXXXXXXXXXX-XXXXXXXXXXXX-XXXXXXXXXXXXXXXXXXXXXXXX
```

### Step 5: Optional - Incoming Webhook (Fallback)

For webhook-based fallback:

1. Go to **"Incoming Webhooks"**
2. Activate Incoming Webhooks: **ON**
3. Click **"Add New Webhook to Workspace"**
4. Select `#critical-alerts` channel
5. Copy the Webhook URL

---

## Configuration

### Local Development

#### 1. Environment Variables

Copy the template and fill in your credentials:

```bash
cd surveillance/
cp .env.template .env
```

Edit `.env`:

```bash
# Slack Configuration
SLACK_BOT_TOKEN=xoxb-your-actual-token-here
SLACK_WEBHOOK_URL=https://hooks.slack.com/services/T.../B.../...  # Optional
SLACK_APP_ID=A09TVLTGDSL
SLACK_CLIENT_ID=580133244533.9947707557904
SLACK_CLIENT_SECRET=78bbc76e4c6e2cb7d591d6b4a359d370
SLACK_SIGNING_SECRET=5ebeb36dadc1df3a68714befc819c615
```

#### 2. Load Environment

```bash
export $(cat .env | grep -v '^#' | xargs)
```

### AWS Lambda Deployment

#### Option A: Automated Script

```bash
cd surveillance/scripts/
./setup_lambda_env.sh
```

This script will:
- Load credentials from `.env`
- Update all Lambda functions
- Verify configuration

#### Option B: Manual Configuration

For each Lambda function (`surveillance-external-collector`, `surveillance-pipeline-listener`):

```bash
aws lambda update-function-configuration \
  --function-name surveillance-external-collector \
  --region ap-northeast-1 \
  --environment "Variables={
    SLACK_BOT_TOKEN=xoxb-your-token,
    SLACK_APP_ID=A09TVLTGDSL,
    SLACK_CLIENT_ID=580133244533.9947707557904,
    SLACK_CLIENT_SECRET=78bbc76e4c6e2cb7d591d6b4a359d370,
    SLACK_SIGNING_SECRET=5ebeb36dadc1df3a68714befc819c615,
    ENABLE_SLACK_NOTIFICATIONS=true
  }"
```

### Configuration Files

The system uses two configuration files:

#### `config/config.yaml`

```yaml
slack:
  enabled: true
  app_id: A09TVLTGDSL
  channels:
    critical: "#critical-alerts"
    high: "#pathogen-alerts"
    medium: "#pathogen-monitoring"
    daily_summary: "#pathogen-monitoring"
  dashboard_url: "http://localhost:8501"  # Update for production
```

#### `config/severity_rules.yaml`

```yaml
alert_routing:
  critical:
    slack:
      - "#critical-alerts"
  high:
    slack:
      - "#pathogen-alerts"
  medium:
    slack:
      - "#pathogen-monitoring"
```

---

## Testing

### Test 1: Connection Test

Verify Slack credentials and connectivity:

```bash
cd surveillance/tests/
python test_slack_integration.py --test-conn
```

Expected output:
```
✅ PASSED: Slack connection successful
```

### Test 2: Send Test Alerts

Send test notifications to all severity channels:

```bash
python test_slack_integration.py --test-alert
```

This sends:
- Critical alert → `#critical-alerts`
- High alert → `#pathogen-alerts`
- Medium alert → `#pathogen-monitoring`

### Test 3: Full Integration Test

Run complete test suite:

```bash
python test_slack_integration.py
```

Tests:
1. Connection test
2. Critical severity alert
3. High severity alert
4. Medium severity alert
5. Daily summary
6. NotificationRouter integration

### Manual Testing

#### Test SlackClient Directly

```python
from surveillance.alerting.slack_client import SlackClient
from datetime import datetime

client = SlackClient()

# Test connection
success = client.test_connection()
print(f"Connection: {success}")

# Send test alert
detection = {
    'detection_id': 'manual-test-001',
    'virus_type': 'spumavirus',
    'severity': 'high',
    'copies_per_ml': 250.0,
    'sample_id': 'TEST-001',
    'run_id': 'RUN-001',
    'timestamp': datetime.now().isoformat(),
    'reason': 'Manual test',
    'source': 'internal_pipeline',
    'validation_required': ['Test validation']
}

result = client.send_alert(detection, '#pathogen-alerts')
print(result)
```

#### Test NotificationRouter

```python
from surveillance.alerting.notification_router import NotificationRouter

router = NotificationRouter()

detection = {
    # ... same as above ...
    'notification_config': {
        'slack': ['#pathogen-alerts']
    }
}

status = router.route_alert(detection)
print(status)
```

---

## Deployment

### 1. Deploy Lambda Functions

Ensure Slack client is included in Lambda package:

```bash
cd surveillance/
zip -r lambda-package.zip \
    alerting/slack_client.py \
    alerting/notification_router.py \
    config/ \
    -x "*.pyc" "__pycache__/*"

# Upload to Lambda
aws lambda update-function-code \
    --function-name surveillance-external-collector \
    --zip-file fileb://lambda-package.zip \
    --region ap-northeast-1
```

### 2. Set Environment Variables

Run setup script:

```bash
./scripts/setup_lambda_env.sh
```

### 3. Test in AWS

Trigger a test event:

```bash
aws lambda invoke \
    --function-name surveillance-pipeline-listener \
    --payload file://test-event.json \
    --region ap-northeast-1 \
    response.json
```

### 4. Monitor Logs

Check CloudWatch Logs:

```bash
aws logs tail /aws/lambda/surveillance-external-collector \
    --follow \
    --region ap-northeast-1
```

Look for:
```
Initialized Notification Router (Slack: enabled)
Slack notification sent to #critical-alerts: sent
```

---

## Message Format Examples

### Critical Alert

```
🚨 CRITICAL ALERT
SPUMAVIRUS detected in MinION surveillance system

Severity: CRITICAL
Viral Load: 650.50 copies/mL
Sample ID: SAMPLE-001
Run ID: RUN-001
Source: internal_pipeline_kraken2
Timestamp: 2025-11-15 10:30:45

Reason:
High viral load - xenotransplantation risk

Validation Required:
• Sequence confirmation (Sanger/consensus)
• Independent bioinformatics verification
• MHLW/PMDA notification preparation

[View Dashboard] [Acknowledge]

MinION 4-Virus Surveillance System | Detection ID: det-12345
```

### High Alert

```
⚠️ HIGH PRIORITY ALERT
HANTAVIRUS detected in MinION surveillance system

Severity: HIGH
Viral Load: 150.00 copies/mL
...
```

### Daily Summary

```
📊 Daily Surveillance Summary - 2025-11-15

Samples Processed: 15
Total Detections: 3
MAFF Reports: 2 new
Academic Papers: 1 new

Alert Summary:
🚨 Critical: 1
⚠️ High: 1
ℹ️ Medium: 1
📝 Low: 0

MinION 4-Virus Surveillance System
```

---

## Troubleshooting

### Issue: "Slack client not enabled"

**Cause**: Missing environment variables

**Solution**:
```bash
# Check variables
echo $SLACK_BOT_TOKEN
echo $SLACK_WEBHOOK_URL

# Set from .env
export $(cat .env | grep -v '^#' | xargs)

# Or set manually
export SLACK_BOT_TOKEN=xoxb-your-token
```

### Issue: "channel_not_found"

**Cause**: Bot not invited to channel or channel doesn't exist

**Solution**:
1. Invite bot to channel: `/invite @4-Virus Surveillance`
2. Or use `chat:write.public` scope (already configured)
3. Verify channel name starts with `#`

### Issue: "invalid_auth"

**Cause**: Expired or incorrect bot token

**Solution**:
1. Go to https://api.slack.com/apps/A09TVLTGDSL
2. Reinstall app to workspace
3. Copy new bot token
4. Update environment variables

### Issue: "not_allowed_token_type"

**Cause**: Using user token instead of bot token

**Solution**: Ensure token starts with `xoxb-` (not `xoxp-`)

### Issue: Messages sent but not visible

**Cause**: Bot not in channel

**Solution**:
- With `chat:write.public` scope, bot can post without joining
- If not working, manually invite bot: `/invite @4-Virus Surveillance`

### Issue: Rate limiting

**Cause**: Too many messages (Slack allows 1 message/second)

**Solution**:
- System implements automatic deduplication (1 hour window)
- Rate limiting handled in `severity_engine.py`
- Check deduplication rules in `config/severity_rules.yaml`

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)

from surveillance.alerting.slack_client import SlackClient
client = SlackClient()
client.test_connection()
```

Check Lambda logs:

```bash
aws logs tail /aws/lambda/surveillance-external-collector \
    --since 10m \
    --filter-pattern "Slack" \
    --region ap-northeast-1
```

---

## Security Best Practices

### 1. Never Commit Secrets

Add to `.gitignore`:

```
.env
.env.local
*.secret
```

### 2. Use Environment Variables

Always load from environment, never hardcode:

```python
# ✅ Good
bot_token = os.environ.get('SLACK_BOT_TOKEN')

# ❌ Bad
bot_token = 'xoxb-hardcoded-token'
```

### 3. Rotate Credentials

Periodically rotate Slack credentials:
1. Generate new bot token in Slack app settings
2. Update all Lambda functions
3. Update local `.env` file
4. Invalidate old token

### 4. Least Privilege

Only grant required scopes:
- `chat:write` - Required
- `chat:write.public` - Recommended
- `channels:read` - Optional

Avoid:
- `admin` scopes
- `files:write` (not needed)
- User-level permissions

---

## Performance & Costs

### Slack API Limits

- **Tier 1**: 1 message/second per channel
- **Burst**: Up to 100 messages in bursts
- **Daily**: Unlimited (within rate limits)

### System Implementation

- Deduplication: 1-hour window prevents duplicate alerts
- Rate limiting: Handled by severity engine
- Queueing: No additional queue needed (low volume)

### Expected Usage

| Scenario | Messages/Day | Cost |
|----------|--------------|------|
| Normal operations | 10-50 | Free |
| Outbreak detection | 100-200 | Free |
| Full surveillance | 500+ | Free |

Slack free tier: **Unlimited messages**

---

## Advanced Configuration

### Custom Channel Routing

Edit `config/severity_rules.yaml`:

```yaml
alert_routing:
  critical:
    slack:
      - "#critical-alerts"
      - "#lab-director"  # Additional channel

  high:
    slack:
      - "#pathogen-alerts"
```

### Custom Message Templates

Edit `surveillance/alerting/slack_client.py:_build_alert_message()`:

```python
# Add custom blocks
blocks.append({
    "type": "section",
    "text": {
        "type": "mrkdwn",
        "text": "*Custom Field:*\nYour custom data"
    }
})
```

### Webhook-Only Mode

If Bot API unavailable, use webhook only:

```bash
# Don't set SLACK_BOT_TOKEN
export SLACK_WEBHOOK_URL=https://hooks.slack.com/services/...
```

Note: Webhooks have limited features (no buttons, single channel)

---

## Support

### Documentation

- [Slack API Documentation](https://api.slack.com/docs)
- [Block Kit Builder](https://app.slack.com/block-kit-builder)
- [Bot Token Scopes](https://api.slack.com/scopes)

### Internal Resources

- Main README: `../README.md`
- Architecture: `ARCHITECTURE.md`
- Quick Reference: `QUICK_REFERENCE.md`

### Contact

For issues with Slack integration:
1. Check CloudWatch logs
2. Run test suite: `python tests/test_slack_integration.py`
3. Verify credentials in `.env`
4. Check this guide's troubleshooting section

---

## Changelog

### v2.2.0 (2025-11-15)

- Initial Slack integration implementation
- Bot API with webhook fallback
- Severity-based channel routing
- Rich Block Kit message formatting
- Daily summary notifications
- Comprehensive test suite
- Setup automation scripts

---

**Last Updated**: 2025-11-15
**Version**: 2.2.0
**Author**: MinION Surveillance Team
