# Slack Integration - Quick Start Guide

## What You Have

Slack App Credentials:
```
App ID: A09TVLTGDSL
Client ID: 580133244533.9947707557904
Client Secret: 78bbc76e4c6e2cb7d591d6b4a359d370
Signing Secret: 5ebeb36dadc1df3a68714befc819c615
```

## What You Need

1. **Bot User OAuth Token** (starts with `xoxb-`)
   - Go to: https://api.slack.com/apps/A09TVLTGDSL
   - Navigate to: **OAuth & Permissions**
   - Under **Bot Token Scopes**, ensure these are added:
     - `chat:write`
     - `chat:write.public`
     - `channels:read`
   - Click **"Install to Workspace"** (or **"Reinstall"**)
   - Copy the **Bot User OAuth Token**

2. **Slack Channels** (create these in your workspace):
   - `#critical-alerts` - For CRITICAL severity alerts
   - `#pathogen-alerts` - For HIGH severity alerts
   - `#pathogen-monitoring` - For MEDIUM severity + daily summaries

## Setup Steps

### 1. Configure Environment Variables

```bash
cd surveillance/

# Copy template
cp .env.template .env

# Edit .env and add your Bot Token:
nano .env
```

Update these lines in `.env`:
```bash
SLACK_BOT_TOKEN=xoxb-YOUR-ACTUAL-TOKEN-HERE
SLACK_APP_ID=A09TVLTGDSL
SLACK_CLIENT_ID=580133244533.9947707557904
SLACK_CLIENT_SECRET=78bbc76e4c6e2cb7d591d6b4a359d370
SLACK_SIGNING_SECRET=5ebeb36dadc1df3a68714befc819c615
```

### 2. Load Environment Variables

```bash
export $(cat .env | grep -v '^#' | xargs)
```

### 3. Test Connection

```bash
cd tests/
python test_slack_integration.py --test-conn
```

Expected output:
```
✅ PASSED: Slack connection successful
```

### 4. Test Alert Sending

```bash
python test_slack_integration.py --test-alert
```

This will send test messages to:
- `#critical-alerts` (Critical severity)
- `#pathogen-alerts` (High severity)
- `#pathogen-monitoring` (Medium severity)

### 5. Deploy to AWS Lambda

```bash
cd ../scripts/
./setup_lambda_env.sh
```

This script will:
- Load credentials from `.env`
- Update Lambda function environment variables
- Verify configuration

## Verify Setup

Check that messages appear in your Slack channels with:
- Rich formatting (colored blocks)
- Detection details (virus type, viral load, etc.)
- Action buttons (for critical alerts)

## Troubleshooting

### "Slack client not enabled"
- Check that `SLACK_BOT_TOKEN` is set: `echo $SLACK_BOT_TOKEN`
- Load from .env: `export $(cat ../. env | grep -v '^#' | xargs)`

### "channel_not_found"
- Ensure channels exist in Slack workspace
- Channel names must match exactly: `#critical-alerts`, `#pathogen-alerts`, `#pathogen-monitoring`
- Bot automatically posts to public channels (no need to invite)

### "invalid_auth"
- Bot token may be expired
- Reinstall app in Slack: https://api.slack.com/apps/A09TVLTGDSL → **Install App**
- Copy new token and update `.env`

## Message Examples

### Critical Alert (Spumavirus)
```
🚨 CRITICAL ALERT
SPUMAVIRUS detected in MinION surveillance system

Severity: CRITICAL
Viral Load: 650.50 copies/mL
Sample ID: SAMPLE-001
...

[View Dashboard] [Acknowledge]
```

### High Alert (Hantavirus)
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
Alerts: 🚨1 ⚠️1 ℹ️1
```

## Next Steps

1. ✅ Test locally with `test_slack_integration.py`
2. ✅ Deploy to AWS Lambda with `setup_lambda_env.sh`
3. ✅ Monitor CloudWatch Logs for notification delivery
4. 📖 Read full documentation: `docs/SLACK_SETUP.md`

## Support

- Full Documentation: `surveillance/docs/SLACK_SETUP.md`
- Main README: `surveillance/README.md`
- Slack API Docs: https://api.slack.com/docs

---

**Quick Reference**

```bash
# Test connection
python tests/test_slack_integration.py --test-conn

# Send test alerts
python tests/test_slack_integration.py --test-alert

# Full test suite
python tests/test_slack_integration.py

# Deploy to Lambda
./scripts/setup_lambda_env.sh

# Check logs
aws logs tail /aws/lambda/surveillance-external-collector --follow
```
