#!/bin/bash
#
# Setup Lambda Environment Variables for Slack Integration
#
# This script updates Lambda function environment variables with Slack credentials
# Run after deploying Lambda functions to enable Slack notifications
#
# Usage: ./setup_lambda_env.sh

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# AWS Configuration
AWS_REGION="ap-northeast-1"

# Lambda functions to update
LAMBDA_FUNCTIONS=(
    "surveillance-external-collector"
    "surveillance-pipeline-listener"
)

echo -e "${GREEN}=== Lambda Environment Variables Setup ===${NC}"
echo ""

# Check if .env file exists
if [ ! -f "../.env" ]; then
    echo -e "${RED}Error: .env file not found${NC}"
    echo "Please copy .env.template to .env and fill in your Slack credentials"
    exit 1
fi

# Load environment variables from .env
echo "Loading credentials from .env file..."
export $(cat ../.env | grep -v '^#' | xargs)

# Verify required variables
REQUIRED_VARS=(
    "SLACK_BOT_TOKEN"
    "SLACK_APP_ID"
    "SLACK_CLIENT_ID"
    "SLACK_CLIENT_SECRET"
    "SLACK_SIGNING_SECRET"
)

echo "Validating required environment variables..."
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo -e "${RED}Error: $var not set in .env${NC}"
        exit 1
    fi
    echo -e "  ${GREEN}✓${NC} $var"
done

echo ""
echo "Updating Lambda functions..."

# Update each Lambda function
for function_name in "${LAMBDA_FUNCTIONS[@]}"; do
    echo ""
    echo -e "${YELLOW}Updating $function_name...${NC}"

    # Check if function exists
    if ! aws lambda get-function --function-name "$function_name" --region "$AWS_REGION" &>/dev/null; then
        echo -e "${RED}  ✗ Function not found: $function_name${NC}"
        echo "  Skipping..."
        continue
    fi

    # Update environment variables
    aws lambda update-function-configuration \
        --function-name "$function_name" \
        --region "$AWS_REGION" \
        --environment "Variables={
            SLACK_BOT_TOKEN=${SLACK_BOT_TOKEN},
            SLACK_WEBHOOK_URL=${SLACK_WEBHOOK_URL:-},
            SLACK_APP_ID=${SLACK_APP_ID},
            SLACK_CLIENT_ID=${SLACK_CLIENT_ID},
            SLACK_CLIENT_SECRET=${SLACK_CLIENT_SECRET},
            SLACK_SIGNING_SECRET=${SLACK_SIGNING_SECRET},
            AWS_REGION=${AWS_REGION},
            DYNAMODB_DETECTIONS_TABLE=${DYNAMODB_DETECTIONS_TABLE:-surveillance-detections},
            DYNAMODB_NOTIFICATIONS_TABLE=${DYNAMODB_NOTIFICATIONS_TABLE:-surveillance-notifications},
            ENABLE_SLACK_NOTIFICATIONS=${ENABLE_SLACK_NOTIFICATIONS:-true}
        }" \
        --output json > /dev/null

    if [ $? -eq 0 ]; then
        echo -e "  ${GREEN}✓${NC} Updated successfully"
    else
        echo -e "  ${RED}✗${NC} Update failed"
    fi
done

echo ""
echo -e "${GREEN}=== Setup Complete ===${NC}"
echo ""
echo "Next steps:"
echo "1. Test Slack notifications with: python ../tests/test_slack_integration.py"
echo "2. Verify channel permissions in Slack workspace"
echo "3. Monitor CloudWatch logs for notification delivery"
echo ""
