#!/bin/bash
##############################################################################
# J-STAGE Data Cleanup Script
#
# Purpose: Delete J-STAGE academic data older than 24 hours from S3
#          to comply with J-STAGE Terms of Service Article 3, Clause 5
#
# Usage:
#   ./cleanup_jstage_data.sh [--dry-run] [--bucket BUCKET_NAME]
#
# Options:
#   --dry-run       Show what would be deleted without actually deleting
#   --bucket NAME   Specify S3 bucket (default: surveillance-data)
#   --region REGION AWS region (default: ap-northeast-1)
#
# Author: Yoichiro Hara
# Date: 2025-01-17
##############################################################################

set -euo pipefail

# Default values
BUCKET="surveillance-data"
REGION="ap-northeast-1"
DRY_RUN=false
PREFIX="external/academic/"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --bucket)
      BUCKET="$2"
      shift 2
      ;;
    --region)
      REGION="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [--dry-run] [--bucket BUCKET_NAME] [--region REGION]"
      echo ""
      echo "Delete J-STAGE academic data older than 24 hours from S3"
      echo ""
      echo "Options:"
      echo "  --dry-run       Show what would be deleted without actually deleting"
      echo "  --bucket NAME   Specify S3 bucket (default: surveillance-data)"
      echo "  --region REGION AWS region (default: ap-northeast-1)"
      echo "  -h, --help      Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check AWS CLI is installed
if ! command -v aws &> /dev/null; then
    echo -e "${RED}Error: AWS CLI is not installed${NC}"
    exit 1
fi

# Check AWS credentials
if ! aws sts get-caller-identity --region "$REGION" &> /dev/null; then
    echo -e "${RED}Error: AWS credentials not configured${NC}"
    exit 1
fi

echo "================================================================"
echo "J-STAGE Data Cleanup - ToS Compliance Tool"
echo "================================================================"
echo "Bucket:    $BUCKET"
echo "Region:    $REGION"
echo "Prefix:    $PREFIX"
echo "Dry Run:   $DRY_RUN"
echo "================================================================"
echo ""

# Calculate cutoff time (24 hours ago)
CUTOFF=$(date -u -d '24 hours ago' +%Y-%m-%dT%H:%M:%S 2>/dev/null || date -u -v-24H +%Y-%m-%dT%H:%M:%S)
echo "Cutoff time: $CUTOFF (24 hours ago)"
echo ""

# List objects older than 24 hours
echo "Scanning for J-STAGE data older than 24 hours..."
OLD_OBJECTS=$(aws s3api list-objects-v2 \
  --bucket "$BUCKET" \
  --prefix "$PREFIX" \
  --region "$REGION" \
  --query "Contents[?LastModified<'$CUTOFF'].Key" \
  --output text 2>/dev/null || echo "")

if [ -z "$OLD_OBJECTS" ]; then
  echo -e "${GREEN}✅ No old J-STAGE data found. System is compliant.${NC}"
  exit 0
fi

# Count objects
OBJECT_COUNT=$(echo "$OLD_OBJECTS" | wc -l | tr -d ' ')
echo -e "${YELLOW}Found $OBJECT_COUNT objects older than 24 hours:${NC}"
echo ""

# Display objects to be deleted
echo "$OLD_OBJECTS" | while read -r key; do
  if [ -n "$key" ]; then
    # Get object details
    LAST_MODIFIED=$(aws s3api head-object \
      --bucket "$BUCKET" \
      --key "$key" \
      --region "$REGION" \
      --query 'LastModified' \
      --output text 2>/dev/null || echo "Unknown")

    SIZE=$(aws s3api head-object \
      --bucket "$BUCKET" \
      --key "$key" \
      --region "$REGION" \
      --query 'ContentLength' \
      --output text 2>/dev/null || echo "0")

    SIZE_KB=$((SIZE / 1024))

    echo "  - $key"
    echo "    Last Modified: $LAST_MODIFIED"
    echo "    Size: ${SIZE_KB} KB"
    echo ""
  fi
done

# Confirm deletion
if [ "$DRY_RUN" = true ]; then
  echo -e "${YELLOW}DRY RUN: No files will be deleted.${NC}"
  exit 0
fi

echo -e "${RED}WARNING: This will permanently delete $OBJECT_COUNT objects from S3.${NC}"
read -p "Do you want to continue? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
  echo "Deletion cancelled."
  exit 0
fi

# Delete objects
echo ""
echo "Deleting old J-STAGE data..."
DELETED_COUNT=0
FAILED_COUNT=0

echo "$OLD_OBJECTS" | while read -r key; do
  if [ -n "$key" ]; then
    if aws s3 rm "s3://$BUCKET/$key" --region "$REGION" > /dev/null 2>&1; then
      echo -e "${GREEN}✓${NC} Deleted: $key"
      ((DELETED_COUNT++))
    else
      echo -e "${RED}✗${NC} Failed to delete: $key"
      ((FAILED_COUNT++))
    fi
  fi
done

echo ""
echo "================================================================"
echo "Cleanup Summary"
echo "================================================================"
echo "Total objects found:    $OBJECT_COUNT"
echo "Successfully deleted:   $DELETED_COUNT"
echo "Failed deletions:       $FAILED_COUNT"
echo "================================================================"

if [ $FAILED_COUNT -gt 0 ]; then
  echo -e "${RED}⚠️  Some deletions failed. Please check AWS permissions.${NC}"
  exit 1
else
  echo -e "${GREEN}✅ Cleanup completed successfully. System is now compliant.${NC}"
  exit 0
fi
