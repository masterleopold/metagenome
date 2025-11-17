#!/usr/bin/env python3
"""
J-STAGE Terms of Service Compliance Verification Script

Checks for violations of J-STAGE ToS Article 3, Clause 5:
- S3: No J-STAGE data older than 24 hours
- DynamoDB: All items have TTL set to 24 hours

Usage:
    python verify_jstage_compliance.py [--bucket BUCKET] [--region REGION]

Exit codes:
    0: Compliant
    1: Violations found
    2: Error during check
"""

import argparse
import sys
import boto3
from datetime import datetime, timedelta
from typing import List, Dict, Tuple
from botocore.exceptions import ClientError


class ComplianceViolation:
    """Represents a compliance violation"""

    def __init__(self, source: str, key: str, issue: str, details: str):
        self.source = source  # 'S3' or 'DynamoDB'
        self.key = key
        self.issue = issue
        self.details = details

    def __str__(self):
        return f"[{self.source}] {self.issue}\n  Key: {self.key}\n  Details: {self.details}"


class JStageComplianceChecker:
    """Check J-STAGE Terms of Service compliance"""

    def __init__(self, bucket: str, region: str = 'ap-northeast-1'):
        self.bucket = bucket
        self.region = region
        self.s3_client = boto3.client('s3', region_name=region)
        self.dynamodb_client = boto3.client('dynamodb', region_name=region)
        self.dynamodb = boto3.resource('dynamodb', region_name=region)
        self.violations: List[ComplianceViolation] = []

    def check_s3_lifecycle_rules(self) -> bool:
        """Verify S3 lifecycle rule is configured correctly"""
        print("Checking S3 lifecycle rules...")

        try:
            response = self.s3_client.get_bucket_lifecycle_configuration(
                Bucket=self.bucket
            )

            # Find the academic data rule
            academic_rule = None
            for rule in response.get('Rules', []):
                if 'academic' in rule.get('ID', '').lower():
                    academic_rule = rule
                    break

            if academic_rule is None:
                self.violations.append(ComplianceViolation(
                    source='S3',
                    key='lifecycle-rule',
                    issue='Missing lifecycle rule for academic data',
                    details='No lifecycle rule found with "academic" in the ID'
                ))
                return False

            # Check expiration is set to 1 day
            expiration_days = academic_rule.get('Expiration', {}).get('Days', 0)
            if expiration_days != 1:
                self.violations.append(ComplianceViolation(
                    source='S3',
                    key='lifecycle-rule',
                    issue='Incorrect expiration period',
                    details=f'Expected 1 day, found {expiration_days} days'
                ))
                return False

            # Check rule is enabled
            if academic_rule.get('Status') != 'Enabled':
                self.violations.append(ComplianceViolation(
                    source='S3',
                    key='lifecycle-rule',
                    issue='Lifecycle rule not enabled',
                    details='Rule exists but is not enabled'
                ))
                return False

            print("  ✓ Lifecycle rule configured correctly (1 day expiration)")
            return True

        except ClientError as e:
            if e.response['Error']['Code'] == 'NoSuchLifecycleConfiguration':
                self.violations.append(ComplianceViolation(
                    source='S3',
                    key='lifecycle-rule',
                    issue='No lifecycle configuration found',
                    details='Bucket has no lifecycle rules configured'
                ))
                return False
            else:
                raise

    def check_s3_old_data(self) -> bool:
        """Check for J-STAGE data older than 24 hours in S3"""
        print("Checking S3 for data older than 24 hours...")

        cutoff = datetime.now(datetime.now().astimezone().tzinfo) - timedelta(hours=24)
        prefix = 'external/academic/'
        found_violations = False

        try:
            paginator = self.s3_client.get_paginator('list_objects_v2')
            pages = paginator.paginate(Bucket=self.bucket, Prefix=prefix)

            object_count = 0
            for page in pages:
                for obj in page.get('Contents', []):
                    object_count += 1
                    last_modified = obj['LastModified']

                    if last_modified < cutoff:
                        age_hours = (datetime.now(datetime.now().astimezone().tzinfo) - last_modified).total_seconds() / 3600
                        self.violations.append(ComplianceViolation(
                            source='S3',
                            key=obj['Key'],
                            issue='Data older than 24 hours',
                            details=f'Age: {age_hours:.1f} hours, Size: {obj["Size"]} bytes'
                        ))
                        found_violations = True

            if not found_violations:
                print(f"  ✓ No old data found ({object_count} objects checked)")
                return True
            else:
                print(f"  ✗ Found {len([v for v in self.violations if v.source == 'S3' and 'older than' in v.issue])} objects older than 24 hours")
                return False

        except ClientError as e:
            print(f"  ✗ Error checking S3: {e}")
            return False

    def check_dynamodb_ttl_enabled(self) -> bool:
        """Verify DynamoDB TTL is enabled"""
        print("Checking DynamoDB TTL configuration...")

        table_name = 'surveillance-external-updates'

        try:
            response = self.dynamodb_client.describe_time_to_live(
                TableName=table_name
            )

            ttl_status = response.get('TimeToLiveDescription', {}).get('TimeToLiveStatus')

            if ttl_status != 'ENABLED':
                self.violations.append(ComplianceViolation(
                    source='DynamoDB',
                    key=table_name,
                    issue='TTL not enabled',
                    details=f'TTL status: {ttl_status}'
                ))
                return False

            ttl_attribute = response.get('TimeToLiveDescription', {}).get('AttributeName')
            if ttl_attribute != 'ttl':
                self.violations.append(ComplianceViolation(
                    source='DynamoDB',
                    key=table_name,
                    issue='Incorrect TTL attribute name',
                    details=f'Expected "ttl", found "{ttl_attribute}"'
                ))
                return False

            print(f"  ✓ TTL enabled on attribute '{ttl_attribute}'")
            return True

        except ClientError as e:
            if e.response['Error']['Code'] == 'ResourceNotFoundException':
                self.violations.append(ComplianceViolation(
                    source='DynamoDB',
                    key=table_name,
                    issue='Table not found',
                    details=f'Table {table_name} does not exist'
                ))
                return False
            else:
                raise

    def check_dynamodb_items_have_ttl(self) -> bool:
        """Check that DynamoDB items have TTL set"""
        print("Checking DynamoDB items for TTL values...")

        table_name = 'surveillance-external-updates'

        try:
            table = self.dynamodb.Table(table_name)
            response = table.scan(
                FilterExpression='begins_with(#pk, :academic)',
                ExpressionAttributeNames={'#pk': 'source#date'},
                ExpressionAttributeValues={':academic': 'academic#'},
                Limit=100  # Check recent items only
            )

            items = response.get('Items', [])
            if not items:
                print("  ℹ No academic items found in DynamoDB")
                return True

            items_without_ttl = []
            items_with_expired_ttl = []
            now_timestamp = int(datetime.now().timestamp())

            for item in items:
                item_key = item.get('source#date', 'unknown')

                # Check if TTL field exists
                if 'ttl' not in item:
                    items_without_ttl.append(item_key)
                    self.violations.append(ComplianceViolation(
                        source='DynamoDB',
                        key=item_key,
                        issue='Missing TTL field',
                        details=f'Item created at {item.get("timestamp", "unknown")}'
                    ))
                else:
                    # Check if TTL is set to >24 hours
                    ttl_timestamp = item['ttl']
                    item_timestamp = item.get('timestamp', now_timestamp)

                    expected_ttl = item_timestamp + (24 * 3600)  # 24 hours
                    ttl_diff_hours = (ttl_timestamp - expected_ttl) / 3600

                    if ttl_diff_hours > 1:  # Allow 1 hour tolerance
                        self.violations.append(ComplianceViolation(
                            source='DynamoDB',
                            key=item_key,
                            issue='TTL set to >24 hours',
                            details=f'TTL is {ttl_diff_hours:.1f} hours longer than expected'
                        ))

            if items_without_ttl:
                print(f"  ✗ Found {len(items_without_ttl)} items without TTL")
                return False
            else:
                print(f"  ✓ All {len(items)} items have TTL set correctly")
                return True

        except ClientError as e:
            print(f"  ✗ Error checking DynamoDB items: {e}")
            return False

    def run_all_checks(self) -> Tuple[bool, List[ComplianceViolation]]:
        """Run all compliance checks"""
        print("=" * 70)
        print("J-STAGE Terms of Service Compliance Check")
        print("=" * 70)
        print(f"Bucket: {self.bucket}")
        print(f"Region: {self.region}")
        print(f"Check time: {datetime.now().isoformat()}")
        print("=" * 70)
        print()

        checks = [
            self.check_s3_lifecycle_rules(),
            self.check_s3_old_data(),
            self.check_dynamodb_ttl_enabled(),
            self.check_dynamodb_items_have_ttl()
        ]

        print()
        print("=" * 70)
        print("Compliance Check Summary")
        print("=" * 70)

        all_passed = all(checks)

        if all_passed:
            print("✅ COMPLIANT: All checks passed")
            print()
            print("Your system complies with J-STAGE Terms of Service Article 3,")
            print("Clause 5 (no machine-readable storage >24 hours).")
        else:
            print(f"❌ VIOLATIONS FOUND: {len(self.violations)} compliance issues")
            print()
            for i, violation in enumerate(self.violations, 1):
                print(f"{i}. {violation}")
                print()

            print("=" * 70)
            print("Required Actions:")
            print("=" * 70)
            print("1. Review violations above")
            print("2. Run cleanup script: ./scripts/cleanup_jstage_data.sh")
            print("3. Verify configuration changes are deployed")
            print("4. Re-run this compliance check")

        print("=" * 70)

        return all_passed, self.violations


def main():
    parser = argparse.ArgumentParser(
        description='Verify J-STAGE Terms of Service compliance'
    )
    parser.add_argument(
        '--bucket',
        default='surveillance-data',
        help='S3 bucket name (default: surveillance-data)'
    )
    parser.add_argument(
        '--region',
        default='ap-northeast-1',
        help='AWS region (default: ap-northeast-1)'
    )

    args = parser.parse_args()

    try:
        checker = JStageComplianceChecker(
            bucket=args.bucket,
            region=args.region
        )
        compliant, violations = checker.run_all_checks()

        sys.exit(0 if compliant else 1)

    except Exception as e:
        print(f"\n❌ ERROR: Compliance check failed: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
