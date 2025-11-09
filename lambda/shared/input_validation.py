#!/usr/bin/env python3
"""
Input validation utilities to prevent command injection attacks.

All Lambda phase trigger functions MUST use these validation functions
to sanitize user input before passing to bash scripts or EC2 UserData.

Security: OWASP A03:2021 - Injection Prevention
"""

import re
from typing import List


def validate_run_id(run_id: str) -> str:
    """Validate run_id format to prevent command injection.

    Accepts formats like: RUN-2024-001, RUN-20241109-ABC123

    Args:
        run_id: Run identifier from event parameter

    Returns:
        Validated run_id string

    Raises:
        ValueError: If run_id contains invalid characters or format

    Security:
        Prevents: Command injection, path traversal
        Allows: Alphanumeric, dash, underscore only
        Max length: 64 characters
    """
    pattern = r'^[A-Z0-9][A-Z0-9_-]{0,63}$'
    if not re.match(pattern, run_id):
        raise ValueError(
            f"Invalid run_id format: {run_id}. "
            "Must be alphanumeric with dashes/underscores, max 64 chars."
        )
    return run_id


def validate_s3_path_component(component: str, name: str = "path") -> str:
    """Validate S3 path component to prevent command injection.

    Args:
        component: S3 path component (prefix, key, etc.)
        name: Parameter name for error messages

    Returns:
        Validated path component string

    Raises:
        ValueError: If path contains invalid characters or traversal

    Security:
        Prevents: Command injection, path traversal
        Allows: alphanumeric, dash, underscore, forward slash, period
        Blocks: .., backticks, semicolons, pipes, etc.
    """
    pattern = r'^[a-zA-Z0-9/_.-]+$'
    if not re.match(pattern, component):
        raise ValueError(
            f"Invalid {name} format: {component}. "
            "Must contain only alphanumeric, dash, underscore, slash, period."
        )
    if '..' in component:
        raise ValueError(f"Path traversal detected in {name}: {component}")
    return component


def validate_bucket_name(bucket: str) -> str:
    """Validate S3 bucket name format to prevent command injection.

    Args:
        bucket: S3 bucket name

    Returns:
        Validated bucket name string

    Raises:
        ValueError: If bucket name format is invalid

    Security:
        Prevents: Command injection
        Enforces: Official S3 bucket naming rules
        Pattern: lowercase alphanumeric, dots, dashes only
    """
    pattern = r'^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$'
    if not re.match(pattern, bucket):
        raise ValueError(f"Invalid bucket name format: {bucket}")
    return bucket


def validate_database_list(databases: List[str]) -> List[str]:
    """Validate database list to prevent command injection.

    Args:
        databases: List of database names (e.g., ['kraken2', 'rvdb', 'pmda'])

    Returns:
        Validated list of database names

    Raises:
        ValueError: If any database name is invalid

    Security:
        Prevents: Command injection through database names
        Allows: lowercase alphanumeric and underscore only
    """
    pattern = r'^[a-z0-9_]+$'
    for db in databases:
        if not re.match(pattern, db):
            raise ValueError(
                f"Invalid database name: {db}. "
                "Must be lowercase alphanumeric with underscores only."
            )
    return databases


def validate_parameter(value: str, name: str, pattern: str) -> str:
    """Generic parameter validation with custom regex pattern.

    Args:
        value: Parameter value to validate
        name: Parameter name for error messages
        pattern: Regex pattern for validation

    Returns:
        Validated parameter value

    Raises:
        ValueError: If value doesn't match pattern

    Security:
        Generic validation function for custom parameters
        Use specific validators when available
    """
    if not re.match(pattern, value):
        raise ValueError(
            f"Invalid {name} format: {value}. "
            f"Must match pattern: {pattern}"
        )
    return value


def validate_numeric_parameter(value, name: str, min_val: float = None,
                               max_val: float = None) -> float:
    """Validate numeric parameter with optional range checking.

    Args:
        value: Numeric value (int or float)
        name: Parameter name for error messages
        min_val: Minimum allowed value (optional)
        max_val: Maximum allowed value (optional)

    Returns:
        Validated numeric value as float

    Raises:
        ValueError: If value is not numeric or out of range

    Security:
        Prevents: Type confusion attacks
        Enforces: Numeric type and range constraints
    """
    try:
        num_value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid {name}: {value}. Must be numeric.")

    if min_val is not None and num_value < min_val:
        raise ValueError(f"{name} must be >= {min_val}, got {num_value}")

    if max_val is not None and num_value > max_val:
        raise ValueError(f"{name} must be <= {max_val}, got {num_value}")

    return num_value
