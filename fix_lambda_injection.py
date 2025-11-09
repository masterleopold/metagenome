#!/usr/bin/env python3
"""
Automated script to fix command injection vulnerabilities in Lambda phase triggers.

Applies fixes to all remaining Lambda functions:
- trigger_qc.py
- trigger_host_removal.py
- trigger_pathogen_detection.py
- trigger_quantification.py
- trigger_reporting.py
"""

import re
from pathlib import Path


VALIDATION_FUNCTIONS = '''

def validate_run_id(run_id: str) -> str:
    """Validate run_id format to prevent command injection."""
    pattern = r'^[A-Z0-9][A-Z0-9_-]{0,63}$'
    if not re.match(pattern, run_id):
        raise ValueError(f"Invalid run_id format: {run_id}. Must be alphanumeric with dashes/underscores, max 64 chars.")
    return run_id


def validate_s3_path_component(component: str, name: str = "path") -> str:
    """Validate S3 path component to prevent command injection."""
    pattern = r'^[a-zA-Z0-9/_.-]+$'
    if not re.match(pattern, component):
        raise ValueError(f"Invalid {name} format: {component}. Must contain only alphanumeric, dash, underscore, slash, period.")
    if '..' in component:
        raise ValueError(f"Path traversal detected in {name}: {component}")
    return component


def validate_bucket_name(bucket: str) -> str:
    """Validate S3 bucket name format to prevent command injection."""
    pattern = r'^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$'
    if not re.match(pattern, bucket):
        raise ValueError(f"Invalid bucket name format: {bucket}")
    return bucket

'''


def fix_imports(content: str) -> str:
    """Add shlex and re imports if missing."""
    if 'import shlex' not in content:
        content = content.replace('from typing import', 'import re\nimport shlex\nfrom typing import')
    elif 'import re' not in content:
        content = content.replace('import shlex', 'import re\nimport shlex')
    return content


def add_validation_functions(content: str) -> str:
    """Add validation functions after environment variable definitions."""
    if 'def validate_run_id' in content:
        return content  # Already added

    # Find where to insert (after SECURITY_GROUP or similar env var)
    insertion_point = None
    for match in re.finditer(r"(SECURITY_GROUP|SNS_TOPIC|ANALYSIS_AMI)\s*=.*\n", content):
        insertion_point = match.end()

    if insertion_point:
        content = content[:insertion_point] + '\n' + VALIDATION_FUNCTIONS + '\n' + content[insertion_point:]

    return content


def fix_bash_interpolations(content: str) -> str:
    """Wrap bash variable assignments with shlex.quote()."""

    # Pattern 1: export VAR='value'  →  export VAR=shlex.quote(value)
    content = re.sub(
        r"export\s+(\w+)='([^']*\{[^}]+\}[^']*)'",
        lambda m: f"export {m.group(1)}={{shlex.quote(f'{m.group(2)}')",
        content
    )

    # Pattern 2: export VAR='{simple_var}'  →  export VAR=shlex.quote(simple_var)
    patterns = [
        (r"export RUN_ID='\{run_id\}'", "export RUN_ID={shlex.quote(run_id)}"),
        (r"export S3_INPUT='s3://\{bucket\}/\{input_prefix\}'", "export S3_INPUT={shlex.quote(f's3://{bucket}/{input_prefix}')}"),
        (r"export S3_INPUT='s3://\{bucket\}/\{prefix\}'", "export S3_INPUT={shlex.quote(f's3://{bucket}/{prefix}')}"),
        (r"export S3_OUTPUT='s3://\{bucket\}/\{output_prefix\}'", "export S3_OUTPUT={shlex.quote(f's3://{bucket}/{output_prefix}')}"),
        (r"export S3_BASE='s3://\{bucket\}/\{base_prefix\}'", "export S3_BASE={shlex.quote(f's3://{bucket}/{base_prefix}')}"),
        (r"export PHASE='\{phase\}'", "export PHASE={shlex.quote(phase)}"),
    ]

    for pattern, replacement in patterns:
        content = re.sub(pattern, replacement, content)

    # Add comment before first user_data
    content = re.sub(
        r'(\n\s+user_data = f""")',
        r'\n    # Use shlex.quote() to prevent command injection in bash scripts\1',
        content,
        count=1
    )

    # Add comment before first command =
    content = re.sub(
        r'(\n\s+command = f""")',
        r'\n    # Use shlex.quote() to prevent command injection\1',
        content,
        count=1
    )

    return content


def add_input_validation(content: str) -> str:
    """Add input validation in lambda_handler."""

    # Find lambda_handler and add validation
    patterns = [
        (r"(\n\s+run_id = event\['run_id'\])", r"\n    # Validate all inputs to prevent command injection\1\n    run_id = validate_run_id(run_id)"),
        (r"(\n\s+bucket = )([^v][^\n]+)", r"\1validate_bucket_name(\2)"),
        (r"(\n\s+)(prefix|input_prefix|output_prefix|base_prefix)( = [^\n]+)",
         lambda m: f"\n    {m.group(2)} = validate_s3_path_component({m.group(2)}, '{m.group(2)}')  # Added after line"),
    ]

    # Simpler approach: just add validation right after extracting run_id
    if "run_id = event['run_id']" in content and "validate_run_id" not in content.split("def lambda_handler")[1].split("def ")[0]:
        content = content.replace(
            "run_id = event['run_id']",
            "# Validate all inputs to prevent command injection\n    run_id = validate_run_id(event['run_id'])"
        )

    return content


def fix_file(filepath: Path):
    """Apply all fixes to a single file."""
    print(f"Fixing {filepath.name}...")

    content = filepath.read_text()
    original_content = content

    # Apply fixes
    content = fix_imports(content)
    content = add_validation_functions(content)
    content = add_input_validation(content)
    content = fix_bash_interpolations(content)

    if content != original_content:
        filepath.write_text(content)
        print(f"  ✓ {filepath.name} fixed")
    else:
        print(f"  - {filepath.name} no changes needed")


def main():
    """Fix all Lambda trigger files."""
    lambda_dir = Path('/home/user/metagenome/lambda/phases')

    files_to_fix = [
        'trigger_qc.py',
        'trigger_host_removal.py',
        'trigger_pathogen_detection.py',
        'trigger_quantification.py',
        'trigger_reporting.py'
    ]

    print("="*60)
    print("FIXING COMMAND INJECTION VULNERABILITIES")
    print("="*60)

    for filename in files_to_fix:
        filepath = lambda_dir / filename
        if filepath.exists():
            fix_file(filepath)
        else:
            print(f"  ! {filename} not found")

    print("="*60)
    print("All files processed!")


if __name__ == '__main__':
    main()
