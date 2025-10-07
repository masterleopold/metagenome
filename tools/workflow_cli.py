#!/usr/bin/env python3
"""
MinION Pipeline CLI Tool
"""

import click
import boto3
import json
import time
from pathlib import Path
from tabulate import tabulate
from datetime import datetime
import sys
sys.path.append('/opt/minion/lib')

from workflow_manager import WorkflowManager
from config_manager import ConfigManager
from s3_utils import S3Utils
from monitoring_client import MonitoringClient

@click.group()
@click.option('--region', default='ap-northeast-1', help='AWS region')
@click.option('--environment', default='production', help='Environment')
@click.pass_context
def cli(ctx, region, environment):
    """MinION Pipeline Management CLI"""
    ctx.ensure_object(dict)
    ctx.obj['region'] = region
    ctx.obj['environment'] = environment
    ctx.obj['workflow_manager'] = WorkflowManager(region)
    ctx.obj['config_manager'] = ConfigManager()
    ctx.obj['s3_utils'] = S3Utils(region)

@cli.command()
@click.option('--run-id', required=True, help='Run ID for the analysis')
@click.option('--bucket', required=True, help='S3 bucket containing input data')
@click.option('--input-prefix', required=True, help='S3 prefix for input FAST5/POD5 files')
@click.option('--config', type=click.Path(exists=True), help='Configuration file')
@click.option('--skip-validation', is_flag=True, help='Skip input validation')
@click.pass_context
def start(ctx, run_id, bucket, input_prefix, config, skip_validation):
    """Start a new workflow execution"""

    workflow_manager = ctx.obj['workflow_manager']
    config_manager = ctx.obj['config_manager']

    # Load configuration
    if config:
        pipeline_config = config_manager.load_config(config)
    else:
        pipeline_config = config_manager.load_config()

    # Validate configuration
    validation = config_manager.validate_config(pipeline_config)
    if not validation['valid']:
        click.echo(click.style('Configuration errors:', fg='red'))
        for error in validation['errors']:
            click.echo(f'  - {error}')
        return

    if validation['warnings']:
        click.echo(click.style('Configuration warnings:', fg='yellow'))
        for warning in validation['warnings']:
            click.echo(f'  - {warning}')

    # Validate input
    if not skip_validation:
        click.echo('Validating input files...')
        s3_utils = ctx.obj['s3_utils']
        validation = workflow_manager.validate_input(bucket, input_prefix)

        if not validation['valid']:
            click.echo(click.style('Input validation failed:', fg='red'))
            for error in validation['errors']:
                click.echo(f'  - {error}')
            return

        click.echo(f"Found {validation['file_count']} files, "
                  f"total size: {validation['total_size'] / 1e9:.2f} GB")

    # Estimate runtime
    estimates = workflow_manager.estimate_runtime(
        validation['file_count'],
        validation['total_size'],
        config_manager._config_to_dict(pipeline_config)
    )

    click.echo('\nEstimated runtime:')
    for phase, hours in estimates.items():
        if phase != 'total':
            click.echo(f'  {phase}: {hours:.1f} hours')
    click.echo(click.style(f'  Total: {estimates["total"]:.1f} hours', fg='cyan'))

    # Start workflow
    click.echo('\nStarting workflow...')
    state_machine_arn = f'arn:aws:states:{ctx.obj["region"]}:123456789012:stateMachine:minion-pipeline-{ctx.obj["environment"]}'

    execution_arn = workflow_manager.start_workflow(
        run_id=run_id,
        state_machine_arn=state_machine_arn,
        input_bucket=bucket,
        input_prefix=input_prefix,
        config=config_manager._config_to_dict(pipeline_config)
    )

    click.echo(click.style(f'✓ Workflow started successfully', fg='green'))
    click.echo(f'Run ID: {run_id}')
    click.echo(f'Execution ARN: {execution_arn}')

@cli.command()
@click.option('--execution-arn', help='Execution ARN')
@click.option('--run-id', help='Run ID')
@click.option('--watch', is_flag=True, help='Watch status updates')
@click.pass_context
def status(ctx, execution_arn, run_id, watch):
    """Check workflow status"""

    workflow_manager = ctx.obj['workflow_manager']

    if not execution_arn and not run_id:
        click.echo('Please provide either --execution-arn or --run-id')
        return

    if run_id:
        # Find execution by run ID
        workflows = workflow_manager.list_workflows(
            state_machine_arn=f'arn:aws:states:{ctx.obj["region"]}:123456789012:stateMachine:minion-pipeline-{ctx.obj["environment"]}'
        )

        matching = [w for w in workflows if run_id in w['name']]
        if not matching:
            click.echo(f'No workflow found for run ID: {run_id}')
            return

        execution_arn = matching[0]['execution_arn']

    while True:
        # Get status
        status_info = workflow_manager.get_workflow_status(execution_arn)

        # Clear screen if watching
        if watch:
            click.clear()

        # Display status
        click.echo(click.style(f'Workflow Status: {status_info["status"]}',
                              fg='green' if status_info['status'] == 'SUCCEEDED' else
                                 'red' if status_info['status'] == 'FAILED' else
                                 'yellow'))

        click.echo(f'Started: {status_info["started_at"]}')
        if 'stopped_at' in status_info:
            click.echo(f'Stopped: {status_info["stopped_at"]}')
            click.echo(f'Duration: {status_info["duration_seconds"] / 3600:.1f} hours')

        # Get phase history
        phases = workflow_manager.get_phase_history(execution_arn)

        if phases:
            click.echo('\nPhase Progress:')
            phase_table = []
            for phase in phases:
                status_icon = '✓' if phase['status'] == 'COMPLETED' else \
                             '✗' if phase['status'] == 'FAILED' else \
                             '⟳' if phase['status'] == 'RUNNING' else '⋯'

                phase_table.append([
                    status_icon,
                    phase['name'],
                    phase['status'],
                    phase.get('started_at', ''),
                    phase.get('completed_at', phase.get('failed_at', ''))
                ])

            click.echo(tabulate(phase_table,
                               headers=['', 'Phase', 'Status', 'Started', 'Completed'],
                               tablefmt='simple'))

        if not watch or status_info['status'] in ['SUCCEEDED', 'FAILED', 'ABORTED']:
            break

        time.sleep(10)

@cli.command()
@click.option('--state-machine-arn', help='State machine ARN')
@click.option('--status-filter', help='Filter by status')
@click.option('--limit', default=10, help='Maximum results')
@click.pass_context
def list(ctx, state_machine_arn, status_filter, limit):
    """List workflow executions"""

    workflow_manager = ctx.obj['workflow_manager']

    if not state_machine_arn:
        state_machine_arn = f'arn:aws:states:{ctx.obj["region"]}:123456789012:stateMachine:minion-pipeline-{ctx.obj["environment"]}'

    workflows = workflow_manager.list_workflows(
        state_machine_arn=state_machine_arn,
        status_filter=status_filter,
        max_results=limit
    )

    if not workflows:
        click.echo('No workflows found')
        return

    table = []
    for w in workflows:
        run_id = w['name'].split('-')[0] if '-' in w['name'] else w['name']
        status_color = 'green' if w['status'] == 'SUCCEEDED' else \
                       'red' if w['status'] == 'FAILED' else \
                       'yellow'

        table.append([
            run_id,
            click.style(w['status'], fg=status_color),
            w['started_at'],
            w.get('stopped_at', 'Running')
        ])

    click.echo(tabulate(table,
                       headers=['Run ID', 'Status', 'Started', 'Stopped'],
                       tablefmt='grid'))

@cli.command()
@click.option('--execution-arn', required=True, help='Execution ARN')
@click.option('--reason', default='User requested', help='Stop reason')
@click.pass_context
def stop(ctx, execution_arn, reason):
    """Stop a running workflow"""

    workflow_manager = ctx.obj['workflow_manager']

    if click.confirm('Are you sure you want to stop this workflow?'):
        success = workflow_manager.stop_workflow(execution_arn, reason)

        if success:
            click.echo(click.style('✓ Workflow stopped', fg='green'))
        else:
            click.echo(click.style('✗ Failed to stop workflow', fg='red'))

@cli.command()
@click.option('--run-id', required=True, help='Run ID')
@click.option('--phase', help='Specific phase metrics')
@click.pass_context
def metrics(ctx, run_id, phase):
    """View workflow metrics"""

    monitoring = MonitoringClient(ctx.obj['region'])

    summary = monitoring.get_workflow_summary(run_id)

    click.echo(f'\nMetrics for Run: {run_id}')
    click.echo('=' * 50)

    # Phase durations
    if summary['phase_durations']:
        click.echo('\nPhase Durations:')
        for metric in summary['phase_durations']:
            timestamp = datetime.fromisoformat(metric['Timestamp'].replace('Z', '+00:00'))
            click.echo(f'  {timestamp.strftime("%H:%M")} - '
                      f'{metric.get("Average", 0) / 60:.1f} minutes')

    # Pathogen counts
    if summary['pathogens_detected']:
        click.echo(f'\nPathogens Detected: {summary["pathogens_detected"][0].get("Maximum", 0):.0f}')

    # Errors
    click.echo(f'Errors: {summary["error_count"]}')

@cli.command()
@click.option('--source', required=True, help='Source path (local or S3)')
@click.option('--destination', required=True, help='Destination path (local or S3)')
@click.option('--exclude', multiple=True, help='Exclude patterns')
@click.pass_context
def sync(ctx, source, destination, exclude):
    """Sync data between local and S3"""

    s3_utils = ctx.obj['s3_utils']

    click.echo(f'Syncing from {source} to {destination}...')

    result = s3_utils.sync_directories(source, destination, delete=False)

    if 'error' in result:
        click.echo(click.style(f'Error: {result["error"]}', fg='red'))
    else:
        click.echo(click.style('✓ Sync completed', fg='green'))
        if 'copied' in result:
            click.echo(f'  Files copied: {result["copied"]}')
        if 'uploaded' in result:
            click.echo(f'  Files uploaded: {result["uploaded"]}')
        if 'downloaded' in result:
            click.echo(f'  Files downloaded: {result["downloaded"]}')

@cli.command()
@click.option('--create-default', is_flag=True, help='Create default configuration')
@click.option('--validate', type=click.Path(exists=True), help='Validate configuration file')
@click.option('--export', type=click.Path(), help='Export current configuration')
@click.pass_context
def config(ctx, create_default, validate, export):
    """Manage pipeline configuration"""

    config_manager = ctx.obj['config_manager']

    if create_default:
        config_manager.create_default_configs()
        click.echo(click.style('✓ Default configuration created', fg='green'))

    elif validate:
        pipeline_config = config_manager.load_config(validate)
        validation = config_manager.validate_config(pipeline_config)

        if validation['valid']:
            click.echo(click.style('✓ Configuration is valid', fg='green'))
        else:
            click.echo(click.style('✗ Configuration errors:', fg='red'))
            for error in validation['errors']:
                click.echo(f'  - {error}')

        if validation['warnings']:
            click.echo(click.style('⚠ Warnings:', fg='yellow'))
            for warning in validation['warnings']:
                click.echo(f'  - {warning}')

    elif export:
        pipeline_config = config_manager.load_config()
        config_manager.save_config(pipeline_config, export)
        click.echo(f'Configuration exported to {export}')

    else:
        # Show current configuration
        pipeline_config = config_manager.load_config()
        config_dict = config_manager._config_to_dict(pipeline_config)
        click.echo(json.dumps(config_dict, indent=2))

@cli.command()
@click.option('--execution-arn', required=True, help='Execution ARN')
@click.pass_context
def cost(ctx, execution_arn):
    """Calculate workflow cost"""

    workflow_manager = ctx.obj['workflow_manager']

    costs = workflow_manager.get_workflow_cost(execution_arn)

    click.echo('\nWorkflow Cost Breakdown:')
    click.echo('=' * 40)
    click.echo(f'Compute:       ${costs["compute"]:.2f}')
    click.echo(f'Storage:       ${costs["storage"]:.2f}')
    click.echo(f'Data Transfer: ${costs["data_transfer"]:.2f}')
    click.echo('-' * 40)
    click.echo(click.style(f'Total:         ${costs["total"]:.2f}', fg='cyan', bold=True))

if __name__ == '__main__':
    cli()