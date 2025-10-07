#!/usr/bin/env python3
"""
Terminate EC2 instances after phase completion
"""

import json
import boto3
import os
from typing import Dict, Any, List
from datetime import datetime, timedelta

ec2 = boto3.client('ec2')
sns = boto3.client('sns')

SNS_TOPIC = os.environ['SNS_TOPIC_ARN']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Terminate EC2 instances based on various criteria.

    Can be called:
    1. After phase completion
    2. On workflow failure
    3. Periodically to clean up orphaned instances
    """

    action = event.get('action', 'terminate_phase')

    if action == 'terminate_phase':
        # Terminate specific instance after phase completion
        instance_id = event.get('instance_id')
        phase = event.get('phase')
        return terminate_instance(instance_id, phase)

    elif action == 'terminate_workflow':
        # Terminate all instances for a workflow
        run_id = event.get('run_id')
        return terminate_workflow_instances(run_id)

    elif action == 'cleanup_orphaned':
        # Clean up instances running longer than expected
        max_age_hours = event.get('max_age_hours', 24)
        return cleanup_orphaned_instances(max_age_hours)

    else:
        return {
            'status': 'ERROR',
            'message': f'Unknown action: {action}'
        }

def terminate_instance(instance_id: str, phase: str) -> Dict[str, Any]:
    """Terminate a specific instance."""

    try:
        # Check instance exists and is running
        response = ec2.describe_instances(InstanceIds=[instance_id])
        if not response['Reservations']:
            return {
                'status': 'NOT_FOUND',
                'instance_id': instance_id
            }

        instance = response['Reservations'][0]['Instances'][0]
        current_state = instance['State']['Name']

        if current_state in ['terminated', 'terminating']:
            return {
                'status': 'ALREADY_TERMINATED',
                'instance_id': instance_id
            }

        # Terminate the instance
        ec2.terminate_instances(InstanceIds=[instance_id])

        print(f'Terminated instance {instance_id} for phase {phase}')

        return {
            'status': 'TERMINATED',
            'instance_id': instance_id,
            'phase': phase,
            'previous_state': current_state
        }

    except Exception as e:
        print(f'Error terminating instance {instance_id}: {e}')
        return {
            'status': 'ERROR',
            'instance_id': instance_id,
            'error': str(e)
        }

def terminate_workflow_instances(run_id: str) -> Dict[str, Any]:
    """Terminate all instances for a workflow."""

    # Find instances with RunID tag
    response = ec2.describe_instances(
        Filters=[
            {'Name': 'tag:RunID', 'Values': [run_id]},
            {'Name': 'instance-state-name', 'Values': ['running', 'pending', 'stopping', 'stopped']}
        ]
    )

    instance_ids = []
    for reservation in response['Reservations']:
        for instance in reservation['Instances']:
            instance_ids.append(instance['InstanceId'])

    if not instance_ids:
        return {
            'status': 'NO_INSTANCES',
            'run_id': run_id
        }

    # Terminate all instances
    ec2.terminate_instances(InstanceIds=instance_ids)

    print(f'Terminated {len(instance_ids)} instances for run {run_id}')

    # Send notification
    send_notification(
        f'Instances Terminated - Run {run_id}',
        f'Terminated {len(instance_ids)} instances for workflow {run_id}:\n' +
        '\n'.join(instance_ids)
    )

    return {
        'status': 'TERMINATED',
        'run_id': run_id,
        'instance_count': len(instance_ids),
        'instance_ids': instance_ids
    }

def cleanup_orphaned_instances(max_age_hours: int) -> Dict[str, Any]:
    """Clean up instances running longer than expected."""

    cutoff_time = datetime.utcnow() - timedelta(hours=max_age_hours)

    # Find old instances with AutoTerminate tag
    response = ec2.describe_instances(
        Filters=[
            {'Name': 'tag:AutoTerminate', 'Values': ['true']},
            {'Name': 'instance-state-name', 'Values': ['running']}
        ]
    )

    orphaned = []
    for reservation in response['Reservations']:
        for instance in reservation['Instances']:
            launch_time = instance['LaunchTime'].replace(tzinfo=None)
            if launch_time < cutoff_time:
                orphaned.append({
                    'instance_id': instance['InstanceId'],
                    'launch_time': launch_time.isoformat(),
                    'age_hours': (datetime.utcnow() - launch_time).total_seconds() / 3600,
                    'name': get_tag_value(instance['Tags'], 'Name'),
                    'run_id': get_tag_value(instance['Tags'], 'RunID')
                })

    if not orphaned:
        return {
            'status': 'NO_ORPHANED',
            'max_age_hours': max_age_hours
        }

    # Terminate orphaned instances
    instance_ids = [i['instance_id'] for i in orphaned]
    ec2.terminate_instances(InstanceIds=instance_ids)

    print(f'Terminated {len(orphaned)} orphaned instances')

    # Send alert
    orphan_details = '\n'.join([
        f"{i['instance_id']} ({i['name']}) - Age: {i['age_hours']:.1f}h, Run: {i['run_id']}"
        for i in orphaned
    ])

    send_notification(
        f'Orphaned Instances Terminated',
        f'Terminated {len(orphaned)} instances older than {max_age_hours} hours:\n\n' +
        orphan_details
    )

    return {
        'status': 'TERMINATED',
        'orphaned_count': len(orphaned),
        'max_age_hours': max_age_hours,
        'instances': orphaned
    }

def get_tag_value(tags: List[Dict], key: str) -> str:
    """Get tag value from tag list."""
    for tag in tags:
        if tag['Key'] == key:
            return tag['Value']
    return ''

def send_notification(subject: str, message: str):
    """Send SNS notification."""
    sns.publish(
        TopicArn=SNS_TOPIC,
        Subject=subject,
        Message=message
    )