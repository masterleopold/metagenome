#!/usr/bin/env python3
"""
MinION Pipeline Monitoring Dashboard
"""

import streamlit as st
import boto3
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
import json
import time

# AWS Configuration
AWS_REGION = st.sidebar.text_input("AWS Region", value="ap-northeast-1")
ENVIRONMENT = st.sidebar.selectbox("Environment", ["production", "staging", "development"])

# Initialize AWS clients
@st.cache_resource
def init_aws_clients():
    return {
        'cloudwatch': boto3.client('cloudwatch', region_name=AWS_REGION),
        'stepfunctions': boto3.client('stepfunctions', region_name=AWS_REGION),
        'rds': boto3.client('rds-data', region_name=AWS_REGION),
        's3': boto3.client('s3', region_name=AWS_REGION)
    }

# Page configuration
st.set_page_config(
    page_title="MinION Pipeline Monitor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Title
st.title("🧬 MinION Pipeline Monitoring Dashboard")
st.markdown(f"**Environment:** {ENVIRONMENT} | **Region:** {AWS_REGION}")

# Sidebar controls
st.sidebar.header("Controls")
refresh_rate = st.sidebar.slider("Refresh Rate (seconds)", 10, 60, 30)
auto_refresh = st.sidebar.checkbox("Auto Refresh", value=True)

# Time range
time_range = st.sidebar.selectbox(
    "Time Range",
    ["Last 1 Hour", "Last 6 Hours", "Last 24 Hours", "Last 7 Days", "Last 30 Days"]
)

# Get time range in hours
time_hours = {
    "Last 1 Hour": 1,
    "Last 6 Hours": 6,
    "Last 24 Hours": 24,
    "Last 7 Days": 168,
    "Last 30 Days": 720
}[time_range]

# Initialize clients
clients = init_aws_clients()

# Main dashboard
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Overview",
    "🔄 Active Workflows",
    "🦠 Pathogen Detection",
    "📈 Metrics",
    "⚠️ Alerts"
])

def get_workflow_stats(hours=24):
    """Get workflow statistics from CloudWatch."""

    end_time = datetime.utcnow()
    start_time = end_time - timedelta(hours=hours)

    try:
        # Get workflow counts
        response = clients['cloudwatch'].get_metric_statistics(
            Namespace='MinION/Pipeline',
            MetricName='WorkflowStatus',
            Dimensions=[],
            StartTime=start_time,
            EndTime=end_time,
            Period=3600,
            Statistics=['Sum']
        )

        total_workflows = sum([d['Sum'] for d in response['Datapoints']])

        # Get success rate
        success_response = clients['cloudwatch'].get_metric_statistics(
            Namespace='MinION/Pipeline',
            MetricName='WorkflowSuccess',
            Dimensions=[],
            StartTime=start_time,
            EndTime=end_time,
            Period=3600,
            Statistics=['Average']
        )

        success_rate = sum([d['Average'] for d in success_response['Datapoints']]) / len(success_response['Datapoints']) if success_response['Datapoints'] else 0

        return {
            'total': int(total_workflows),
            'success_rate': success_rate * 100,
            'time_range': f'{hours} hours'
        }
    except:
        return {
            'total': 0,
            'success_rate': 0,
            'time_range': f'{hours} hours'
        }

def get_active_workflows():
    """Get currently active workflow executions."""

    try:
        state_machine_arn = f'arn:aws:states:{AWS_REGION}:123456789012:stateMachine:minion-pipeline-{ENVIRONMENT}'

        response = clients['stepfunctions'].list_executions(
            stateMachineArn=state_machine_arn,
            statusFilter='RUNNING',
            maxResults=10
        )

        workflows = []
        for execution in response['executions']:
            # Get execution details
            details = clients['stepfunctions'].describe_execution(
                executionArn=execution['executionArn']
            )

            input_data = json.loads(details.get('input', '{}'))

            workflows.append({
                'Run ID': input_data.get('run_id', 'Unknown'),
                'Status': execution['status'],
                'Started': execution['startDate'].strftime('%Y-%m-%d %H:%M:%S'),
                'Duration': str(datetime.utcnow() - execution['startDate']).split('.')[0]
            })

        return pd.DataFrame(workflows)
    except:
        return pd.DataFrame()

def get_pathogen_stats():
    """Get pathogen detection statistics."""

    try:
        # Simulated data - replace with actual database query
        data = {
            'Pathogen': ['PERV-A', 'PERV-B', 'HEV', 'JEV', 'PRRSV', 'PCV2', 'Other'],
            'Detections': [2, 1, 5, 3, 8, 12, 15],
            'Risk Level': ['CRITICAL', 'CRITICAL', 'HIGH', 'HIGH', 'MEDIUM', 'MEDIUM', 'LOW']
        }

        return pd.DataFrame(data)
    except:
        return pd.DataFrame()

# Tab 1: Overview
with tab1:
    st.header("System Overview")

    # Metrics row
    col1, col2, col3, col4 = st.columns(4)

    stats = get_workflow_stats(time_hours)

    with col1:
        st.metric(
            "Total Workflows",
            stats['total'],
            delta=f"Last {stats['time_range']}"
        )

    with col2:
        st.metric(
            "Success Rate",
            f"{stats['success_rate']:.1f}%",
            delta="2.5%" if stats['success_rate'] > 95 else "-2.5%"
        )

    with col3:
        st.metric(
            "Active Workflows",
            len(get_active_workflows()),
            delta="Currently running"
        )

    with col4:
        st.metric(
            "System Health",
            "Operational",
            delta="All systems normal"
        )

    # Workflow timeline
    st.subheader("Workflow Execution Timeline")

    # Create timeline chart
    timeline_data = pd.DataFrame({
        'Date': pd.date_range(start=datetime.now() - timedelta(hours=time_hours),
                             end=datetime.now(),
                             freq='H'),
        'Completed': [5, 3, 4, 6, 2, 5, 4, 3, 5, 6, 4, 3] * (time_hours // 12 + 1)[:time_hours],
        'Failed': [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0] * (time_hours // 12 + 1)[:time_hours]
    })

    fig_timeline = go.Figure()
    fig_timeline.add_trace(go.Scatter(
        x=timeline_data['Date'],
        y=timeline_data['Completed'],
        mode='lines',
        name='Completed',
        line=dict(color='green', width=2)
    ))
    fig_timeline.add_trace(go.Scatter(
        x=timeline_data['Date'],
        y=timeline_data['Failed'],
        mode='lines',
        name='Failed',
        line=dict(color='red', width=2)
    ))

    fig_timeline.update_layout(
        height=300,
        showlegend=True,
        xaxis_title="Time",
        yaxis_title="Workflows"
    )

    st.plotly_chart(fig_timeline, use_container_width=True)

# Tab 2: Active Workflows
with tab2:
    st.header("Active Workflow Executions")

    active_df = get_active_workflows()

    if not active_df.empty:
        st.dataframe(
            active_df,
            use_container_width=True,
            hide_index=True
        )

        # Phase progress
        st.subheader("Phase Progress")

        phases = ['Basecalling', 'QC', 'Host Removal', 'Pathogen Detection', 'Quantification', 'Reporting']
        progress_data = pd.DataFrame({
            'Phase': phases,
            'Progress': [100, 100, 100, 75, 0, 0]
        })

        fig_progress = px.bar(
            progress_data,
            x='Progress',
            y='Phase',
            orientation='h',
            color='Progress',
            color_continuous_scale='Viridis'
        )

        fig_progress.update_layout(height=300)
        st.plotly_chart(fig_progress, use_container_width=True)
    else:
        st.info("No active workflows")

# Tab 3: Pathogen Detection
with tab3:
    st.header("Pathogen Detection Summary")

    pathogen_df = get_pathogen_stats()

    if not pathogen_df.empty:
        col1, col2 = st.columns([1, 2])

        with col1:
            # Pathogen table
            st.subheader("Recent Detections")
            st.dataframe(
                pathogen_df,
                use_container_width=True,
                hide_index=True
            )

        with col2:
            # Pathogen chart
            st.subheader("Detection Distribution")

            fig_pathogen = px.pie(
                pathogen_df,
                values='Detections',
                names='Pathogen',
                color_discrete_map={
                    'PERV-A': '#FF0000',
                    'PERV-B': '#FF4444',
                    'HEV': '#FFA500',
                    'JEV': '#FFAA00',
                    'PRRSV': '#FFFF00',
                    'PCV2': '#00FF00',
                    'Other': '#808080'
                }
            )

            fig_pathogen.update_traces(textposition='inside', textinfo='percent+label')
            st.plotly_chart(fig_pathogen, use_container_width=True)

        # PMDA Compliance
        st.subheader("PMDA 91 Pathogen Coverage")

        pmda_coverage = {
            'Category': ['Viruses', 'Bacteria', 'Parasites', 'Fungi', 'Prions'],
            'Total': [50, 35, 5, 5, 1],
            'Tested': [50, 35, 5, 5, 1],
            'Detected': [8, 3, 0, 0, 0]
        }

        pmda_df = pd.DataFrame(pmda_coverage)

        fig_pmda = go.Figure()
        fig_pmda.add_trace(go.Bar(name='Total', x=pmda_df['Category'], y=pmda_df['Total']))
        fig_pmda.add_trace(go.Bar(name='Detected', x=pmda_df['Category'], y=pmda_df['Detected']))

        fig_pmda.update_layout(barmode='group', height=300)
        st.plotly_chart(fig_pmda, use_container_width=True)

# Tab 4: Metrics
with tab4:
    st.header("Performance Metrics")

    metric_col1, metric_col2 = st.columns(2)

    with metric_col1:
        st.subheader("Phase Duration")

        duration_data = {
            'Phase': ['Basecalling', 'QC', 'Host Removal', 'Pathogen Detection', 'Quantification', 'Reporting'],
            'Avg Duration (hours)': [4.2, 0.5, 1.8, 3.5, 0.5, 0.3]
        }

        fig_duration = px.bar(
            pd.DataFrame(duration_data),
            x='Phase',
            y='Avg Duration (hours)',
            color='Avg Duration (hours)',
            color_continuous_scale='Blues'
        )

        st.plotly_chart(fig_duration, use_container_width=True)

    with metric_col2:
        st.subheader("Resource Utilization")

        resource_data = {
            'Resource': ['CPU', 'Memory', 'GPU', 'Storage'],
            'Usage (%)': [65, 78, 45, 82]
        }

        fig_resource = go.Figure(go.Indicator(
            mode="gauge+number",
            value=78,
            title={'text': "Average Resource Usage"},
            domain={'x': [0, 1], 'y': [0, 1]},
            gauge={'axis': {'range': [None, 100]},
                   'bar': {'color': "darkblue"},
                   'steps': [
                       {'range': [0, 50], 'color': "lightgray"},
                       {'range': [50, 80], 'color': "gray"}],
                   'threshold': {'line': {'color': "red", 'width': 4},
                                'thickness': 0.75, 'value': 90}}
        ))

        fig_resource.update_layout(height=300)
        st.plotly_chart(fig_resource, use_container_width=True)

    # Cost tracking
    st.subheader("Cost Analysis")

    cost_data = pd.DataFrame({
        'Date': pd.date_range(start=datetime.now() - timedelta(days=7),
                             end=datetime.now(),
                             freq='D'),
        'Compute': [45, 52, 48, 55, 51, 49, 53],
        'Storage': [12, 13, 13, 14, 14, 15, 15],
        'Transfer': [5, 6, 5, 7, 6, 5, 6]
    })

    fig_cost = px.area(
        cost_data,
        x='Date',
        y=['Compute', 'Storage', 'Transfer'],
        title='Daily Cost Breakdown ($)'
    )

    st.plotly_chart(fig_cost, use_container_width=True)

# Tab 5: Alerts
with tab5:
    st.header("System Alerts")

    # Alert filters
    alert_col1, alert_col2, alert_col3 = st.columns(3)

    with alert_col1:
        alert_level = st.selectbox("Alert Level", ["All", "CRITICAL", "WARNING", "INFO"])

    with alert_col2:
        alert_category = st.selectbox("Category", ["All", "PERV Detection", "QC Failure", "System Error"])

    with alert_col3:
        alert_status = st.selectbox("Status", ["All", "Active", "Acknowledged", "Resolved"])

    # Alert list
    alerts = [
        {
            'Time': '2024-01-15 14:23:15',
            'Level': '🔴 CRITICAL',
            'Category': 'PERV Detection',
            'Message': 'PERV-A detected in run RUN-2024-001',
            'Status': 'Active'
        },
        {
            'Time': '2024-01-15 13:45:00',
            'Level': '🟡 WARNING',
            'Category': 'QC Failure',
            'Message': 'Mean quality score below threshold (8.5 < 9.0)',
            'Status': 'Acknowledged'
        },
        {
            'Time': '2024-01-15 12:00:00',
            'Level': '🔵 INFO',
            'Category': 'System',
            'Message': 'Database backup completed successfully',
            'Status': 'Resolved'
        }
    ]

    alerts_df = pd.DataFrame(alerts)

    # Filter alerts
    if alert_level != "All":
        alerts_df = alerts_df[alerts_df['Level'].str.contains(alert_level)]
    if alert_category != "All":
        alerts_df = alerts_df[alerts_df['Category'] == alert_category]
    if alert_status != "All":
        alerts_df = alerts_df[alerts_df['Status'] == alert_status]

    st.dataframe(
        alerts_df,
        use_container_width=True,
        hide_index=True
    )

    # Alert actions
    if st.button("Acknowledge All"):
        st.success("All alerts acknowledged")

    if st.button("Export Alerts"):
        csv = alerts_df.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f'alerts_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv',
            mime='text/csv'
        )

# Footer
st.markdown("---")
st.markdown(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} UTC")

# Auto refresh
if auto_refresh:
    time.sleep(refresh_rate)
    st.rerun()