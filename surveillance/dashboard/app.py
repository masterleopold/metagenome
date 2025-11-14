"""
Real-time 4-Virus Surveillance Dashboard

Streamlit-based dashboard for monitoring virus detections from both
external sources and internal MinION pipeline.

Run with: streamlit run surveillance/dashboard/app.py
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
from pathlib import Path
import sys
import json

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent.parent))

import boto3
from boto3.dynamodb.conditions import Key, Attr


# Page configuration
st.set_page_config(
    page_title="4-Virus Surveillance Dashboard",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# AWS Configuration
AWS_REGION = "ap-northeast-1"
DYNAMODB_DETECTIONS_TABLE = "surveillance-detections"
DYNAMODB_EXTERNAL_TABLE = "surveillance-external-updates"


@st.cache_resource
def get_dynamodb_client():
    """Get cached DynamoDB client"""
    return boto3.resource('dynamodb', region_name=AWS_REGION)


def load_recent_detections(hours=24):
    """
    Load recent virus detections from DynamoDB

    Args:
        hours: Number of hours to look back

    Returns:
        DataFrame of detections
    """
    try:
        dynamodb = get_dynamodb_client()
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        # Calculate timestamp threshold
        threshold_time = datetime.now() - timedelta(hours=hours)
        threshold_timestamp = int(threshold_time.timestamp())

        # Scan table for recent detections
        response = table.scan(
            FilterExpression=Attr('timestamp_sort').gte(threshold_timestamp)
        )

        detections = response.get('Items', [])

        if not detections:
            return pd.DataFrame()

        df = pd.DataFrame(detections)
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        df = df.sort_values('timestamp', ascending=False)

        return df

    except Exception as e:
        st.error(f"Failed to load detections: {e}")
        return pd.DataFrame()


def load_external_updates(days=7):
    """
    Load external source updates

    Args:
        days: Number of days to look back

    Returns:
        DataFrame of external updates
    """
    try:
        dynamodb = get_dynamodb_client()
        table = dynamodb.Table(DYNAMODB_EXTERNAL_TABLE)

        # Scan for recent updates
        response = table.scan()
        updates = response.get('Items', [])

        if not updates:
            return pd.DataFrame()

        df = pd.DataFrame(updates)
        df['update_date'] = pd.to_datetime(df['update_date'])
        df = df[df['update_date'] >= datetime.now() - timedelta(days=days)]
        df = df.sort_values('update_date', ascending=False)

        return df

    except Exception as e:
        st.error(f"Failed to load external updates: {e}")
        return pd.DataFrame()


def main():
    """Main dashboard application"""

    # Header
    st.title("🔬 4-Virus Surveillance Dashboard")
    st.markdown("Real-time monitoring of Hantavirus, Polyomavirus, Spumavirus, and EEEV")

    # Sidebar configuration
    st.sidebar.header("Dashboard Settings")

    time_range = st.sidebar.selectbox(
        "Time Range",
        ["Last 24 hours", "Last 7 days", "Last 30 days", "All time"],
        index=0
    )

    auto_refresh = st.sidebar.checkbox("Auto-refresh (30s)", value=True)

    if auto_refresh:
        import time
        time.sleep(30)
        st.rerun()

    # Map time range to hours
    hours_map = {
        "Last 24 hours": 24,
        "Last 7 days": 24 * 7,
        "Last 30 days": 24 * 30,
        "All time": 24 * 365
    }
    hours = hours_map[time_range]

    # Load data
    with st.spinner("Loading surveillance data..."):
        detections_df = load_recent_detections(hours=hours)
        external_df = load_external_updates(days=hours//24)

    # Metrics row
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        total_detections = len(detections_df)
        st.metric("Total Detections", total_detections)

    with col2:
        if not detections_df.empty:
            critical_count = len(detections_df[detections_df['severity'] == 'critical'])
            st.metric("Critical Alerts", critical_count, delta_color="inverse")
        else:
            st.metric("Critical Alerts", 0)

    with col3:
        if not detections_df.empty:
            high_count = len(detections_df[detections_df['severity'] == 'high'])
            st.metric("High Priority", high_count)
        else:
            st.metric("High Priority", 0)

    with col4:
        external_count = len(external_df)
        st.metric("External Sources Updated", external_count)

    st.divider()

    # Tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs([
        "🚨 Active Alerts",
        "📊 Analytics",
        "🌐 External Sources",
        "📈 Trends"
    ])

    # Tab 1: Active Alerts
    with tab1:
        st.header("Active Detections")

        if detections_df.empty:
            st.success("✅ No detections in selected time range")
        else:
            # Filter by severity
            severity_filter = st.multiselect(
                "Filter by Severity",
                options=['critical', 'high', 'medium', 'low'],
                default=['critical', 'high']
            )

            filtered_df = detections_df[detections_df['severity'].isin(severity_filter)]

            # Display detections
            for _, detection in filtered_df.iterrows():
                severity = detection.get('severity', 'unknown')

                # Color coding
                severity_colors = {
                    'critical': '🔴',
                    'high': '🟠',
                    'medium': '🟡',
                    'low': '🟢'
                }

                with st.expander(
                    f"{severity_colors.get(severity, '⚪')} {detection['virus_type'].upper()} - "
                    f"{detection.get('timestamp', 'Unknown time')} - {severity.upper()}",
                    expanded=(severity == 'critical')
                ):
                    col_a, col_b = st.columns(2)

                    with col_a:
                        st.write("**Detection Details**")
                        st.write(f"- **Virus**: {detection.get('virus_type', 'Unknown')}")
                        st.write(f"- **Source**: {detection.get('source', 'Unknown')}")
                        st.write(f"- **Severity**: {severity.upper()}")
                        st.write(f"- **Sample ID**: {detection.get('sample_id', 'N/A')}")
                        st.write(f"- **Run ID**: {detection.get('run_id', 'N/A')}")

                    with col_b:
                        st.write("**Quantification**")
                        copies = detection.get('copies_per_ml', detection.get('reads', 0))
                        st.write(f"- **Estimated Copies/mL**: {copies:,.0f}")
                        st.write(f"- **Reads**: {detection.get('reads', 'N/A')}")
                        st.write(f"- **Percentage**: {detection.get('percentage', 'N/A')}")

                    st.write("**Reason**")
                    st.info(detection.get('reason', 'No reason provided'))

                    if detection.get('metadata'):
                        st.write("**Metadata**")
                        st.json(detection['metadata'])

    # Tab 2: Analytics
    with tab2:
        st.header("Detection Analytics")

        if not detections_df.empty:
            col_chart1, col_chart2 = st.columns(2)

            with col_chart1:
                st.subheader("Detections by Virus Type")
                virus_counts = detections_df['virus_type'].value_counts()
                fig_virus = px.pie(
                    values=virus_counts.values,
                    names=virus_counts.index,
                    title="Distribution by Virus"
                )
                st.plotly_chart(fig_virus, use_container_width=True)

            with col_chart2:
                st.subheader("Detections by Severity")
                severity_counts = detections_df['severity'].value_counts()
                colors = {'critical': '#DC3545', 'high': '#FFC107', 'medium': '#17A2B8', 'low': '#28A745'}
                fig_severity = px.bar(
                    x=severity_counts.index,
                    y=severity_counts.values,
                    title="Severity Distribution",
                    color=severity_counts.index,
                    color_discrete_map=colors
                )
                st.plotly_chart(fig_severity, use_container_width=True)

            # Timeline
            st.subheader("Detection Timeline")
            detections_df['date'] = detections_df['timestamp'].dt.date
            timeline_data = detections_df.groupby(['date', 'virus_type']).size().reset_index(name='count')

            fig_timeline = px.line(
                timeline_data,
                x='date',
                y='count',
                color='virus_type',
                title="Detections Over Time",
                markers=True
            )
            st.plotly_chart(fig_timeline, use_container_width=True)

        else:
            st.info("No detections to analyze in selected time range")

    # Tab 3: External Sources
    with tab3:
        st.header("External Information Sources")

        if not external_df.empty:
            # Group by source
            maff_df = external_df[external_df['source'] == 'maff']
            estat_df = external_df[external_df['source'] == 'estat']
            academic_df = external_df[external_df['source'] == 'academic']

            col_ext1, col_ext2, col_ext3 = st.columns(3)

            with col_ext1:
                st.subheader("🏛️ MAFF Updates")
                if not maff_df.empty:
                    latest_maff = maff_df.iloc[0]
                    st.metric("Latest Check", latest_maff['update_date'].strftime('%Y-%m-%d'))
                    st.metric("New Reports", latest_maff.get('new_reports', 0))
                    st.metric("Keywords Found", len(latest_maff.get('virus_keywords_found', [])))

                    if latest_maff.get('virus_keywords_found'):
                        st.warning(f"Keywords: {', '.join(latest_maff['virus_keywords_found'])}")
                else:
                    st.info("No MAFF updates")

            with col_ext2:
                st.subheader("📊 E-Stat Updates")
                if not estat_df.empty:
                    latest_estat = estat_df.iloc[0]
                    st.metric("Latest Check", latest_estat['update_date'].strftime('%Y-%m-%d'))
                    st.metric("Tables Checked", latest_estat.get('stats_checked', 0))
                    st.metric("Keywords Found", len(latest_estat.get('virus_keywords_found', [])))
                else:
                    st.info("No E-Stat updates")

            with col_ext3:
                st.subheader("📚 Academic Papers")
                if not academic_df.empty:
                    latest_academic = academic_df.iloc[0]
                    st.metric("Latest Check", latest_academic['update_date'].strftime('%Y-%m-%d'))
                    st.metric("New Publications", latest_academic.get('total_articles', 0))

                    if latest_academic.get('articles_by_virus'):
                        st.write("**By Virus:**")
                        for virus, count in latest_academic['articles_by_virus'].items():
                            st.write(f"- {virus}: {count}")
                else:
                    st.info("No academic updates")

        else:
            st.info("No external source updates in selected time range")

    # Tab 4: Trends
    with tab4:
        st.header("Surveillance Trends")

        if not detections_df.empty:
            # Add trend analysis
            st.subheader("Detection Rate Trends")

            detections_df['hour'] = detections_df['timestamp'].dt.floor('H')
            hourly_counts = detections_df.groupby('hour').size().reset_index(name='count')

            fig_trend = go.Figure()
            fig_trend.add_trace(go.Scatter(
                x=hourly_counts['hour'],
                y=hourly_counts['count'],
                mode='lines+markers',
                name='Detections per Hour',
                line=dict(color='#17A2B8', width=2)
            ))

            fig_trend.update_layout(
                title="Hourly Detection Rate",
                xaxis_title="Time",
                yaxis_title="Detections",
                hovermode='x unified'
            )

            st.plotly_chart(fig_trend, use_container_width=True)

            # Severity trends
            st.subheader("Severity Trends")
            severity_over_time = detections_df.groupby([
                pd.Grouper(key='timestamp', freq='D'),
                'severity'
            ]).size().reset_index(name='count')

            fig_severity_trend = px.area(
                severity_over_time,
                x='timestamp',
                y='count',
                color='severity',
                title="Severity Trends Over Time"
            )

            st.plotly_chart(fig_severity_trend, use_container_width=True)

        else:
            st.info("Insufficient data for trend analysis")

    # Footer
    st.divider()
    st.caption(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | "
               f"Auto-refresh: {'ON' if auto_refresh else 'OFF'}")


if __name__ == "__main__":
    main()
