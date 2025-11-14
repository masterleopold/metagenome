"""
REST API for 4-Virus Surveillance System

FastAPI-based API for accessing surveillance data programmatically.
Provides endpoints for detections, alerts, external sources, and statistics.

Run with: uvicorn surveillance.api.main:app --reload --port 8000
"""

from fastapi import FastAPI, HTTPException, Query, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from typing import List, Optional, Dict, Any
from datetime import datetime, timedelta
from pydantic import BaseModel
import boto3
from boto3.dynamodb.conditions import Key, Attr
from decimal import Decimal
import json


# FastAPI app
app = FastAPI(
    title="4-Virus Surveillance API",
    description="Real-time virus surveillance data access for MinION pipeline",
    version="1.0.0"
)

# CORS middleware for dashboard access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# AWS Configuration
AWS_REGION = "ap-northeast-1"
DYNAMODB_DETECTIONS_TABLE = "surveillance-detections"
DYNAMODB_EXTERNAL_TABLE = "surveillance-external-updates"


# Pydantic models
class Detection(BaseModel):
    detection_id: str
    timestamp: str
    virus_type: str
    source: str
    severity: str
    copies_per_ml: Optional[float] = None
    reads: Optional[int] = None
    percentage: Optional[float] = None
    sample_id: Optional[str] = None
    run_id: Optional[str] = None
    reason: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class ExternalUpdate(BaseModel):
    source: str
    update_date: str
    new_items_count: int
    virus_keywords_found: List[str] = []
    errors: List[str] = []


class AlertSummary(BaseModel):
    total_alerts: int
    by_severity: Dict[str, int]
    by_virus: Dict[str, int]
    critical_alerts: List[Detection]


# Dependency for DynamoDB client
def get_dynamodb():
    """Get DynamoDB resource"""
    return boto3.resource('dynamodb', region_name=AWS_REGION)


# Helper to convert DynamoDB Decimal to float
def decimal_to_float(obj):
    """Convert Decimal to float for JSON serialization"""
    if isinstance(obj, list):
        return [decimal_to_float(i) for i in obj]
    elif isinstance(obj, dict):
        return {k: decimal_to_float(v) for k, v in obj.items()}
    elif isinstance(obj, Decimal):
        return float(obj)
    else:
        return obj


# Root endpoint
@app.get("/")
def read_root():
    """API root - health check"""
    return {
        "service": "4-Virus Surveillance API",
        "version": "1.0.0",
        "status": "operational",
        "timestamp": datetime.now().isoformat()
    }


# Health check
@app.get("/health")
def health_check():
    """Health check endpoint"""
    try:
        dynamodb = get_dynamodb()
        # Test DynamoDB connection
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)
        table.table_status

        return {
            "status": "healthy",
            "dynamodb": "connected",
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        return JSONResponse(
            status_code=503,
            content={"status": "unhealthy", "error": str(e)}
        )


# Detections endpoints
@app.get("/api/v1/detections", response_model=List[Detection])
def get_detections(
    virus: Optional[str] = Query(None, description="Filter by virus type"),
    severity: Optional[str] = Query(None, description="Filter by severity"),
    from_date: Optional[str] = Query(None, description="Start date (ISO format)"),
    to_date: Optional[str] = Query(None, description="End date (ISO format)"),
    limit: int = Query(100, le=1000, description="Maximum results"),
    dynamodb = Depends(get_dynamodb)
):
    """
    Get virus detections with optional filters

    - **virus**: Filter by virus type (hantavirus, polyomavirus, spumavirus, eeev)
    - **severity**: Filter by severity (critical, high, medium, low)
    - **from_date**: Start date in ISO format
    - **to_date**: End date in ISO format
    - **limit**: Maximum number of results (default 100, max 1000)
    """
    try:
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        # Build filter expression
        filter_expressions = []

        if virus:
            filter_expressions.append(Attr('virus_type').eq(virus))

        if severity:
            filter_expressions.append(Attr('severity').eq(severity))

        if from_date:
            from_timestamp = int(datetime.fromisoformat(from_date).timestamp())
            filter_expressions.append(Attr('timestamp_sort').gte(from_timestamp))

        if to_date:
            to_timestamp = int(datetime.fromisoformat(to_date).timestamp())
            filter_expressions.append(Attr('timestamp_sort').lte(to_timestamp))

        # Combine filters
        if filter_expressions:
            from functools import reduce
            filter_expr = reduce(lambda a, b: a & b, filter_expressions)
            response = table.scan(FilterExpression=filter_expr, Limit=limit)
        else:
            response = table.scan(Limit=limit)

        items = response.get('Items', [])

        # Convert Decimal to float
        items = decimal_to_float(items)

        # Sort by timestamp descending
        items.sort(key=lambda x: x.get('timestamp', ''), reverse=True)

        return items

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch detections: {str(e)}")


@app.get("/api/v1/detections/{detection_id}", response_model=Detection)
def get_detection_by_id(detection_id: str, dynamodb = Depends(get_dynamodb)):
    """
    Get specific detection by ID

    - **detection_id**: Unique detection identifier
    """
    try:
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        response = table.scan(
            FilterExpression=Attr('detection_id').eq(detection_id)
        )

        items = response.get('Items', [])

        if not items:
            raise HTTPException(status_code=404, detail="Detection not found")

        return decimal_to_float(items[0])

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch detection: {str(e)}")


@app.get("/api/v1/detections/realtime", response_model=List[Detection])
def get_realtime_detections(
    hours: int = Query(24, le=168, description="Hours to look back"),
    dynamodb = Depends(get_dynamodb)
):
    """
    Get real-time detections from internal pipeline

    - **hours**: Number of hours to look back (default 24, max 168)
    """
    try:
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        threshold_time = datetime.now() - timedelta(hours=hours)
        threshold_timestamp = int(threshold_time.timestamp())

        response = table.scan(
            FilterExpression=Attr('timestamp_sort').gte(threshold_timestamp) &
                           Attr('source').begins_with('internal_pipeline')
        )

        items = decimal_to_float(response.get('Items', []))
        items.sort(key=lambda x: x.get('timestamp', ''), reverse=True)

        return items

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch real-time detections: {str(e)}")


# Alerts endpoints
@app.get("/api/v1/alerts/active", response_model=AlertSummary)
def get_active_alerts(
    hours: int = Query(24, description="Hours to look back"),
    dynamodb = Depends(get_dynamodb)
):
    """
    Get active alerts summary

    - **hours**: Number of hours to look back for active alerts
    """
    try:
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        threshold_time = datetime.now() - timedelta(hours=hours)
        threshold_timestamp = int(threshold_time.timestamp())

        response = table.scan(
            FilterExpression=Attr('timestamp_sort').gte(threshold_timestamp)
        )

        items = decimal_to_float(response.get('Items', []))

        # Calculate summary
        by_severity = {}
        by_virus = {}
        critical_alerts = []

        for item in items:
            severity = item.get('severity', 'unknown')
            virus = item.get('virus_type', 'unknown')

            by_severity[severity] = by_severity.get(severity, 0) + 1
            by_virus[virus] = by_virus.get(virus, 0) + 1

            if severity == 'critical':
                critical_alerts.append(item)

        return {
            "total_alerts": len(items),
            "by_severity": by_severity,
            "by_virus": by_virus,
            "critical_alerts": critical_alerts
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch alerts: {str(e)}")


# External sources endpoints
@app.get("/api/v1/external/daily-updates")
def get_daily_external_updates(
    date: Optional[str] = Query(None, description="Date (YYYY-MM-DD)"),
    source: Optional[str] = Query(None, description="Source (maff, estat, academic)"),
    dynamodb = Depends(get_dynamodb)
):
    """
    Get daily external source updates

    - **date**: Specific date in YYYY-MM-DD format (defaults to today)
    - **source**: Filter by source (maff, estat, academic)
    """
    try:
        table = dynamodb.Table(DYNAMODB_EXTERNAL_TABLE)

        if date is None:
            date = datetime.now().strftime("%Y-%m-%d")

        # Build filter
        if source:
            filter_expr = Attr('update_date').eq(date) & Attr('source').eq(source)
        else:
            filter_expr = Attr('update_date').eq(date)

        response = table.scan(FilterExpression=filter_expr)

        items = decimal_to_float(response.get('Items', []))

        return {
            "date": date,
            "updates": items,
            "count": len(items)
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch external updates: {str(e)}")


# Statistics endpoints
@app.get("/api/v1/statistics/trends")
def get_trends(
    days: int = Query(30, le=365, description="Days to analyze"),
    virus: Optional[str] = Query(None, description="Filter by virus"),
    dynamodb = Depends(get_dynamodb)
):
    """
    Get detection trends and statistics

    - **days**: Number of days to analyze (default 30, max 365)
    - **virus**: Optional virus type filter
    """
    try:
        table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)

        threshold_time = datetime.now() - timedelta(days=days)
        threshold_timestamp = int(threshold_time.timestamp())

        if virus:
            filter_expr = Attr('timestamp_sort').gte(threshold_timestamp) & Attr('virus_type').eq(virus)
        else:
            filter_expr = Attr('timestamp_sort').gte(threshold_timestamp)

        response = table.scan(FilterExpression=filter_expr)

        items = decimal_to_float(response.get('Items', []))

        # Calculate trends
        daily_counts = {}
        severity_distribution = {}
        virus_distribution = {}

        for item in items:
            # Daily counts
            date = item.get('timestamp', '')[:10]
            daily_counts[date] = daily_counts.get(date, 0) + 1

            # Severity distribution
            severity = item.get('severity', 'unknown')
            severity_distribution[severity] = severity_distribution.get(severity, 0) + 1

            # Virus distribution
            virus_type = item.get('virus_type', 'unknown')
            virus_distribution[virus_type] = virus_distribution.get(virus_type, 0) + 1

        return {
            "period_days": days,
            "total_detections": len(items),
            "daily_counts": daily_counts,
            "severity_distribution": severity_distribution,
            "virus_distribution": virus_distribution,
            "average_per_day": len(items) / days if days > 0 else 0
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to calculate trends: {str(e)}")


@app.get("/api/v1/statistics/summary")
def get_statistics_summary(dynamodb = Depends(get_dynamodb)):
    """
    Get overall surveillance statistics summary
    """
    try:
        detections_table = dynamodb.Table(DYNAMODB_DETECTIONS_TABLE)
        external_table = dynamodb.Table(DYNAMODB_EXTERNAL_TABLE)

        # Get recent detections (last 7 days)
        threshold_time = datetime.now() - timedelta(days=7)
        threshold_timestamp = int(threshold_time.timestamp())

        detections_response = detections_table.scan(
            FilterExpression=Attr('timestamp_sort').gte(threshold_timestamp)
        )

        external_response = external_table.scan()

        detections = decimal_to_float(detections_response.get('Items', []))
        external_updates = decimal_to_float(external_response.get('Items', []))

        # Calculate summary
        return {
            "last_7_days": {
                "total_detections": len(detections),
                "critical_count": sum(1 for d in detections if d.get('severity') == 'critical'),
                "high_count": sum(1 for d in detections if d.get('severity') == 'high'),
                "by_virus": {
                    virus: sum(1 for d in detections if d.get('virus_type') == virus)
                    for virus in ['hantavirus', 'polyomavirus', 'spumavirus', 'eeev']
                }
            },
            "external_sources": {
                "total_updates": len(external_updates),
                "by_source": {
                    source: sum(1 for u in external_updates if u.get('source') == source)
                    for source in ['maff', 'estat', 'academic']
                }
            },
            "generated_at": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to generate summary: {str(e)}")


# Webhook endpoint for external integrations
@app.post("/api/v1/webhooks")
def receive_webhook(payload: Dict[str, Any]):
    """
    Receive webhook notifications from external systems

    Accepts arbitrary JSON payloads and logs them for processing
    """
    try:
        # Log webhook receipt
        print(f"Webhook received at {datetime.now().isoformat()}: {json.dumps(payload)}")

        # TODO: Process webhook payload
        # This could trigger additional monitoring, alerts, etc.

        return {
            "status": "received",
            "timestamp": datetime.now().isoformat(),
            "payload_keys": list(payload.keys())
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Webhook processing failed: {str(e)}")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
