# MinION Pipeline API Documentation

## Base URL

```
https://api.minion-pipeline.com/{environment}/
```

Environments:
- `production`: Production environment
- `staging`: Staging environment
- `development`: Development environment

## Authentication

All API requests require an API key in the header:

```
x-api-key: YOUR_API_KEY
```

## Endpoints

### 1. Start Workflow

**POST** `/workflows`

Start a new analysis workflow.

#### Request Body

```json
{
  "run_id": "RUN-2024-001",
  "bucket": "minion-data",
  "input_prefix": "runs/RUN-2024-001/fast5/",
  "config": {
    "phases": {
      "basecalling": {
        "enabled": true,
        "skip_duplex": false
      },
      "qc": {
        "enabled": true,
        "min_quality": 9
      },
      "host_removal": {
        "enabled": true,
        "reference": "sus_scrofa_11.1"
      },
      "pathogen_detection": {
        "enabled": true,
        "databases": ["kraken2", "rvdb", "pmda"]
      },
      "quantification": {
        "enabled": true,
        "spike_in": "PhiX174"
      },
      "reporting": {
        "enabled": true,
        "formats": ["pdf", "json"]
      }
    }
  }
}
```

#### Response

```json
{
  "run_id": "RUN-2024-001",
  "workflow_id": "123456",
  "execution_arn": "arn:aws:states:...",
  "status": "STARTED",
  "estimated_duration_hours": 8.5
}
```

#### Status Codes

- `200`: Workflow started successfully
- `400`: Invalid input parameters
- `401`: Authentication failed
- `429`: Rate limit exceeded
- `500`: Internal server error

---

### 2. Get Workflow Status

**GET** `/workflows/{workflow_id}`

Get the status of a workflow execution.

#### Path Parameters

- `workflow_id`: Unique workflow identifier

#### Response

```json
{
  "workflow_id": "123456",
  "run_id": "RUN-2024-001",
  "status": "RUNNING",
  "created_at": "2024-01-15T10:00:00Z",
  "updated_at": "2024-01-15T14:30:00Z",
  "phases": [
    {
      "name": "basecalling",
      "status": "COMPLETED",
      "started_at": "2024-01-15T10:05:00Z",
      "completed_at": "2024-01-15T12:00:00Z"
    },
    {
      "name": "qc",
      "status": "COMPLETED",
      "started_at": "2024-01-15T12:05:00Z",
      "completed_at": "2024-01-15T12:30:00Z"
    },
    {
      "name": "host_removal",
      "status": "RUNNING",
      "started_at": "2024-01-15T12:35:00Z"
    }
  ],
  "metrics": {
    "basecalling": {
      "total_reads": 50000,
      "mean_quality": 10.5,
      "total_bases": 150000000
    }
  },
  "pathogens": []
}
```

#### Status Values

- `INITIATED`: Workflow created but not started
- `RUNNING`: Currently executing
- `COMPLETED`: Successfully completed
- `FAILED`: Execution failed
- `CANCELLED`: Manually cancelled

---

### 3. List Workflows

**GET** `/workflows`

List workflow executions with optional filtering.

#### Query Parameters

- `status`: Filter by status (RUNNING, COMPLETED, FAILED)
- `run_id`: Filter by run ID (partial match)
- `limit`: Maximum results (default: 20, max: 100)
- `offset`: Pagination offset
- `days`: Get workflows from last N days (default: 7)

#### Response

```json
{
  "workflows": [
    {
      "workflow_id": "123456",
      "run_id": "RUN-2024-001",
      "status": "COMPLETED",
      "created_at": "2024-01-15T10:00:00Z",
      "updated_at": "2024-01-15T18:00:00Z",
      "phases_completed": 6,
      "error_count": 0,
      "has_pathogens": true,
      "duration_seconds": 28800
    }
  ],
  "pagination": {
    "total": 45,
    "limit": 20,
    "offset": 0,
    "has_more": true
  },
  "filters": {
    "status": null,
    "run_id": null,
    "days": 7
  }
}
```

---

### 4. Get Pathogen Results

**GET** `/pathogens/{run_id}`

Get pathogen detection results for a specific run.

#### Path Parameters

- `run_id`: Run identifier

#### Query Parameters

- `pmda_only`: Show only PMDA 91 pathogens (default: false)
- `min_reads`: Minimum read count threshold (default: 0)
- `include_quantification`: Include absolute quantification data (default: false)

#### Response

```json
{
  "run_id": "RUN-2024-001",
  "workflow_id": "123456",
  "total_pathogens": 8,
  "pmda_pathogens_detected": 5,
  "pmda_compliance": {
    "total": 91,
    "detected": 5,
    "not_detected": 86,
    "detection_rate": 0.055,
    "detected_list": ["HEV", "PCV2", "SS", "EC", "CA"],
    "critical_detected": false
  },
  "perv_analysis": {
    "detected": false
  },
  "pathogens": [
    {
      "code": "PCV2",
      "name": "Porcine circovirus 2",
      "read_count": 1250,
      "confidence": 0.98,
      "method": "kraken2",
      "is_pmda": true,
      "risk_level": "MEDIUM",
      "quantification": {
        "copies_per_ml": 5.2e4,
        "log10_copies": 4.72,
        "ci_lower": 4.8e4,
        "ci_upper": 5.6e4
      }
    }
  ]
}
```

---

### 5. Stop Workflow

**DELETE** `/workflows/{workflow_id}`

Stop a running workflow execution.

#### Path Parameters

- `workflow_id`: Workflow to stop

#### Request Body (optional)

```json
{
  "reason": "User requested termination"
}
```

#### Response

```json
{
  "workflow_id": "123456",
  "status": "STOPPING",
  "message": "Workflow termination initiated"
}
```

---

### 6. Get Report

**GET** `/reports/{run_id}`

Get analysis report URLs.

#### Path Parameters

- `run_id`: Run identifier

#### Query Parameters

- `format`: Report format (pdf, json, html)

#### Response

```json
{
  "run_id": "RUN-2024-001",
  "reports": {
    "pdf": "https://s3.amazonaws.com/...",
    "json": "https://s3.amazonaws.com/...",
    "html": "https://s3.amazonaws.com/...",
    "checklist": "https://s3.amazonaws.com/..."
  },
  "expires_at": "2024-01-22T10:00:00Z"
}
```

---

## Error Responses

All error responses follow this format:

```json
{
  "error": {
    "code": "VALIDATION_ERROR",
    "message": "Invalid input parameters",
    "details": {
      "field": "run_id",
      "issue": "Required field missing"
    }
  },
  "timestamp": "2024-01-15T10:00:00Z",
  "request_id": "abc123"
}
```

### Error Codes

- `VALIDATION_ERROR`: Input validation failed
- `AUTHENTICATION_ERROR`: Authentication failed
- `AUTHORIZATION_ERROR`: Not authorized for resource
- `NOT_FOUND`: Resource not found
- `RATE_LIMIT_EXCEEDED`: Too many requests
- `INTERNAL_ERROR`: Server error

---

## Rate Limiting

- **Default limit**: 100 requests per minute
- **Burst limit**: 200 requests
- Headers returned:
  - `X-RateLimit-Limit`: Request limit
  - `X-RateLimit-Remaining`: Remaining requests
  - `X-RateLimit-Reset`: Reset timestamp

---

## Webhooks

Configure webhooks to receive notifications:

### Webhook Events

1. **workflow.started**
2. **workflow.completed**
3. **workflow.failed**
4. **pathogen.detected**
5. **perv.detected** (critical alert)

### Webhook Payload

```json
{
  "event": "workflow.completed",
  "timestamp": "2024-01-15T18:00:00Z",
  "data": {
    "run_id": "RUN-2024-001",
    "workflow_id": "123456",
    "status": "COMPLETED",
    "duration_hours": 8.0,
    "pathogens_detected": 5,
    "report_url": "https://..."
  }
}
```

---

## SDK Examples

### Python

```python
import requests

API_KEY = 'your-api-key'
BASE_URL = 'https://api.minion-pipeline.com/production'

headers = {'x-api-key': API_KEY}

# Start workflow
response = requests.post(
    f'{BASE_URL}/workflows',
    headers=headers,
    json={
        'run_id': 'RUN-2024-001',
        'bucket': 'minion-data',
        'input_prefix': 'runs/RUN-2024-001/fast5/'
    }
)

workflow = response.json()
print(f"Workflow started: {workflow['workflow_id']}")

# Check status
status = requests.get(
    f"{BASE_URL}/workflows/{workflow['workflow_id']}",
    headers=headers
).json()

print(f"Status: {status['status']}")
```

### Node.js

```javascript
const axios = require('axios');

const API_KEY = 'your-api-key';
const BASE_URL = 'https://api.minion-pipeline.com/production';

const client = axios.create({
  baseURL: BASE_URL,
  headers: {'x-api-key': API_KEY}
});

// Start workflow
async function startWorkflow(runId, bucket, prefix) {
  const response = await client.post('/workflows', {
    run_id: runId,
    bucket: bucket,
    input_prefix: prefix
  });

  return response.data;
}

// Usage
startWorkflow('RUN-2024-001', 'minion-data', 'runs/RUN-2024-001/fast5/')
  .then(workflow => console.log(`Started: ${workflow.workflow_id}`))
  .catch(error => console.error(error));
```

---

## Support

For API support:
- Email: api-support@minion-pipeline.com
- Documentation: https://docs.minion-pipeline.com
- Status Page: https://status.minion-pipeline.com