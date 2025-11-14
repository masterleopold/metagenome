# Architecture Documentation

## Pipeline Orchestration Pattern

**Lambda-triggered EC2 pattern** (no containers):
1. S3 upload triggers Lambda orchestrator (`lambda/orchestration/pipeline_orchestrator.py`)
2. Lambda launches Step Functions workflow
3. Each phase runs on dedicated EC2 instance with custom AMI
4. EC2 instances use UserData scripts and auto-terminate
5. Phase completion triggers next Lambda → EC2 cycle

## Phase Structure (7-Phase Pipeline)

```
scripts/
├── phase0_sample_prep/     # NEW: Sample routing & workflow determination
│   ├── sample_router.py    # Determines DNA vs RNA extraction workflow
│   └── README.md           # Phase 0 usage documentation
├── phase1_basecalling/     # FAST5→FASTQ (GPU: g4dn.xlarge, Dorado)
├── phase2_qc/              # Quality metrics (t3.large, NanoPlot/PycoQC)
│                           # RNA integrity (RIN) for RNA viruses
├── phase3_host_removal/    # Host depletion (r5.4xlarge, Minimap2)
│                           # CpG methylation (DNA), rRNA depletion (RNA), poly(A) selection
├── phase4_pathogen/        # Multi-DB screening (4 parallel EC2, Kraken2/BLAST)
│   ├── perv_typing.py      # PERV-A/B/C subtype detection
│   ├── detect_recombinants.py
│   ├── perv_phylogenetics.py
│   ├── pmda_targeted_search.py
│   └── detect_pmda_4viruses.py  # NEW: Polyoma/Hantavirus/EEEV/Spumavirus
├── phase5_quantification/  # Abundance calc (t3.large, copies/mL)
└── phase6_reports/         # PMDA-compliant reports (t3.large)
```

## Lambda Functions Organization

```
lambda/
├── orchestration/          # Pipeline coordination
│   ├── pipeline_orchestrator.py  # Main S3 trigger handler
│   └── state_machine_handler.py  # Step Functions logic
├── monitoring/             # EC2 lifecycle and metrics
│   ├── check_phase_status.py
│   ├── collect_metrics.py
│   └── terminate_instances.py
├── api/                    # REST API endpoints
│   ├── get_workflow_status.py
│   └── list_workflows.py
└── shared/                 # Common utilities
    └── input_validation.py # Reusable validation logic
```

## Key Directories

```
tools/              # workflow_cli.py (start/status/metrics), monitoring_dashboard.py
tests/              # pytest suite (uses moto for AWS mocking)
infrastructure/     # Terraform IaC, CloudFormation templates
templates/          # Config YAML, Step Functions JSON, report HTML
docs/               # Technical docs (see TECHNICAL_DETAILS.md, TROUBLESHOOTING.md)
docs-portal/        # Next.js documentation site (port 3003)
md/                 # Japanese strategic documents and protocols
```

## AWS Infrastructure

- **Region**: ap-northeast-1 (Tokyo)
- **Database**: RDS PostgreSQL (Aurora Serverless v2)
- **Shared Storage**: EFS for reference databases (Kraken2, BLAST, PERV)
- **Cost Optimization**: Spot Instances (70% savings), auto-terminating EC2
- **Monitoring**: CloudWatch + SNS alerts, Prometheus metrics

## Data Flow

```
S3 Upload (FAST5) → Lambda Orchestrator → Step Functions → EC2 Phase Execution → Results to S3 → SNS Notifications
```

## EC2 Instance Types by Phase

| Phase | Instance Type | Purpose | Key Features |
|-------|--------------|---------|--------------|
| Phase 1 | g4dn.xlarge | GPU Basecalling | NVIDIA T4 GPU for Dorado |
| Phase 2 | t3.large | QC Analysis | CPU-based quality metrics |
| Phase 3 | r5.4xlarge | Host Removal | High memory for alignment |
| Phase 4 | Multiple t3.xlarge | Pathogen Detection | Parallel processing |
| Phase 5 | t3.large | Quantification | Calculation intensive |
| Phase 6 | t3.large | Report Generation | Document creation |

## Database Structure

### Reference Databases (EFS Mount)
- **Kraken2 DB**: `/mnt/efs/databases/kraken2/pmda_2024/`
- **BLAST DB**: `/mnt/efs/databases/blast/nt_2024/`
- **PERV DB**: `/mnt/efs/databases/perv/references/`
- **Custom PMDA DB**: `/mnt/efs/databases/pmda/2024.1/`

### RDS PostgreSQL Schema
- Pipeline runs metadata
- Sample information
- Pathogen detection results
- Quality metrics
- Report generation logs

## Security Architecture

1. **VPC Configuration**: Private subnets for EC2, public for NAT Gateway
2. **IAM Roles**: Least privilege per Lambda/EC2 function
3. **S3 Encryption**: AES-256 for all data at rest
4. **Secrets Manager**: Database credentials and API keys
5. **CloudTrail**: Audit logging for all API calls

## Cost Optimization Strategies

1. **Spot Instances**: 70% cost reduction for non-critical phases
2. **Auto-termination**: EC2 instances terminate after task completion
3. **EFS Intelligent Tiering**: Automatic cost optimization for reference DBs
4. **Lambda Concurrency Limits**: Prevent runaway costs
5. **S3 Lifecycle Policies**: Archive old run data to Glacier

## Monitoring and Alerting

### CloudWatch Metrics
- Pipeline execution time per phase
- EC2 instance utilization
- Failed run counts
- Cost per analysis

### SNS Alert Topics
- PERV detection (CRITICAL)
- Pipeline failures
- Cost threshold exceeded ($400)
- Quality metrics below threshold

## Scalability Considerations

- **Horizontal Scaling**: Phase 4 runs 4 parallel EC2 instances
- **Queue-based Processing**: SQS for batch sample processing
- **Database Scaling**: Aurora Serverless v2 auto-scales
- **EFS Performance Mode**: Max I/O for high throughput
- **Step Functions**: Handles complex workflow orchestration

---

## 4-Virus Surveillance System Architecture

### Overview

The surveillance system provides **dual-source monitoring** for 4 target viruses in Japanese pig populations:
- **Hantavirus** (ハンタウイルス)
- **Polyomavirus** (ポリオーマウイルス)
- **Spumavirus** (スピューマウイルス) - MHLW Special Management Pathogen #5
- **EEEV** (東部ウマ脳炎ウイルス)

**Architecture Pattern**: Standalone monitoring system with API integration to main pipeline

```
External Sources (Daily)        Internal Pipeline (Real-time)
      ↓                                ↓
  Lambda Collector              S3 Event Trigger
      ↓                                ↓
  DynamoDB  ←─────────────────  Lambda Listener
      ↓
Severity Engine
      ↓
Notification Router
   ├─ SNS/SES (Email/SMS)
   ├─ Streamlit Dashboard
   └─ REST API
```

### Directory Structure

```
surveillance/
├── external/              # External information collection
│   ├── estat_client.py    # E-Stat API (gov statistics)
│   ├── maff_scraper.py    # MAFF surveillance reports
│   └── academic_monitor.py # PubMed + J-STAGE
├── internal/              # Internal pipeline integration
│   └── pipeline_listener.py # Phase 4 result monitoring
├── alerting/              # Notification system
│   ├── severity_engine.py      # 4-level classification
│   └── notification_router.py  # Multi-channel alerts
├── dashboard/             # Streamlit real-time UI
│   └── app.py
├── api/                   # REST API (FastAPI)
│   └── main.py
├── lambda/                # Lambda functions
│   ├── external_collector/ # Daily external data collection
│   └── pipeline_listener/  # Real-time internal monitoring
└── config/
    ├── severity_rules.yaml # Classification rules
    └── config.yaml         # System configuration
```

### Data Flow

#### External Information Collection (Daily @ 11:00 JST)

```
EventBridge Schedule (cron: 0 2 * * ? *)
    ↓
Lambda: external_collector
    ├─ MAFF Scraper (Web scraping)
    ├─ E-Stat API (APP_ID: bae1f981a6d093a9676b03c8eea37324b8de421b)
    └─ Academic Monitor (PubMed API + J-STAGE scraping)
        ↓
    DynamoDB: surveillance-external-updates
    S3: s3://surveillance-data/external/
        ↓
    Keyword Detection → Severity Engine → Alerts
```

#### Internal Pipeline Integration (Real-time)

```
Phase 4 Completion → S3: phase4_results/
    ↓
S3 Event Trigger
    ↓
Lambda: pipeline_listener
    ↓
Parse Kraken2/BLAST results
    ↓
Filter 4 target viruses
    ↓
DynamoDB: surveillance-detections
    ↓
Severity Engine → Notification Router
```

### AWS Infrastructure Components

#### DynamoDB Tables

| Table | Purpose | Key Schema | TTL |
|-------|---------|-----------|-----|
| surveillance-detections | Virus detections | PK: detection_id, SK: timestamp | No |
| surveillance-external-updates | Daily source updates | PK: source#date, SK: update_id | 90d |
| surveillance-notifications | Notification tracking | PK: notification_id, SK: timestamp | 90d |

**Indexes**:
- `virus-type-index`: Query by virus type
- `severity-index`: Query by severity level
- `source-timestamp-index`: Query by source and time

#### S3 Bucket Structure

```
s3://surveillance-data/
├── external/
│   ├── maff/YYYY/MM/DD/reports/*.pdf
│   ├── estat/YYYY/MM/DD/stats.json
│   └── academic/YYYY/MM/DD/papers/*.json
└── internal/
    └── detections/YYYY/MM/DD/HH/detection_*.json
```

**Lifecycle Policies**:
- External data: 365 days retention
- Internal detections: 730 days retention

#### Lambda Functions

| Function | Trigger | Timeout | Memory | Purpose |
|----------|---------|---------|--------|---------|
| external_collector | EventBridge (daily) | 15m | 512MB | Collect from MAFF/E-Stat/Academic |
| pipeline_listener | S3 Event | 5m | 256MB | Monitor Phase 4 results |
| alert_processor | DynamoDB Stream | 3m | 256MB | Process severity & route alerts |

#### SNS Topics

- **4virus-critical-alerts**: CRITICAL severity (Spumavirus >500, ANY EEEV)
- **4virus-high-alerts**: HIGH severity (Hantavirus >100, Polyomavirus >100)
- **4virus-daily-summary**: Daily monitoring summary (09:00 JST)

### Severity Classification System

**4-Level Classification** (defined in `surveillance/config/severity_rules.yaml`):

| Level | Response Time | Notification Channels | Example Criteria |
|-------|---------------|----------------------|------------------|
| CRITICAL | < 5 min | SNS + SMS + Dashboard Flash | Spumavirus >500 copies/mL, ANY EEEV |
| HIGH | < 30 min | SNS + Email + Dashboard | Hantavirus >100, Polyomavirus >100 |
| MEDIUM | < 2 hours | Email + Dashboard | External keyword match |
| LOW | < 24 hours | Dashboard only | Academic publications |

**Engine Features**:
- YAML-based rule configuration
- Compound detection logic (multiple viruses)
- External validation boosting
- Deduplication (1-hour window)

### User Interfaces

#### Streamlit Dashboard

- **Port**: 8501
- **Auto-refresh**: 30 seconds
- **Features**:
  - Real-time detection alerts
  - External source timeline
  - Trend analysis (Plotly charts)
  - Active alert management

**Launch**: `streamlit run surveillance/dashboard/app.py`

#### REST API (FastAPI)

- **Port**: 8000
- **Documentation**: http://localhost:8000/docs

**Key Endpoints**:
```
GET  /api/v1/detections           # List detections (filterable)
GET  /api/v1/alerts/active        # Active alerts summary
GET  /api/v1/external/daily-updates  # External source updates
GET  /api/v1/statistics/trends    # Detection trends
POST /api/v1/webhooks             # External system integration
```

### Integration with Main Pipeline

#### Phase 4 Integration

New detection module: `scripts/phase4_pathogen/detect_4viruses.py`

**Usage**:
```bash
./scripts/phase4_pathogen/detect_4viruses.py \
  --bam results/phase4/RUN-001/aligned.bam \
  --run-id RUN-001 \
  --sample-id SAMPLE-001 \
  --output-dir results/surveillance/
```

**Integration Pattern**: Follows `perv_typing.py` pattern
- BAM-based detection
- Kraken2 report parsing
- SNS immediate alert (on CRITICAL)
- DynamoDB result storage

#### Notification Integration

Reuses existing SNS infrastructure with new topics:
- Same IAM roles pattern as PERV alerts
- Same SES sender configuration
- Additional topic ARNs in `surveillance/config/config.yaml`

### Security Architecture

#### IAM Policies

```
surveillance-external-collector-role:
  - s3:PutObject (surveillance-data/*)
  - dynamodb:PutItem (surveillance-*)
  - sns:Publish (4virus-*)

surveillance-pipeline-listener-role:
  - s3:GetObject (minion-data/results/phase4/*)
  - dynamodb:PutItem, Query (surveillance-*)
  - sns:Publish (4virus-critical)
```

#### Data Encryption

- **S3**: AES-256 server-side encryption
- **DynamoDB**: Encryption at rest enabled
- **SNS/SES**: TLS in transit
- **API Keys**: Stored in `config.yaml` (E-Stat APP_ID)

#### Network Architecture

- **Lambda Functions**: VPC optional (external APIs require internet)
- **Dashboard/API**: Local development or EC2 deployment
- **No public endpoints**: Dashboard accessed via SSH tunnel or VPC

### Monitoring and Observability

#### CloudWatch Logs

- `/aws/lambda/surveillance-external-collector`
- `/aws/lambda/surveillance-pipeline-listener`
- `/aws/lambda/surveillance-alert-processor`

**Retention**: 30 days

#### CloudWatch Metrics

Custom metrics:
- `Surveillance/DailyCollections` (count by source)
- `Surveillance/Detections` (count by virus, severity)
- `Surveillance/AlertsSent` (count by channel)

#### Dashboard Metrics

- Total detections (last 24h/7d/30d)
- Detection rate by virus type
- External source update status
- Alert response times

### Cost Estimates

**Monthly Operating Costs** (estimated):

| Component | Usage | Cost |
|-----------|-------|------|
| Lambda Executions | ~30 daily + events | ~$2 |
| DynamoDB | PAY_PER_REQUEST | ~$5 |
| S3 Storage | ~10GB | ~$0.30 |
| SNS/SES | ~50 notifications/month | ~$0.10 |
| **Total** | | **~$7.40/month** |

*Main pipeline cost: ~$50-200 per run (separate)*

### Scalability Considerations

- **Lambda Concurrent Executions**: 10 reserved (prevent throttling)
- **DynamoDB**: On-demand pricing, auto-scales
- **S3**: Unlimited scalability
- **External API Rate Limits**:
  - E-Stat: 5 req/sec
  - PubMed: 3 req/sec
  - J-STAGE: 1 req/sec (self-imposed)

### Deployment

**Terraform**:
```bash
cd infrastructure/surveillance
terraform init
terraform apply \
  -var="estat_app_id=${E_STAT_APP_ID}" \
  -var="pubmed_email=${PUBMED_EMAIL}"
```

**Manual Setup**:
1. Create DynamoDB tables (see `infrastructure/surveillance/dynamodb_schemas.json`)
2. Create S3 bucket: `surveillance-data`
3. Deploy Lambda functions (see `surveillance/lambda/*/`)
4. Configure EventBridge schedule (daily @ 11:00 JST)
5. Set up SNS topics and subscriptions

### Integration Points with Main Pipeline

| Integration Point | Method | Direction |
|------------------|---------|-----------|
| Phase 4 Detection | S3 Event → Lambda | Pipeline → Surveillance |
| PERV Alert Reuse | SNS Topic ARN | Surveillance → Pipeline SNS |
| Report Generation | API Call | Pipeline → Surveillance API |
| Dashboard Metrics | DynamoDB Query | Dashboard ← Both Systems |

### Related Documentation

- **User Guide**: `/surveillance/README.md`
- **API Documentation**: `/surveillance/api/main.py` (FastAPI auto-docs)
- **Configuration**: `/surveillance/config/config.yaml`
- **Severity Rules**: `/surveillance/config/severity_rules.yaml`
- **Web Portal**: `/docs-portal/src/app/surveillance/page.tsx`