# Appendix: Architecture Diagrams and Visual Materials

## For NVIDIA Academic Grant Program Application
**Institution:** Meiji University
**Project:** GPU-Accelerated AI Inference for Pathogen Surveillance

---- 

## Diagram 1: Overall System Architecture (DGX Spark + A100 Hybrid)

	┌─────────────────────────────────────────────────────────────────────────────┐
	│                         MEIJI UNIVERSITY DEPLOYMENT                          │
	│                    (On-Premises + Cloud Hybrid Architecture)                 │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  SAMPLE INPUT                                                                 │
	│  ┌──────────────────────────────────────────────────────────────┐           │
	│  │  Oxford Nanopore MinION Sequencer                            │           │
	│  │  • 24-48h sequencing run                                     │           │
	│  │  • Output: 30-50 GB FAST5/POD5 files                        │           │
	│  │  • ~4 million reads per run                                  │           │
	│  └──────────────────────────────────────────────────────────────┘           │
	│                              ▼                                                │
	│  ┌──────────────────────────────────────────────────────────────┐           │
	│  │  Local Storage (Meiji University Network-Attached Storage)   │           │
	│  │  • 2TB NVMe SSD buffer                                       │           │
	│  │  • Encrypted at rest (PMDA compliance)                       │           │
	│  └──────────────────────────────────────────────────────────────┘           │
	└──────────────────────────────────────────────────────────────────────────────┘
	                                    │
	                    ┌───────────────┴────────────────┐
	                    ▼                                 ▼
	┌─────────────────────────────────┐   ┌──────────────────────────────────────┐
	│   PRIMARY PATH                  │   │   BENCHMARKING PATH                  │
	│   DGX Spark (On-Premises)       │   │   DGX Cloud A100 (x86 Baseline)      │
	│   ARM Blackwell Architecture    │   │   For Comparative Study              │
	└─────────────────────────────────┘   └──────────────────────────────────────┘
	              │                                       │
	              ▼                                       ▼
	
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                      PHASE 1: GPU-ACCELERATED BASECALLING                    │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark #1 (Primary)          │  DGX Cloud A100 (Benchmark)              │
	│  ┌─────────────────────────────┐ │  ┌──────────────────────────────────┐   │
	│  │ Hardware:                   │ │  │ Hardware:                        │   │
	│  │ • ARM Cortex (20 cores)     │ │  │ • x86 Xeon (16+ vCPUs)          │   │
	│  │ • Blackwell GPU (6,144 CUDA)│ │  │ • A100 GPU (6,912 CUDA)         │   │
	│  │ • 128GB unified memory      │ │  │ • 80GB HBM2e + 512GB RAM        │   │
	│  │                             │ │  │                                  │   │
	│  │ Software:                   │ │  │ Software:                        │   │
	│  │ • Dorado duplex (ARM build?)│ │  │ • Dorado duplex (x86 build)     │   │
	│  │ • CUDA on ARM               │ │  │ • CUDA on x86                   │   │
	│  │                             │ │  │                                  │   │
	│  │ Input:  30-50 GB FAST5      │ │  │ Input:  30-50 GB FAST5          │   │
	│  │ Output: 0.1-0.5 GB FASTQ    │ │  │ Output: 0.1-0.5 GB FASTQ        │   │
	│  │ Runtime: 4-6 hours (est.)   │ │  │ Runtime: 1-2 hours              │   │
	│  └─────────────────────────────┘ │  └──────────────────────────────────┘   │
	│                                   │                                          │
	│  DGX Spark #2 (Parallel Sample)  │  DGX Cloud A100 (Parallel Burst)        │
	│  • Processes even-numbered runs  │  • 4-8× instances during outbreaks      │
	│  • Same configuration as #1      │  • Elastic scaling                       │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              └───────────────┬───────────────────────┘
	                              ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                      PHASE 2: QUALITY CONTROL (CPU-Based)                    │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark (CPU cores)           │  A100 Instance (CPU)                     │
	│  • NanoPlot QC metrics           │  • Same QC pipeline                      │
	│  • PycoQC visualization          │  • Baseline validation                   │
	│  • Read length/quality stats     │  │                                        │
	│  • Runtime: 20-30 min            │  • Runtime: 20-30 min                    │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              └───────────────┬───────────────────────┘
	                              ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│               PHASE 3: HOST GENOME REMOVAL (Memory + CPU Intensive)          │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark (128GB Unified Memory) │  A100 (512GB System RAM)                │
	│  • Minimap2 alignment            │  • Same Minimap2 pipeline               │
	│  • Sus scrofa genome (3GB)       │  • x86 optimization                     │
	│  • 20 ARM cores (16 threads)     │  • 16 x86 vCPUs                         │
	│  • Runtime: 2.5-3.5 hours (est.) │  • Runtime: 1.8 hours                   │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              └───────────────┬───────────────────────┘
	                              ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│            PHASE 4: PATHOGEN DETECTION (CRITICAL - Memory Bound)             │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark (128GB Unified)       │  A100 (80GB GPU + 512GB RAM)            │
	│  ┌─────────────────────────────┐ │  ┌──────────────────────────────────┐   │
	│  │ Kraken2 Database:           │ │  │ Kraken2 Database:                │   │
	│  │ • PMDA 2024.1 (20-30GB)     │ │  │ • Same database                  │   │
	│  │ • Loaded into 128GB pool    │ │  │ • Loaded into system RAM         │   │
	│  │ • ARM compilation           │ │  │ • x86 native binary              │   │
	│  │                             │ │  │                                  │   │
	│  │ BLAST RVDB:                 │ │  │ BLAST RVDB:                      │   │
	│  │ • Viral database (10-15GB)  │ │  │ • Same database                  │   │
	│  │                             │ │  │                                  │   │
	│  │ PERV Typing:                │ │  │ PERV Typing:                     │   │
	│  │ • Python script (ARM)       │ │  │ • Python script (x86)            │   │
	│  │ • SNS alerts                │ │  │ • SNS alerts                     │   │
	│  │                             │ │  │                                  │   │
	│  │ Runtime: 4-5 hours (est.)   │ │  │ Runtime: 2.5-3 hours             │   │
	│  └─────────────────────────────┘ │  └──────────────────────────────────┘   │
	│                                   │                                          │
	│  **91 PMDA Pathogens Detected**  │  **Same 91 Pathogens**                  │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              └───────────────┬───────────────────────┘
	                              ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│           PHASE 5: AI RISK PREDICTION (GPU Inference - NEW!)                 │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark Blackwell GPU         │  A100 GPU (Training + Inference)        │
	│  ┌─────────────────────────────┐ │  ┌──────────────────────────────────┐   │
	│  │ Transformer Model:          │ │  │ Model Training:                  │   │
	│  │ • BERT-based (768-dim)      │ │  │ • 2-4× A100 distributed         │   │
	│  │ • Pathogen risk scoring     │ │  │ • 10,000+ training samples      │   │
	│  │ • Inference: <10 sec        │ │  │ • 500 GPU hours                 │   │
	│  │                             │ │  │                                  │   │
	│  │ Output:                     │ │  │ Inference (validation):         │   │
	│  │ • Risk Level: LOW/MED/HIGH  │ │  │ • Baseline performance          │   │
	│  │ • Confidence: 87.3%         │ │  │ • x86 optimization              │   │
	│  │ • Recommendation            │ │  │                                  │   │
	│  └─────────────────────────────┘ │  └──────────────────────────────────┘   │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              └───────────────┬───────────────────────┘
	                              ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                    PHASE 6: REPORT GENERATION (CPU)                          │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark (CPU cores)           │  A100 Instance (CPU)                     │
	│  • PDF report (PMDA format)      │  • Same report generation               │
	│  • JSON output (API)             │  • Format validation                    │
	│  • HTML dashboard                │  │                                        │
	│  • Runtime: 15-20 min            │  • Runtime: 15-20 min                    │
	└─────────────────────────────────────────────────────────────────────────────┘
	              │                                       │
	              ▼                                       ▼
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                           FINAL OUTPUT & STORAGE                             │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  Meiji University Archive        │  Cloud Backup (Optional)                 │
	│  • NAS: 7-year retention         │  • Encrypted S3 bucket                  │
	│  • Encrypted at rest             │  • De-identified only                   │
	│  • PMDA audit logs               │  • Disaster recovery                    │
	│  • No internet upload required   │  │                                        │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                         PERFORMANCE COMPARISON                                │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  Metric                     │  DGX Spark (ARM)  │  A100 Cloud (x86)         │
	│─────────────────────────────┼───────────────────┼───────────────────────────│
	│  Total Pipeline Runtime     │  10-13 hours      │  6-7 hours                │
	│  Phase 1 (Basecalling)      │  4-6 hours (est.) │  1-2 hours                │
	│  Phase 4 (Pathogen)         │  4-5 hours (est.) │  2.5-3 hours              │
	│  Power Consumption          │  240W × 2 = 480W  │  500W+ per instance       │
	│  Data Sovereignty           │  ✅ On-premises   │  ⚠️ Cloud upload          │
	│  PMDA Compliance            │  ✅ Compliant     │  ⚠️ Requires approval     │
	│  Parallel Capacity          │  2 samples/batch  │  4-8 samples (elastic)    │
	│  Cost Model                 │  Capital (1-time) │  Operational (per-run)    │
	│  Expected Lifespan          │  5+ years         │  Credit exhaustion        │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	┌─────────────────────────────────────────────────────────────────────────────┐
	│                      RESEARCH CONTRIBUTION MATRIX                             │
	├─────────────────────────────────────────────────────────────────────────────┤
	│  DGX Spark Path                  │  A100 Cloud Path                         │
	├──────────────────────────────────┼──────────────────────────────────────────┤
	│  ✅ Novel ARM genomics benchmark │  ✅ x86 baseline (proven platform)      │
	│  ✅ On-premises deployment demo  │  ✅ Cloud-native best practices          │
	│  ✅ PMDA compliance validation   │  ✅ Elastic scaling demonstration        │
	│  ✅ Power efficiency study       │  ✅ Multi-GPU training (AI model)        │
	│  ✅ ARM CUDA maturity assessment │  ✅ Burst capacity during outbreaks      │
	│  ✅ Student training platform    │  ✅ Rapid iteration for optimization     │
	└─────────────────────────────────────────────────────────────────────────────┘

---- 

## Diagram 2: DGX Spark Deployment at Meiji University

	┌─────────────────────────────────────────────────────────────────────────────┐
	│              MEIJI UNIVERSITY BIOINFORMATICS CORE FACILITY                   │
	│                        (Ikuta Campus, Kawasaki, Japan)                       │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  SECURE BSL-2 LABORATORY                                                      │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  WET LAB AREA                                                          │  │
	│  │  ┌──────────────────────┐    ┌──────────────────────┐                 │  │
	│  │  │ Oxford Nanopore      │    │ Sample Preparation   │                 │  │
	│  │  │ MinION Sequencer     │───▶│ • DNA/RNA Extraction │                 │  │
	│  │  │ • GridION (5 devices)│    │ • Library Prep       │                 │  │
	│  │  │ • PromethION (1)     │    │ • Quality Control    │                 │  │
	│  │  └──────────────────────┘    └──────────────────────┘                 │  │
	│  │              │                                                          │  │
	│  │              ▼ (USB/Network Transfer)                                  │  │
	│  │  ┌──────────────────────────────────────────────────────────────────┐ │  │
	│  │  │  DATA ACQUISITION STATION                                        │ │  │
	│  │  │  • Windows 11 PC with MinKNOW software                          │ │  │
	│  │  │  • Real-time basecalling OFF (save raw FAST5)                   │ │  │
	│  │  │  • Transfer to NAS: 10 Gbps ethernet                            │ │  │
	│  │  └──────────────────────────────────────────────────────────────────┘ │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                              ▼ (10 Gbps Network)                              │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  NETWORK-ATTACHED STORAGE (NAS)                                        │  │
	│  │  • Capacity: 100 TB (RAID 6)                                           │  │
	│  │  • Raw FAST5 storage: 50 TB allocated                                  │  │
	│  │  • Encrypted at rest (AES-256)                                         │  │
	│  │  • Automated backup to tape (LTO-9)                                    │  │
	│  │  • 7-year retention (PMDA requirement)                                 │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                              ▼ (10 Gbps Network)                              │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  COMPUTATION AREA (Secure Server Room)                                 │  │
	│  │  ┌──────────────────────────────────────────────────────────────────┐  │  │
	│  │  │  RACK CONFIGURATION                                              │  │  │
	│  │  │                                                                  │  │  │
	│  │  │  ┌─────────────────────────────────────────────────────────┐    │  │  │
	│  │  │  │  DGX Spark #1 (Primary Production System)              │    │  │  │
	│  │  │  │  ────────────────────────────────────────────────────── │    │  │  │
	│  │  │  │  Hostname: dgx-spark-01.meiji.ac.jp                    │    │  │  │
	│  │  │  │  IP: 192.168.100.101 (private network)                 │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Hardware:                                              │    │  │  │
	│  │  │  │  • ARM Cortex (10× X925 + 10× A725)                    │    │  │  │
	│  │  │  │  • NVIDIA Blackwell GPU (6,144 CUDA)                   │    │  │  │
	│  │  │  │  • 128GB LPDDR5x unified memory                        │    │  │  │
	│  │  │  │  • 4TB NVMe M.2 SSD                                    │    │  │  │
	│  │  │  │  • Power: 240W (1U form factor)                        │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Software Stack:                                        │    │  │  │
	│  │  │  │  • DGX OS (Ubuntu 22.04 LTS for ARM)                   │    │  │  │
	│  │  │  │  • CUDA 12.4 (ARM build)                               │    │  │  │
	│  │  │  │  • Dorado 0.5+ (ARM-compiled or official build)        │    │  │  │
	│  │  │  │  • Kraken2 2.1.3 (compiled from source for ARM)        │    │  │  │
	│  │  │  │  • Python 3.11 + BioPython stack                       │    │  │  │
	│  │  │  │  • Docker + NVIDIA Container Runtime                   │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Role: Process odd-numbered sample runs                │    │  │  │
	│  │  │  │  Status: 24/7 operation, auto-restart on power loss    │    │  │  │
	│  │  │  └─────────────────────────────────────────────────────────┘    │  │  │
	│  │  │                                                                  │  │  │
	│  │  │  ┌─────────────────────────────────────────────────────────┐    │  │  │
	│  │  │  │  DGX Spark #2 (Development + Parallel System)          │    │  │  │
	│  │  │  │  ────────────────────────────────────────────────────── │    │  │  │
	│  │  │  │  Hostname: dgx-spark-02.meiji.ac.jp                    │    │  │  │
	│  │  │  │  IP: 192.168.100.102 (private network)                 │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Hardware: Same as DGX Spark #1                        │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Software Stack: Same as DGX Spark #1                  │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Role:                                                  │    │  │  │
	│  │  │  │  • Process even-numbered sample runs (parallel)        │    │  │  │
	│  │  │  │  • Software development and testing                    │    │  │  │
	│  │  │  │  • Backup for DGX Spark #1 during maintenance          │    │  │  │
	│  │  │  │  • Student training and workshops                      │    │  │  │
	│  │  │  │                                                         │    │  │  │
	│  │  │  │  Status: On-demand activation                          │    │  │  │
	│  │  │  └─────────────────────────────────────────────────────────┘    │  │  │
	│  │  │                                                                  │  │  │
	│  │  │  ┌─────────────────────────────────────────────────────────┐    │  │  │
	│  │  │  │  Network Switch (10 Gbps)                              │    │  │  │
	│  │  │  │  • Connects DGX Spark to NAS                           │    │  │  │
	│  │  │  │  • VLAN isolation for security                         │    │  │  │
	│  │  │  └─────────────────────────────────────────────────────────┘    │  │  │
	│  │  │                                                                  │  │  │
	│  │  │  ┌─────────────────────────────────────────────────────────┐    │  │  │
	│  │  │  │  UPS (Uninterruptible Power Supply)                    │    │  │  │
	│  │  │  │  • 2000VA capacity                                     │    │  │  │
	│  │  │  │  • 30-minute runtime at full load                      │    │  │  │
	│  │  │  │  • Graceful shutdown on extended outage                │    │  │  │
	│  │  │  └─────────────────────────────────────────────────────────┘    │  │  │
	│  │  │                                                                  │  │  │
	│  │  │  ┌─────────────────────────────────────────────────────────┐    │  │  │
	│  │  │  │  Environmental Monitoring                              │    │  │  │
	│  │  │  │  • Temperature: 20-22°C (68-72°F)                      │    │  │  │
	│  │  │  │  • Humidity: 40-60% RH                                 │    │  │  │
	│  │  │  │  • Precision AC unit (redundant)                       │    │  │  │
	│  │  │  └─────────────────────────────────────────────────────────┘    │  │  │
	│  │  └──────────────────────────────────────────────────────────────────┘  │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  NETWORK CONNECTIVITY                                                  │  │
	│  │  ┌──────────────────────┐                                             │  │
	│  │  │  Internet Firewall   │  ⚠️ DGX Spark ISOLATED (No internet access) │  │
	│  │  │  • Air-gapped mode   │     for maximum PMDA compliance             │  │
	│  │  │  • Updates via USB   │                                             │  │
	│  │  └──────────────────────┘                                             │  │
	│  │           ▲                                                            │  │
	│  │           │ (Optional: DGX Cloud access for benchmarking only)        │  │
	│  │           ▼                                                            │  │
	│  │  ┌──────────────────────┐                                             │  │
	│  │  │  Management PC       │  • Remote monitoring                        │  │
	│  │  │  • SSH access        │  • Log aggregation                          │  │
	│  │  │  • Monitoring tools  │  • Performance dashboards                   │  │
	│  │  └──────────────────────┘                                             │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  WORKFLOW ORCHESTRATION                                                       │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Sample Arrives → Load to DGX Spark #1 or #2 (round-robin)            │  │
	│  │       ▼                                                                 │  │
	│  │  Automated Pipeline Execution:                                         │  │
	│  │  1. Watchdog detects new FAST5 in /data/incoming/                     │  │
	│  │  2. Triggers pipeline script: /opt/minion/run_pipeline.sh             │  │
	│  │  3. Phases 1-6 execute sequentially on assigned DGX Spark             │  │
	│  │  4. Results written to /data/results/{run_id}/                        │  │
	│  │  5. Automatic sync to NAS for archival                                │  │
	│  │  6. Email/Slack notification to PI on completion                      │  │
	│  │       ▼                                                                 │  │
	│  │  Manual Review → Approve Report → Send to PMDA                        │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  MONITORING & ALERTING                                                        │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Real-Time Monitoring Dashboard                                        │  │
	│  │  • GPU utilization (nvidia-smi)                                        │  │
	│  │  • CPU temperature and frequency                                       │  │
	│  │  • Memory usage (unified 128GB pool)                                   │  │
	│  │  • Disk I/O (NVMe throughput)                                          │  │
	│  │  • Network throughput (NAS → DGX Spark)                                │  │
	│  │  • Pipeline stage progress                                             │  │
	│  │                                                                         │  │
	│  │  Alerts:                                                                │  │
	│  │  • PERV detection → Immediate Slack notification to PI                │  │
	│  │  • Pipeline failure → Email + SMS                                     │  │
	│  │  • Hardware issues (temperature >80°C) → Admin alert                  │  │
	│  │  • Disk space <10% → Warning                                          │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  POWER & ENVIRONMENTAL                                                        │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Power Consumption Analysis                                            │  │
	│  │  ────────────────────────────────────────────────────────────────────  │  │
	│  │  DGX Spark #1: 240W × 24h × 30 days = 172.8 kWh/month                │  │
	│  │  DGX Spark #2: 240W × 8h/day × 20 days = 38.4 kWh/month (part-time)  │  │
	│  │  NAS + Networking: ~500W × 24h × 30 days = 360 kWh/month             │  │
	│  │  Cooling (AC): ~1000W × 24h × 30 days = 720 kWh/month                │  │
	│  │  ────────────────────────────────────────────────────────────────────  │  │
	│  │  Total: ~1,291 kWh/month                                              │  │
	│  │  Cost (Japan rate: ¥30/kWh): ¥38,730/month (~$265 USD)               │  │
	│  │                                                                         │  │
	│  │  Carbon Footprint:                                                     │  │
	│  │  • 1,291 kWh × 0.5 kg CO₂/kWh = 645 kg CO₂/month                     │  │
	│  │  • Annual: 7.74 tons CO₂                                              │  │
	│  │                                                                         │  │
	│  │  Comparison to Cloud:                                                  │  │
	│  │  • AWS A100 instances: ~2,500 kWh/month for equivalent workload       │  │
	│  │  • DGX Spark 48% more energy-efficient                                │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘

---- 

## Diagram 3: Data Flow and Security Architecture

	┌─────────────────────────────────────────────────────────────────────────────┐
	│                     DATA SECURITY & PMDA COMPLIANCE                          │
	│                    (Meiji University On-Premises Model)                      │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  TIER 1: SAMPLE COLLECTION (Veterinary Facility)                             │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Porcine Sample                                                        │  │
	│  │  • Blood, tissue, or organ biopsy                                      │  │
	│  │  • Unique ID: MJ-XEN-2025-{NNN}                                        │  │
	│  │  • Metadata: Age, sex, health status, geographic origin               │  │
	│  │  • ⚠️ NO patient data (pigs only, not humans)                         │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                              ▼ (Physical Transport)                           │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Meiji University BSL-2 Laboratory                                     │  │
	│  │  • DNA/RNA extraction in biosafety cabinet                             │  │
	│  │  • Library preparation for MinION                                      │  │
	│  │  • ⚠️ Samples destroyed after sequencing (no long-term biobank)       │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	                                    │
	                                    ▼
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  TIER 2: SEQUENCING (Raw Data Generation)                                    │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  MinION Sequencer → FAST5/POD5 Files                                   │  │
	│  │  • File naming: {run_id}_{read_id}.fast5                               │  │
	│  │  • Metadata embedded: timestamp, flowcell ID, device serial           │  │
	│  │  • ⚠️ NO PHI/PII: Only sample ID, no names/locations                  │  │
	│  │  • Transfer: Direct cable → Data acquisition PC → NAS                 │  │
	│  │  • ✅ Never uploaded to cloud (PMDA compliance)                        │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	                                    │
	                                    ▼
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  TIER 3: STORAGE (Encrypted Data at Rest)                                    │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Network-Attached Storage (NAS)                                        │  │
	│  │  ─────────────────────────────────────────────────────────────────────  │  │
	│  │  Encryption:                                                            │  │
	│  │  • Algorithm: AES-256-XTS                                              │  │
	│  │  • Key management: Hardware Security Module (HSM)                     │  │
	│  │  • Keys rotated annually                                               │  │
	│  │  • Encrypted at rest, decrypted only during pipeline execution        │  │
	│  │                                                                         │  │
	│  │  Access Control:                                                        │  │
	│  │  • Authentication: LDAP integration with Meiji University directory   │  │
	│  │  • Authorization: Role-based (PI, students, technicians)              │  │
	│  │  • Audit logging: Every file access logged with timestamp + user      │  │
	│  │  • Network isolation: VLAN 100 (research only, no internet)           │  │
	│  │                                                                         │  │
	│  │  Retention Policy:                                                      │  │
	│  │  • Raw FAST5: 7 years (PMDA requirement for clinical trials)          │  │
	│  │  • Processed FASTQ: 5 years                                            │  │
	│  │  • Reports: 10 years (permanent archival)                              │  │
	│  │  • Backup: LTO-9 tape (offline, fireproof safe)                       │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	                                    │
	                                    ▼
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  TIER 4: PROCESSING (DGX Spark - Isolated Environment)                       │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  DGX Spark Security Architecture                                       │  │
	│  │  ─────────────────────────────────────────────────────────────────────  │  │
	│  │  Network Isolation:                                                     │  │
	│  │  • ✅ Connected to NAS (10 Gbps private network)                       │  │
	│  │  • ❌ NO internet access (air-gapped)                                  │  │
	│  │  • ❌ NO SSH from outside Meiji network                                │  │
	│  │  • ✅ USB updates only (software patches via USB drive)                │  │
	│  │                                                                         │  │
	│  │  OS Security:                                                           │  │
	│  │  • SELinux: Enforcing mode                                             │  │
	│  │  • Firewall: UFW (deny all except NAS)                                │  │
	│  │  • User accounts: PI + 2 authorized students only                     │  │
	│  │  • Sudo access: Restricted to PI only                                 │  │
	│  │                                                                         │  │
	│  │  Data Handling:                                                         │  │
	│  │  • Working directory: /data/{run_id}/ (temporary)                     │  │
	│  │  • Encrypted tmpfs for intermediate files                              │  │
	│  │  • Automatic deletion after pipeline completion + NAS sync            │  │
	│  │  • No data persists on DGX Spark after run (ephemeral)                │  │
	│  │                                                                         │  │
	│  │  Audit Trail:                                                           │  │
	│  │  • All pipeline executions logged                                      │  │
	│  │  • Syslog forwarded to central server                                 │  │
	│  │  • Immutable logs (append-only, cannot be modified)                   │  │
	│  │  • Retention: 3 years for PMDA audit                                  │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	                                    │
	                                    ▼
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  TIER 5: RESULTS & REPORTING (De-Identification for Publication)             │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Internal Use (Meiji University - Identified)                          │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  Report includes:                                                       │  │
	│  │  • Sample ID: MJ-XEN-2025-042                                          │  │
	│  │  • Metadata: Age (18 months), Sex (Female), Origin (Farm XYZ)         │  │
	│  │  • Pathogen results: PCV2 detected (247 reads)                         │  │
	│  │  • PERV status: NEGATIVE (critical)                                    │  │
	│  │  • AI risk: LOW (87.3% confidence)                                     │  │
	│  │  • Recommendation: APPROVED for transplantation                        │  │
	│  │                                                                         │  │
	│  │  Storage: NAS encrypted folder, access restricted to PI + PMDA        │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Public Dataset (Zenodo - De-Identified)                               │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  De-identification process:                                             │  │
	│  │  • Sample ID → Anonymous hash (SHA-256)                                │  │
	│  │  • Remove: Farm name, exact coordinates, dates                         │  │
	│  │  • Generalize: Age → Age bin (e.g., "12-24 months")                   │  │
	│  │  • Keep: Pathogen counts, AI predictions, QC metrics                  │  │
	│  │                                                                         │  │
	│  │  Public dataset format:                                                 │  │
	│  │  {                                                                      │  │
	│  │    "sample_id_hash": "a3f8b9c...",                                     │  │
	│  │    "age_bin": "12-24_months",                                          │  │
	│  │    "sex": "F",                                                          │  │
	│  │    "region": "Kanto",  # Generalized                                   │  │
	│  │    "pathogens": {"PCV2": 247, "TTV": 89},                              │  │
	│  │    "perv_status": "negative",                                           │  │
	│  │    "ai_risk": "low",                                                    │  │
	│  │    "qc_metrics": {...}                                                  │  │
	│  │  }                                                                      │  │
	│  │                                                                         │  │
	│  │  Published: Zenodo with DOI, CC BY 4.0 license                        │  │
	│  │  Use case: Train AI models, method validation, meta-analyses          │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  PMDA COMPLIANCE CHECKLIST                                                    │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  ✅ Data Residency: All clinical data on-premises (never cloud)        │  │
	│  │  ✅ Encryption: AES-256 at rest, TLS 1.3 in transit                    │  │
	│  │  ✅ Access Control: LDAP authentication, role-based authorization      │  │
	│  │  ✅ Audit Trail: Immutable logs, 3-year retention                      │  │
	│  │  ✅ Retention: 7-year raw data, 10-year reports                        │  │
	│  │  ✅ Backup: Automated daily, offsite tape backup                       │  │
	│  │  ✅ Incident Response: 24h notification protocol for PERV detection    │  │
	│  │  ✅ Validation: 91 PMDA pathogens, PPA >95%, NPA >98%                  │  │
	│  │  ✅ Documentation: Standard Operating Procedures (SOPs) documented     │  │
	│  │  ✅ Training: All personnel trained on data handling protocols         │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘

---- 

## Diagram 4: Benchmarking Study Design (DGX Spark vs A100)

	┌─────────────────────────────────────────────────────────────────────────────┐
	│              COMPREHENSIVE BENCHMARKING STUDY DESIGN                         │
	│                    (18-Month Research Timeline)                              │
	└─────────────────────────────────────────────────────────────────────────────┘
	
	PHASE 1: BASELINE ESTABLISHMENT (Months 1-3)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  Test Dataset: 30 Porcine Samples (Identical inputs for both platforms)      │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Sample Characteristics:                                               │  │
	│  │  • FAST5 size: 30-50 GB per sample                                     │  │
	│  │  • Read count: 3-5 million reads                                       │  │
	│  │  • Quality range: Q9-Q15 (MinION typical)                              │  │
	│  │  • Pathogen diversity: 0-8 pathogens per sample                        │  │
	│  │  • PERV status: 10 positive, 20 negative (known truth)                │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│         ┌──────────────────────────────┐   ┌───────────────────────────────┐ │
	│         │   DGX Cloud A100 (x86)       │   │   Current AWS (Baseline)      │ │
	│         │   ─────────────────────────  │   │   ──────────────────────────  │ │
	│         │   Platform: NVIDIA DGX Cloud │   │   Platform: AWS EC2           │ │
	│         │   GPU: A100 80GB             │   │   GPU: g4dn.xlarge (T4 16GB)  │ │
	│         │   CPU: x86 Xeon              │   │   CPU: Intel Xeon (x86)       │ │
	│         │   RAM: 512 GB                │   │   RAM: 128 GB (r5.4xlarge)    │ │
	│         │                              │   │                               │ │
	│         │   Process all 30 samples     │   │   Process all 30 samples      │ │
	│         │   Collect metrics:           │   │   Collect metrics:            │ │
	│         │   • Runtime per phase        │   │   • Runtime per phase         │ │
	│         │   • GPU utilization          │   │   • GPU utilization           │ │
	│         │   • Memory bandwidth         │   │   • Memory bandwidth          │ │
	│         │   • Accuracy (PPA/NPA)       │   │   • Accuracy (PPA/NPA)        │ │
	│         │   • Cost per sample          │   │   • Cost per sample           │ │
	│         └──────────────────────────────┘   └───────────────────────────────┘ │
	│                                                                               │
	│  Deliverable: Benchmark Report #1 (DGX Cloud A100 vs AWS T4)                │
	│  • A100 expected to be 4-6× faster than T4 for basecalling                  │
	│  • Establish x86 performance baseline for ARM comparison                     │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	PHASE 2: DGX SPARK ARM PORTING & OPTIMIZATION (Months 4-9)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  Month 4-6: Software Compatibility Testing                                    │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  DGX Spark #1 (ARM Software Stack Development)                         │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  Week 1-2: Dorado Basecalling on ARM                                   │  │
	│  │  • Check Oxford Nanopore for official aarch64 build                    │  │
	│  │  • If unavailable: Compile from source (if possible)                   │  │
	│  │  • Test with 5 samples, measure GPU utilization                        │  │
	│  │  • Compare basecalling accuracy vs x86 (Q-score distribution)          │  │
	│  │  • Result: GO/NO-GO decision for ARM basecalling                       │  │
	│  │                                                                         │  │
	│  │  Week 3-4: Kraken2 Compilation for ARM                                 │  │
	│  │  • Clone Kraken2 GitHub repo                                           │  │
	│  │  • Compile for aarch64 architecture                                    │  │
	│  │  • Load PMDA database (20-30 GB) into 128GB unified memory            │  │
	│  │  • Test classification speed (reads/minute)                            │  │
	│  │  • Validate taxonomy assignments vs x86 baseline                       │  │
	│  │  • Result: Performance ratio (ARM / x86)                               │  │
	│  │                                                                         │  │
	│  │  Week 5-6: BLAST & Supporting Tools                                    │  │
	│  │  • Install BLAST+ for ARM (conda or source compilation)                │  │
	│  │  • Test PERV typing script (Python, should work on ARM)                │  │
	│  │  • Validate Minimap2 (host removal) on ARM                             │  │
	│  │  • End-to-end test: 10 samples through full pipeline                  │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  Month 7-9: Performance Benchmarking (50 Samples)                            │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Test Matrix (Same 50 samples on both platforms)                       │  │
	│  │                                                                         │  │
	│  │      DGX Spark #1+#2 (ARM)           vs       DGX Cloud A100 (x86)     │  │
	│  │      ──────────────────────                   ─────────────────────     │  │
	│  │      25 samples each (parallel)               50 samples (sequential)  │  │
	│  │                                                                         │  │
	│  │  Metrics Collected:                                                     │  │
	│  │  ┌──────────────────┬──────────────────┬──────────────────────────┐    │  │
	│  │  │ Metric           │ DGX Spark (ARM)  │ A100 Cloud (x86)         │    │  │
	│  │  ├──────────────────┼──────────────────┼──────────────────────────┤    │  │
	│  │  │ Phase 1 Runtime  │ ? hours          │ 1-2 hours (baseline)     │    │  │
	│  │  │ Phase 4 Runtime  │ ? hours          │ 2.5-3 hours (baseline)   │    │  │
	│  │  │ Total Runtime    │ ? hours          │ 6-7 hours (baseline)     │    │  │
	│  │  │ GPU Utilization  │ ?%               │ 90-95%                   │    │  │
	│  │  │ Memory Bandwidth │ 273 GB/s (spec)  │ 2,039 GB/s (spec)        │    │  │
	│  │  │ Power/Sample     │ 240W × ?h        │ 500W × 6h                │    │  │
	│  │  │ Accuracy (PPA)   │ ?%               │ 95.2% (baseline)         │    │  │
	│  │  │ Accuracy (NPA)   │ ?%               │ 98.1% (baseline)         │    │  │
	│  │  │ PERV Detection   │ 100% or fail     │ 100% (critical)          │    │  │
	│  │  │ Cost/Sample      │ $0 (amortized)   │ $8-12 (cloud credits)    │    │  │
	│  │  └──────────────────┴──────────────────┴──────────────────────────┘    │  │
	│  │                                                                         │  │
	│  │  Statistical Analysis:                                                  │  │
	│  │  • Paired t-test for runtime comparison                                │  │
	│  │  • Cohen's d for effect size (ARM vs x86)                              │  │
	│  │  • 95% confidence intervals                                            │  │
	│  │  • Power analysis (sample size justification)                          │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  Deliverable: Benchmark Report #2 (DGX Spark ARM vs A100 x86)               │
	│  • First comprehensive ARM genomics benchmark published                      │
	│  • Recommendation: When to use ARM vs x86 for bioinformatics                │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	PHASE 3: PRODUCTION VALIDATION (Months 10-15)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  Clinical Sample Processing: 350 Samples Total                               │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Deployment Strategy:                                                   │  │
	│  │                                                                         │  │
	│  │  If DGX Spark Performance ≥ 80% of A100:                               │  │
	│  │  ✅ Use DGX Spark for ALL 350 samples (on-premises production)         │  │
	│  │  ✅ Use A100 cloud only for burst capacity during outbreaks            │  │
	│  │  ✅ Publish: "DGX Spark is viable for production genomics"             │  │
	│  │                                                                         │  │
	│  │  If DGX Spark Performance < 80% of A100:                               │  │
	│  │  ⚠️ Hybrid: DGX Spark for Phases 2-3, 5-6 (CPU tasks)                 │  │
	│  │  ⚠️ A100 cloud for Phase 1 (basecalling) + Phase 4 (if needed)        │  │
	│  │  ⚠️ Publish: "ARM shows promise but not yet ready for GPU genomics"   │  │
	│  │                                                                         │  │
	│  │  If Critical Tools Fail on ARM (e.g., Dorado unavailable):             │  │
	│  │  ❌ Use A100 cloud for ALL 350 samples                                 │  │
	│  │  ❌ DGX Spark relegated to non-GPU tasks or student training           │  │
	│  │  ❌ Publish: "Negative results - ARM compatibility barriers"           │  │
	│  │  ❌ Still valuable: Informs community about ARM limitations            │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  Data Collection (All 350 Samples):                                          │
	│  • Pathogen abundance matrix (91 × 350)                                      │
	│  • AI training dataset (features + outbreak labels)                          │
	│  • Performance logs (runtime, resource usage)                                │
	│  • Cost accounting (DGX Spark power vs A100 cloud credits)                   │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	PHASE 4: AI MODEL TRAINING & ANALYSIS (Months 10-15)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  Transformer Model for Pathogen Risk Prediction                              │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Training Platform: DGX Cloud (2-4× A100 GPUs)                         │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  Dataset: 10,000+ samples (350 from this study + 9,650 from literature)│  │
	│  │  Model: BERT-based transformer (768-dimensional embeddings)            │  │
	│  │  Training time: 6.5 days on 2× A100 GPUs (316 GPU hours)               │  │
	│  │  Hyperparameter tuning: 10 experiments (300 GPU hours)                 │  │
	│  │  Total GPU hours: ~500 hours                                           │  │
	│  │                                                                         │  │
	│  │  Why A100 (not DGX Spark) for training:                                │  │
	│  │  • Distributed training requires 2-4× GPUs (DGX Spark = 1 GPU each)    │  │
	│  │  • Faster iteration (A100 3-5× faster than Blackwell integrated GPU)   │  │
	│  │  • Proven PyTorch/TensorFlow support on x86                            │  │
	│  │                                                                         │  │
	│  │  Inference Validation (Both Platforms):                                 │  │
	│  │  • Deploy trained model to DGX Spark (ARM) for inference              │  │
	│  │  • Compare inference latency: ARM vs x86                               │  │
	│  │  • Validate accuracy is identical (model architecture-agnostic)        │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	PHASE 5: DISSEMINATION (Months 16-18)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  Publications & Presentations                                                 │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  Paper 1: Technical Implementation (Target: Bioinformatics)            │  │
	│  │  "GPU-Accelerated Pathogen Detection for Xenotransplantation Safety"  │  │
	│  │  • Focus: Pipeline architecture, A100 optimization, clinical results  │  │
	│  │  • Submit: Month 12, Accept: Month 15, Publish: Month 18              │  │
	│  │                                                                         │  │
	│  │  Paper 2: ARM Benchmark (Target: Nature Comp Sci or SC '26)           │  │
	│  │  "Evaluating ARM-based DGX Spark for Production Genomics Workloads"   │  │
	│  │  • Focus: DGX Spark vs A100 comparison, performance analysis          │  │
	│  │  • Submit: Month 15, Accept: Month 18, Publish: Month 20              │  │
	│  │  • ⭐ Flagship deliverable - first ARM genomics benchmark             │  │
	│  │                                                                         │  │
	│  │  Paper 3: AI Model (Target: PLOS Comp Bio)                            │  │
	│  │  "Transformer-Based Pathogen Outbreak Prediction"                     │  │
	│  │  • Focus: AI methodology, 80%+ accuracy, clinical validation          │  │
	│  │  • Submit: Month 15, Accept: Month 18, Publish: Month 20              │  │
	│  │                                                                         │  │
	│  │  Conference Presentations:                                             │  │
	│  │  • NVIDIA GTC 2026 (March): Poster on DGX Spark benchmark             │  │
	│  │  • ISMB 2026 (July): Workshop on "GPU Genomics"                       │  │
	│  │  • SC '26 (November): Talk on ARM scientific computing (if accepted)  │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘
	
	EXPECTED OUTCOMES (Best, Middle, Worst Case)
	┌──────────────────────────────────────────────────────────────────────────────┐
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  BEST CASE (30% probability)                                           │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  • DGX Spark performance ≥ 90% of A100                                 │  │
	│  │  • All tools (Dorado, Kraken2, BLAST) work flawlessly on ARM          │  │
	│  │  • Pipeline runtime: 7-8 hours (vs 6-7h on A100)                       │  │
	│  │  • Power efficiency: 50% better than A100 cloud                        │  │
	│  │  • Publication: Nature Computational Science (high impact)             │  │
	│  │  • Conclusion: "ARM is ready for production genomics"                  │  │
	│  │  • Impact: Accelerates ARM adoption in bioinformatics community        │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  MIDDLE CASE (50% probability) ⭐ MOST LIKELY                          │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  • DGX Spark performance 70-85% of A100                                │  │
	│  │  • Dorado: Works but 20-30% slower than x86                            │  │
	│  │  • Kraken2: Compiles successfully, similar performance                 │  │
	│  │  • Pipeline runtime: 9-11 hours (vs 6-7h on A100)                      │  │
	│  │  • Hybrid deployment: DGX Spark for most tasks, A100 for GPU bursts   │  │
	│  │  • Publication: IEEE TCBB or SC '26 (good journals)                    │  │
	│  │  • Conclusion: "ARM shows promise, optimization needed"                │  │
	│  │  • Impact: Provides roadmap for ARM genomics maturation                │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	│                                                                               │
	│  ┌────────────────────────────────────────────────────────────────────────┐  │
	│  │  WORST CASE (20% probability)                                          │  │
	│  │  ────────────────────────────────────────────────────────────────────   │  │
	│  │  • Dorado has no ARM build, cannot compile from source                │  │
	│  │  • DGX Spark relegated to CPU-only tasks (Phases 2-3, 5-6)            │  │
	│  │  • Main pipeline runs on A100 cloud (all 350 samples)                 │  │
	│  │  • DGX Spark used for: Student training, non-GPU preprocessing        │  │
	│  │  • Publication: Negative results paper (still valuable)                │  │
	│  │  • Conclusion: "ARM not yet ready, community should wait"             │  │
	│  │  • Impact: Saves other researchers from wasting time on ARM porting   │  │
	│  │  • NVIDIA value: Identifies DGX Spark software ecosystem gaps         │  │
	│  └────────────────────────────────────────────────────────────────────────┘  │
	└──────────────────────────────────────────────────────────────────────────────┘

---- 

## Notes for Grant Application

**Include these diagrams in:**
- **Appendix B:** "System Architecture Diagrams"
- **Section 3 (Technical Approach):** Reference specific diagrams when describing methods

**Conversion to Visual Diagrams:**
- Use **draw.io** (https://app.diagrams.net/) to convert ASCII art to professional diagrams
- Use **Lucidchart** for flowcharts
- Use **Microsoft PowerPoint** with SmartArt for quick conversions
- Export as high-resolution PNG or vector PDF for submission

**Key Messages These Diagrams Convey:**
1. ✅ **Professional deployment plan** (not just "we'll try it and see")
2. ✅ **PMDA compliance** is central to architecture
3. ✅ **Rigorous benchmarking** methodology
4. ✅ **Contingency planning** for ARM compatibility risks
5. ✅ **Real-world production** deployment (not just academic exercise)

---- 

**Total Diagrams:** 4 comprehensive architecture visualizations
**Pages:** \~8-10 pages of detailed technical diagrams
**Format:** ASCII art (easily convertible to visual diagrams)
