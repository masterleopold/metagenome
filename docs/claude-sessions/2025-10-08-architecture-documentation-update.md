# Architecture Documentation Update - Lambda + EC2 Custom AMI

**Date**: 2025-10-08
**Session Type**: Documentation Update
**Status**: Completed

## Overview

Updated all core documentation files to accurately reflect the Lambda + EC2 custom AMI architecture, removing all references to Docker containers and AWS Batch.

## Changes Made

### 1. MinION_AWS_クラウド解析パイプライン詳細設計.md
**Version**: Updated to reflect containerless architecture

**Key Updates**:
- Updated architecture section from "AWS Batch並列解析" to "Lambda（16個の関数、オーケストレーション）+ EC2カスタムAMI"
- Replaced section 2.2 "Dockerコンテナ戦略" with "カスタムAMI戦略"
- Added AMI management details with Terraform examples
- Workflow section changed from Nextflow to Step Functions + Lambda
- Updated all instance types and specifications
- Added EFS mount configuration details

**Technical Details Added**:
```yaml
AMI_Strategy:
  Custom_AMIs:
    Basecalling_AMI:
      Name: "minion-basecalling-ami-v1.0"
      Base: "Deep Learning AMI GPU PyTorch（Amazon Linux 2023）"
      Pre-installed:
        - "NVIDIA CUDA 12.x + cuDNN"
        - "Dorado 0.5.0+"
        - "Analysis scripts at /opt/minion/scripts/"
```

### 2. MinION_メタゲノム解析パイプライン完全仕様書.md
**Version**: 2.0 → 2.1

**Key Updates**:
- Added header note explaining containerless architecture
- Updated technology stack table (Lambda count: 10 → 16, added Step Functions, SSM)
- Added new section 1.0 "実装アプローチ" explaining containerless architecture
- Replaced ALL "AWS Batch Job定義" sections with "Lambda関数によるEC2起動"
- Updated basecalling section with actual Python Lambda code example
- Added reference notes to remaining JSON definitions
- Updated terminology throughout

**Code Example Added**:
```python
# lambda/phases/trigger_basecalling.py
def lambda_handler(event, context):
    """Launch GPU EC2 instance for basecalling with Dorado"""
    run_id = event['run_id']
    user_data = f"""#!/bin/bash
    export RUN_ID='{run_id}'
    mount -t efs {os.environ['EFS_DNS']}:/dorado_models /mnt/dorado_models
    /opt/minion/scripts/phase1_basecalling/basecall_duplex.sh
    echo "sudo shutdown -h now" | at now + 12 hours
    """
    response = ec2.run_instances(
        ImageId=os.environ['BASECALLING_AMI_ID'],
        InstanceType='g4dn.xlarge',
        UserData=user_data,
        ...
    )
```

### 3. MinION_メタゲノム解析パイプライン開発計画書.md
**Version**: 1.0 → 2.0

**Key Updates**:
- Added implementation note at header: "Lambda関数によるEC2インスタンス起動とカスタムAMIによる解析"
- Updated technology stack table:
  - Lambda functions: 10 → 16
  - Added Step Functions and SSM
  - Updated EC2 instance types (g4dn.xlarge, r5/c5)
- Updated all phase descriptions with specific EC2 instance types
- Expanded Lambda function details with actual file paths
- Complete rewrite of UserData template section
- Updated cost model:
  - Monthly cost: ¥15,600 → ¥11,400 (with Spot Instances)
  - Per sample: ¥2,400 → ¥830 (with Spot Instances)
- Added Spot Instance interruption to risk management
- Updated next steps with specific AMI creation details

### 4. MinION_解析パイプライン_クイックリファレンスガイド.md
**Version**: 1.0 → 2.0

**Key Updates**:
- Added implementation note: "Lambda + EC2カスタムAMI（Dockerコンテナ不使用）"
- Updated system architecture diagram to show Lambda orchestration
- Added specific EC2 instance types for each phase:
  - Phase 1: g4dn.xlarge (GPU)
  - Phase 2: t3.large
  - Phase 3: r5.4xlarge
  - Phase 4: 4x parallel EC2
  - Phase 5-6: t3.large
- Updated phase descriptions to clarify Lambda triggers EC2
- Added architecture features section:
  - Auto-termination after completion
  - 70% cost savings with Spot Instances
  - EFS for shared reference databases
  - Serverless + on-demand EC2
- Updated troubleshooting section with Lambda EC2 launch code example

### 5. CLAUDE.md
**Updates**:
- Added architecture implementation details to Pipeline Architecture section
- Noted Lambda + EC2 custom AMI approach (no Docker)
- Added specific instance types for each phase
- Added key architecture decisions section
- Updated "Recently Updated" section with today's changes

### 6. README.md
**Updates**:
- Updated Software requirements: removed "Docker (optional)", added note about custom AMIs
- Expanded System Architecture section with detailed containerless explanation
- Added ASCII diagram showing Lambda orchestration flow
- Added key features list explaining the architecture approach

## Protocol Documents Status

### No Updates Required:
The following laboratory protocol documents (Protocol 00-10 series) were checked and **do not require updates**:
- MinION_Protocol_00_目次とマスタードキュメント.md
- MinION_Protocol_01-10 (all 10 protocol documents)

**Reason**: These documents cover wet-lab procedures (sample collection, DNA/RNA extraction, library prep, sequencing) and do not contain AWS architecture or computational pipeline details.

## Technical Architecture Summary

### Implementation Approach
- **Containerless**: No Docker, no AWS Batch, no container registries
- **Serverless Orchestration**: 16 Lambda functions control workflow
- **On-Demand Compute**: EC2 instances launched per phase
- **Custom AMIs**: Pre-installed tools at `/opt/minion/scripts/`
- **Shared Storage**: EFS for reference databases (Kraken2, BLAST, PERV)
- **Cost Optimization**: Spot Instances for 70% savings

### EC2 Instance Types by Phase
```yaml
Phase 1 (Basecalling):    g4dn.xlarge   # GPU, NVIDIA T4, 16GB VRAM
Phase 2 (QC):             t3.large      # General purpose
Phase 3 (Host Removal):   r5.4xlarge    # Memory optimized, 128GB RAM
Phase 4 (Pathogen):       4x parallel   # Multiple instance types
  - Kraken2:              r5.xlarge     # 32GB RAM
  - BLAST:                c5.4xlarge    # 16 vCPU
  - De novo:              c5.8xlarge    # 32 vCPU
  - PERV:                 t3.xlarge     # 4 vCPU
Phase 5 (Quantification): t3.large      # General purpose
Phase 6 (Reporting):      t3.large      # General purpose
```

### Lambda Functions Structure
```
lambda/
├── orchestrator/      (3 files) - S3 trigger, workflow start
├── phases/            (7 files) - EC2 launch per phase
├── validators/        (3 files) - Input validation, QC checks
└── metadata/          (3 files) - RDS updates, SNS notifications
```

### Cost Impact
- **Monthly Cost**: ¥11,400 (down from ¥15,600)
- **Per Sample Cost**: ¥830 with Spot Instances (down from ¥2,400)
- **Annual Savings**: Spot Instances provide 70% reduction on EC2 costs
- **Lambda Costs**: $0 (within free tier)

## Files Modified
```
modified:   CLAUDE.md
modified:   README.md
modified:   md/MinION_AWS_クラウド解析パイプライン詳細設計.md
modified:   md/MinION_メタゲノム解析パイプライン完全仕様書.md
modified:   md/MinION_メタゲノム解析パイプライン開発計画書.md
modified:   md/MinION_解析パイプライン_クイックリファレンスガイド.md
```

## Verification
- ✅ All Docker/container references removed
- ✅ AWS Batch references replaced with Lambda + EC2
- ✅ Nextflow references replaced with Step Functions
- ✅ Instance types specified for each phase
- ✅ Cost models updated with Spot Instance savings
- ✅ Lambda function structure documented
- ✅ EFS mount strategy documented
- ✅ Auto-termination behavior documented

## Next Steps
1. Update codebase implementation to match documentation
2. Create custom AMI build scripts
3. Implement Lambda functions for each phase
4. Set up EFS with reference databases
5. Configure Step Functions state machine
6. Test end-to-end with sample data

## Notes
- All documents now consistently describe the containerless architecture
- Protocol documents (00-10) remain unchanged as they cover lab procedures only
- Documentation versions updated appropriately (1.0 → 2.0 or 2.1)
- Revision history added to all modified documents
