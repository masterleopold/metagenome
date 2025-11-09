#!/bin/bash
# MinION Pipeline Deployment Script

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ENVIRONMENT="${ENVIRONMENT:-production}"
AWS_REGION="${AWS_REGION:-ap-northeast-1}"
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

usage() {
    echo "Usage: $0 [COMMAND] [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  deploy       Full deployment of the pipeline"
    echo "  update       Update existing deployment"
    echo "  rollback     Rollback to previous version"
    echo "  destroy      Remove all resources"
    echo "  validate     Validate deployment configuration"
    echo "  status       Check deployment status"
    echo ""
    echo "Options:"
    echo "  -e, --environment   Environment (development/staging/production)"
    echo "  -r, --region        AWS region"
    echo "  -p, --profile       AWS profile"
    echo "  -y, --yes           Skip confirmation prompts"
    echo "  -h, --help          Show this help message"
    echo ""
    echo "Environment Variables:"
    echo "  ENVIRONMENT         Default: production"
    echo "  AWS_REGION          Default: ap-northeast-1"
    echo "  AWS_PROFILE         AWS profile to use"
    exit 0
}

# Parse arguments
COMMAND=""
SKIP_CONFIRM=false

while [[ $# -gt 0 ]]; do
    case $1 in
        deploy|update|rollback|destroy|validate|status)
            COMMAND=$1
            shift
            ;;
        -e|--environment)
            ENVIRONMENT="$2"
            shift 2
            ;;
        -r|--region)
            AWS_REGION="$2"
            shift 2
            ;;
        -p|--profile)
            export AWS_PROFILE="$2"
            shift 2
            ;;
        -y|--yes)
            SKIP_CONFIRM=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            error "Unknown option: $1"
            ;;
    esac
done

[[ -z "$COMMAND" ]] && usage

# Confirmation prompt
confirm() {
    if [[ "$SKIP_CONFIRM" == "false" ]]; then
        read -p "Are you sure? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Aborted"
            exit 1
        fi
    fi
}

# Validate prerequisites
validate_prerequisites() {
    log "Validating prerequisites..."

    # Check required tools
    local required_tools=("aws" "terraform" "python3" "npm" "docker")
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error "$tool is required but not installed"
        fi
    done

    # Check AWS credentials
    if ! aws sts get-caller-identity &> /dev/null; then
        error "AWS credentials not configured"
    fi

    # Check Python packages
    python3 -c "import boto3, click, jinja2" 2>/dev/null || {
        warn "Required Python packages not installed. Installing..."
        pip3 install boto3 click jinja2 tabulate
    }

    log "Prerequisites validated"
}

# Deploy infrastructure
deploy_infrastructure() {
    log "Deploying infrastructure with Terraform..."

    cd "$PROJECT_ROOT/infrastructure/terraform"

    # Initialize Terraform
    terraform init -backend-config="key=minion-pipeline-${ENVIRONMENT}"

    # Plan deployment
    terraform plan -var="environment=${ENVIRONMENT}" -var="region=${AWS_REGION}" -out=tfplan

    # Apply
    info "Terraform plan:"
    terraform show tfplan

    if [[ "$SKIP_CONFIRM" == "false" ]]; then
        confirm
    fi

    terraform apply tfplan

    # Save outputs
    terraform output -json > "$PROJECT_ROOT/infrastructure/outputs.json"

    log "Infrastructure deployed"
}

# Build and push Docker images (if using containers)
build_images() {
    log "Building AMIs..."

    cd "$PROJECT_ROOT/ec2_setup"

    # Build basecalling AMI
    ./build_basecalling_ami.sh

    # Build analysis AMI
    ./build_analysis_ami.sh

    log "AMIs built successfully"
}

# Deploy Lambda functions
deploy_lambda_functions() {
    log "Deploying Lambda functions..."

    cd "$PROJECT_ROOT/lambda"

    # Package functions
    for dir in orchestration phases monitoring api; do
        cd "$dir"
        zip -r "../${dir}.zip" *.py
        cd ..
    done

    # Deploy each function
    for func in orchestration phases monitoring api; do
        aws lambda update-function-code \
            --function-name "minion-${func}-${ENVIRONMENT}" \
            --zip-file "fileb://${func}.zip" \
            --region "$AWS_REGION"
    done

    log "Lambda functions deployed"
}

# Deploy Step Functions
deploy_step_functions() {
    log "Deploying Step Functions state machine..."

    # Replace variables in template
    sed -e "s/\${AWS_REGION}/${AWS_REGION}/g" \
        -e "s/\${AWS_ACCOUNT}/${AWS_ACCOUNT_ID}/g" \
        "$PROJECT_ROOT/templates/stepfunctions/workflow_definition.json" > /tmp/workflow.json

    # Create or update state machine
    aws stepfunctions update-state-machine \
        --state-machine-arn "arn:aws:states:${AWS_REGION}:${AWS_ACCOUNT_ID}:stateMachine:minion-pipeline-${ENVIRONMENT}" \
        --definition file:///tmp/workflow.json \
        --region "$AWS_REGION" || \
    aws stepfunctions create-state-machine \
        --name "minion-pipeline-${ENVIRONMENT}" \
        --definition file:///tmp/workflow.json \
        --role-arn "arn:aws:iam::${AWS_ACCOUNT_ID}:role/StepFunctionsRole-${ENVIRONMENT}" \
        --region "$AWS_REGION"

    log "Step Functions deployed"
}

# Setup databases
setup_databases() {
    log "Setting up reference databases..."

    # Get EFS mount point from Terraform output
    EFS_ID=$(cd "$PROJECT_ROOT/infrastructure/terraform" && terraform output -raw efs_id)

    # Mount EFS on bastion or setup instance
    # This would typically be done on an EC2 instance
    info "EFS ID: $EFS_ID"
    info "Run database setup on an EC2 instance with EFS mounted"

    log "Database setup instructions provided"
}

# Deploy monitoring
deploy_monitoring() {
    log "Deploying monitoring stack..."

    # Deploy CloudFormation stack for monitoring
    aws cloudformation deploy \
        --template-file "$PROJECT_ROOT/templates/cloudformation/monitoring_stack.yaml" \
        --stack-name "minion-monitoring-${ENVIRONMENT}" \
        --parameter-overrides \
            Environment="${ENVIRONMENT}" \
            AlertEmail="${ALERT_EMAIL:-admin@example.com}" \
        --capabilities CAPABILITY_IAM \
        --region "$AWS_REGION"

    log "Monitoring deployed"
}

# Validate deployment
validate_deployment() {
    log "Validating deployment..."

    # Check infrastructure
    info "Checking infrastructure..."
    cd "$PROJECT_ROOT/infrastructure/terraform"
    terraform validate

    # Check Lambda functions
    info "Checking Lambda functions..."
    for func in orchestration phases monitoring api; do
        aws lambda get-function \
            --function-name "minion-${func}-${ENVIRONMENT}" \
            --region "$AWS_REGION" &>/dev/null || \
            warn "Lambda function minion-${func}-${ENVIRONMENT} not found"
    done

    # Check Step Functions
    info "Checking Step Functions..."
    aws stepfunctions describe-state-machine \
        --state-machine-arn "arn:aws:states:${AWS_REGION}:${AWS_ACCOUNT_ID}:stateMachine:minion-pipeline-${ENVIRONMENT}" \
        --region "$AWS_REGION" &>/dev/null || \
        warn "State machine not found"

    # Check databases
    info "Checking EFS..."
    EFS_ID=$(cd "$PROJECT_ROOT/infrastructure/terraform" && terraform output -raw efs_id 2>/dev/null || echo "")
    if [[ -n "$EFS_ID" ]]; then
        aws efs describe-file-systems --file-system-id "$EFS_ID" --region "$AWS_REGION" &>/dev/null || \
            warn "EFS not accessible"
    fi

    log "Validation complete"
}

# Get deployment status
get_status() {
    log "Deployment Status for ${ENVIRONMENT} environment"
    echo ""

    # Infrastructure status
    echo "Infrastructure:"
    cd "$PROJECT_ROOT/infrastructure/terraform"
    terraform output &>/dev/null && echo "  ✓ Terraform: Deployed" || echo "  ✗ Terraform: Not deployed"

    # Lambda status
    echo ""
    echo "Lambda Functions:"
    for func in orchestration phases monitoring api; do
        aws lambda get-function \
            --function-name "minion-${func}-${ENVIRONMENT}" \
            --region "$AWS_REGION" &>/dev/null && \
            echo "  ✓ ${func}: Active" || \
            echo "  ✗ ${func}: Not found"
    done

    # Step Functions status
    echo ""
    echo "Step Functions:"
    aws stepfunctions describe-state-machine \
        --state-machine-arn "arn:aws:states:${AWS_REGION}:${AWS_ACCOUNT_ID}:stateMachine:minion-pipeline-${ENVIRONMENT}" \
        --region "$AWS_REGION" &>/dev/null && \
        echo "  ✓ State Machine: Active" || \
        echo "  ✗ State Machine: Not found"

    # Recent executions
    echo ""
    echo "Recent Executions:"
    aws stepfunctions list-executions \
        --state-machine-arn "arn:aws:states:${AWS_REGION}:${AWS_ACCOUNT_ID}:stateMachine:minion-pipeline-${ENVIRONMENT}" \
        --max-results 5 \
        --region "$AWS_REGION" 2>/dev/null | \
        jq -r '.executions[] | "  \(.name): \(.status)"' || \
        echo "  No recent executions"
}

# Rollback deployment
rollback_deployment() {
    log "Rolling back deployment..."

    warn "This will revert to the previous Terraform state"
    confirm

    cd "$PROJECT_ROOT/infrastructure/terraform"

    # Get previous state
    terraform state pull > current_state.json

    # Rollback (simplified - in production, use proper versioning)
    warn "Manual rollback required. Use terraform state management commands."

    log "Rollback instructions provided"
}

# Destroy deployment
destroy_deployment() {
    log "Destroying deployment..."

    error "This will delete all resources for ${ENVIRONMENT} environment!"
    confirm

    # Destroy infrastructure
    cd "$PROJECT_ROOT/infrastructure/terraform"
    terraform destroy -var="environment=${ENVIRONMENT}" -var="region=${AWS_REGION}" -auto-approve

    # Delete CloudFormation stacks
    aws cloudformation delete-stack \
        --stack-name "minion-monitoring-${ENVIRONMENT}" \
        --region "$AWS_REGION"

    log "Deployment destroyed"
}

# Main execution
main() {
    case $COMMAND in
        deploy)
            validate_prerequisites
            deploy_infrastructure
            build_images
            deploy_lambda_functions
            deploy_step_functions
            setup_databases
            deploy_monitoring
            validate_deployment
            get_status
            log "Deployment completed successfully!"
            ;;
        update)
            validate_prerequisites
            deploy_lambda_functions
            deploy_step_functions
            validate_deployment
            log "Update completed successfully!"
            ;;
        rollback)
            rollback_deployment
            ;;
        destroy)
            destroy_deployment
            ;;
        validate)
            validate_prerequisites
            validate_deployment
            ;;
        status)
            get_status
            ;;
        *)
            error "Unknown command: $COMMAND"
            ;;
    esac
}

# Run main function
main