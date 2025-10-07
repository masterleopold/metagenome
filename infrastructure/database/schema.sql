-- MinION Metagenomics Pipeline Database Schema
-- PostgreSQL 15.x
-- PMDA 91 Pathogen Screening System
-- Created: 2025-01-08

-- Enable required extensions
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- ===== PMDA 91 Pathogens Master Table =====
CREATE TABLE IF NOT EXISTS pmda_pathogens (
    id SERIAL PRIMARY KEY,
    pathogen_code VARCHAR(20) UNIQUE NOT NULL,
    pathogen_name_ja VARCHAR(255) NOT NULL,
    pathogen_name_en VARCHAR(255) NOT NULL,
    pathogen_type VARCHAR(50) NOT NULL CHECK (pathogen_type IN ('virus', 'bacteria', 'fungus', 'parasite', 'prion')),
    risk_level VARCHAR(20) CHECK (risk_level IN ('BSL-1', 'BSL-2', 'BSL-3', 'BSL-4')),
    reference_genome_id VARCHAR(100),
    reference_genome_length INTEGER,
    detection_priority INTEGER CHECK (detection_priority BETWEEN 1 AND 5),
    is_zoonotic BOOLEAN DEFAULT true,
    clinical_significance TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT unique_pathogen_code UNIQUE (pathogen_code)
);

CREATE INDEX idx_pathogen_type ON pmda_pathogens(pathogen_type);
CREATE INDEX idx_risk_level ON pmda_pathogens(risk_level);

-- ===== Workflow Executions Table =====
CREATE TABLE IF NOT EXISTS workflow_executions (
    run_id VARCHAR(255) PRIMARY KEY,
    sample_id VARCHAR(255),
    donor_id VARCHAR(255),
    status VARCHAR(50) NOT NULL DEFAULT 'initiated' CHECK (status IN (
        'initiated', 'uploading', 'basecalling', 'qc_analysis', 'host_removal',
        'pathogen_detection', 'quantification', 'report_generation', 'completed', 'failed'
    )),
    started_at TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP WITH TIME ZONE,
    operator_name VARCHAR(255),
    operator_email VARCHAR(255),

    -- Sample metadata
    sample_type VARCHAR(50) DEFAULT 'plasma',
    plasma_volume_ml DECIMAL(10,2),
    collection_date DATE,
    processing_date DATE,

    -- Sequencing metadata
    flowcell_id VARCHAR(100),
    sequencing_kit VARCHAR(100),
    sequencing_start TIMESTAMP WITH TIME ZONE,
    sequencing_duration_hours DECIMAL(10,2),

    -- QC metrics
    total_reads BIGINT,
    total_bases_gb DECIMAL(10,2),
    mean_qscore DECIMAL(5,2),
    median_qscore DECIMAL(5,2),
    n50 INTEGER,
    median_read_length INTEGER,
    active_channels INTEGER,

    -- Host removal metrics
    host_removal_rate DECIMAL(5,2),
    non_host_reads BIGINT,
    non_host_bases_gb DECIMAL(10,2),

    -- Analysis metadata
    pipeline_version VARCHAR(50),
    dorado_version VARCHAR(50),
    kraken2_db_version VARCHAR(50),

    -- S3 paths
    s3_raw_data_path TEXT,
    s3_analysis_path TEXT,
    s3_report_path TEXT,

    -- Error tracking
    error_message TEXT,
    error_phase VARCHAR(50),
    error_timestamp TIMESTAMP WITH TIME ZONE,

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_workflow_status ON workflow_executions(status);
CREATE INDEX idx_workflow_donor ON workflow_executions(donor_id);
CREATE INDEX idx_workflow_dates ON workflow_executions(started_at, completed_at);
CREATE INDEX idx_workflow_sample ON workflow_executions(sample_id);

-- ===== Pathogen Detection Results =====
CREATE TABLE IF NOT EXISTS pathogen_detections (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,
    pathogen_id INTEGER REFERENCES pmda_pathogens(id),
    detection_method VARCHAR(50) NOT NULL CHECK (detection_method IN (
        'kraken2', 'blast', 'assembly', 'perv_specific'
    )),

    -- Detection metrics
    is_detected BOOLEAN NOT NULL DEFAULT false,
    confidence_score DECIMAL(5,4) CHECK (confidence_score BETWEEN 0 AND 1),
    read_count INTEGER,
    unique_kmers INTEGER,
    coverage DECIMAL(10,2),
    coverage_breadth DECIMAL(5,2),
    coverage_depth DECIMAL(10,2),

    -- BLAST specific metrics
    blast_evalue TEXT,
    blast_identity DECIMAL(5,2),
    blast_alignment_length INTEGER,

    -- Assembly specific metrics
    contig_count INTEGER,
    longest_contig INTEGER,
    total_contig_length INTEGER,

    detection_timestamp TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT unique_detection UNIQUE (run_id, pathogen_id, detection_method)
);

CREATE INDEX idx_detection_run ON pathogen_detections(run_id);
CREATE INDEX idx_detection_pathogen ON pathogen_detections(pathogen_id);
CREATE INDEX idx_detection_method ON pathogen_detections(detection_method);
CREATE INDEX idx_detection_positive ON pathogen_detections(is_detected) WHERE is_detected = true;

-- ===== Pathogen Quantification Results =====
CREATE TABLE IF NOT EXISTS pathogen_quantifications (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,
    pathogen_id INTEGER REFERENCES pmda_pathogens(id),

    -- Quantification metrics
    copies_per_ml DECIMAL(15,2),
    copies_per_ml_ci_lower DECIMAL(15,2),
    copies_per_ml_ci_upper DECIMAL(15,2),
    rpm DECIMAL(10,2),  -- Reads per million
    tpm DECIMAL(10,2),  -- Transcripts per million
    fpkm DECIMAL(10,2), -- Fragments per kilobase per million

    -- Normalization factors
    spike_in_recovery DECIMAL(5,2),
    normalization_factor DECIMAL(10,4),

    -- Quality metrics
    quantification_method VARCHAR(50),
    confidence_level VARCHAR(20) CHECK (confidence_level IN ('high', 'medium', 'low')),

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT unique_quantification UNIQUE (run_id, pathogen_id)
);

CREATE INDEX idx_quantification_run ON pathogen_quantifications(run_id);
CREATE INDEX idx_quantification_pathogen ON pathogen_quantifications(pathogen_id);
CREATE INDEX idx_quantification_copies ON pathogen_quantifications(copies_per_ml);

-- ===== PERV Analysis Results =====
CREATE TABLE IF NOT EXISTS perv_analysis (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,

    -- PERV typing
    perv_a_detected BOOLEAN DEFAULT false,
    perv_b_detected BOOLEAN DEFAULT false,
    perv_c_detected BOOLEAN DEFAULT false,
    perv_a_reads INTEGER,
    perv_b_reads INTEGER,
    perv_c_reads INTEGER,

    -- Full-length analysis
    full_length_reads INTEGER,
    full_length_consensus TEXT,
    consensus_quality DECIMAL(5,2),

    -- Recombination detection
    recombinant_detected BOOLEAN DEFAULT false,
    recombination_type VARCHAR(50),
    recombination_breakpoint INTEGER,

    -- Phylogenetic analysis
    phylogenetic_clade VARCHAR(100),
    closest_reference VARCHAR(100),
    distance_to_reference DECIMAL(10,6),

    -- Integration analysis
    integration_sites_detected INTEGER,
    integration_evidence TEXT,

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_perv_run ON perv_analysis(run_id);
CREATE INDEX idx_perv_detected ON perv_analysis(perv_a_detected, perv_b_detected, perv_c_detected);

-- ===== EC2 Instance Tracking =====
CREATE TABLE IF NOT EXISTS ec2_instances (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,
    instance_id VARCHAR(50) NOT NULL,
    phase VARCHAR(50) NOT NULL,
    instance_type VARCHAR(50),

    launch_time TIMESTAMP WITH TIME ZONE NOT NULL,
    termination_time TIMESTAMP WITH TIME ZONE,
    runtime_hours DECIMAL(10,2),

    status VARCHAR(50),
    cost_estimate DECIMAL(10,2),

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_ec2_run ON ec2_instances(run_id);
CREATE INDEX idx_ec2_instance ON ec2_instances(instance_id);
CREATE INDEX idx_ec2_phase ON ec2_instances(phase);

-- ===== Quality Control Metrics =====
CREATE TABLE IF NOT EXISTS qc_metrics (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,
    phase VARCHAR(50) NOT NULL,

    metric_name VARCHAR(100) NOT NULL,
    metric_value DECIMAL(20,4),
    metric_unit VARCHAR(50),

    pass_threshold DECIMAL(20,4),
    passed BOOLEAN,

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT unique_qc_metric UNIQUE (run_id, phase, metric_name)
);

CREATE INDEX idx_qc_run ON qc_metrics(run_id);
CREATE INDEX idx_qc_phase ON qc_metrics(phase);
CREATE INDEX idx_qc_passed ON qc_metrics(passed);

-- ===== Audit Logs =====
CREATE TABLE IF NOT EXISTS audit_logs (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255),
    phase VARCHAR(50),
    action VARCHAR(100) NOT NULL,

    user_id VARCHAR(100),
    user_name VARCHAR(255),
    user_role VARCHAR(50),

    ip_address INET,
    user_agent TEXT,

    old_value JSONB,
    new_value JSONB,
    details JSONB,

    timestamp TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_audit_run ON audit_logs(run_id);
CREATE INDEX idx_audit_timestamp ON audit_logs(timestamp);
CREATE INDEX idx_audit_user ON audit_logs(user_id);
CREATE INDEX idx_audit_action ON audit_logs(action);

-- ===== Notifications Log =====
CREATE TABLE IF NOT EXISTS notifications (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,

    notification_type VARCHAR(50) NOT NULL CHECK (notification_type IN (
        'info', 'warning', 'critical', 'error', 'success'
    )),

    subject TEXT NOT NULL,
    message TEXT NOT NULL,

    recipient_email VARCHAR(255),
    recipient_phone VARCHAR(50),

    sent_via VARCHAR(50) CHECK (sent_via IN ('email', 'sms', 'sns', 'slack')),
    sent_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    delivery_status VARCHAR(50),
    error_message TEXT
);

CREATE INDEX idx_notification_run ON notifications(run_id);
CREATE INDEX idx_notification_type ON notifications(notification_type);
CREATE INDEX idx_notification_sent ON notifications(sent_at);

-- ===== Data Integrity Checks =====
CREATE TABLE IF NOT EXISTS data_integrity (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,

    file_path TEXT NOT NULL,
    file_type VARCHAR(50),
    file_size_bytes BIGINT,

    md5_checksum VARCHAR(32),
    sha256_checksum VARCHAR(64),

    validation_status VARCHAR(50),
    validation_timestamp TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT unique_file_integrity UNIQUE (run_id, file_path)
);

CREATE INDEX idx_integrity_run ON data_integrity(run_id);
CREATE INDEX idx_integrity_status ON data_integrity(validation_status);

-- ===== Analysis Parameters =====
CREATE TABLE IF NOT EXISTS analysis_parameters (
    id SERIAL PRIMARY KEY,
    run_id VARCHAR(255) REFERENCES workflow_executions(run_id) ON DELETE CASCADE,

    parameter_category VARCHAR(100),
    parameter_name VARCHAR(100),
    parameter_value TEXT,

    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT unique_parameter UNIQUE (run_id, parameter_category, parameter_name)
);

CREATE INDEX idx_params_run ON analysis_parameters(run_id);
CREATE INDEX idx_params_category ON analysis_parameters(parameter_category);

-- ===== Functions and Triggers =====

-- Auto-update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER update_workflow_executions_updated_at BEFORE UPDATE
    ON workflow_executions FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_pmda_pathogens_updated_at BEFORE UPDATE
    ON pmda_pathogens FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Function to calculate runtime hours for EC2 instances
CREATE OR REPLACE FUNCTION calculate_runtime_hours()
RETURNS TRIGGER AS $$
BEGIN
    IF NEW.termination_time IS NOT NULL AND NEW.launch_time IS NOT NULL THEN
        NEW.runtime_hours = EXTRACT(EPOCH FROM (NEW.termination_time - NEW.launch_time)) / 3600;
    END IF;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER calculate_ec2_runtime BEFORE INSERT OR UPDATE
    ON ec2_instances FOR EACH ROW
    EXECUTE FUNCTION calculate_runtime_hours();

-- ===== Views for Reporting =====

-- Summary view for pipeline runs
CREATE OR REPLACE VIEW v_pipeline_summary AS
SELECT
    w.run_id,
    w.sample_id,
    w.donor_id,
    w.status,
    w.started_at,
    w.completed_at,
    w.total_reads,
    w.total_bases_gb,
    w.mean_qscore,
    w.host_removal_rate,
    COUNT(DISTINCT CASE WHEN pd.is_detected THEN pd.pathogen_id END) as pathogens_detected,
    COUNT(DISTINCT CASE WHEN pq.copies_per_ml > 0 THEN pq.pathogen_id END) as pathogens_quantified,
    CASE WHEN pa.perv_a_detected OR pa.perv_b_detected OR pa.perv_c_detected THEN true ELSE false END as perv_detected
FROM workflow_executions w
LEFT JOIN pathogen_detections pd ON w.run_id = pd.run_id
LEFT JOIN pathogen_quantifications pq ON w.run_id = pq.run_id
LEFT JOIN perv_analysis pa ON w.run_id = pa.run_id
GROUP BY w.run_id, w.sample_id, w.donor_id, w.status, w.started_at, w.completed_at,
         w.total_reads, w.total_bases_gb, w.mean_qscore, w.host_removal_rate,
         pa.perv_a_detected, pa.perv_b_detected, pa.perv_c_detected;

-- Critical pathogen detection view
CREATE OR REPLACE VIEW v_critical_pathogens AS
SELECT
    pd.run_id,
    w.sample_id,
    w.donor_id,
    p.pathogen_name_en,
    p.pathogen_name_ja,
    p.risk_level,
    pd.detection_method,
    pd.confidence_score,
    pq.copies_per_ml,
    pd.detection_timestamp
FROM pathogen_detections pd
JOIN pmda_pathogens p ON pd.pathogen_id = p.id
JOIN workflow_executions w ON pd.run_id = w.run_id
LEFT JOIN pathogen_quantifications pq ON pd.run_id = pq.run_id AND pd.pathogen_id = pq.pathogen_id
WHERE pd.is_detected = true
  AND (p.risk_level IN ('BSL-3', 'BSL-4') OR p.pathogen_code = 'PERV')
ORDER BY pd.detection_timestamp DESC;

-- ===== Permissions =====
-- These should be executed with appropriate role names after deployment

-- GRANT SELECT ON ALL TABLES IN SCHEMA public TO readonly_role;
-- GRANT INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO readwrite_role;
-- GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO admin_role;
-- GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO readwrite_role;
-- GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO admin_role;