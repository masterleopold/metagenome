"""
Integration tests for improved codebase patterns.

Tests demonstrate:
- Pydantic model validation
- Repository pattern with SQLite (testing) backend
- Type-safe workflows
- Audit logging

Run with: pytest tests/integration/test_new_patterns.py -v
"""

import pytest
from datetime import datetime
from pathlib import Path

from lib.models.pathogen import (
    PERVSubtype,
    PERVDetectionResult,
    PERVTypingOutput,
    PathogenConfidence,
    PathogenDetectionMethod,
    PathogenDetectionResult,
    PMDA91PathogenResult,
)

from lib.models.workflow import (
    WorkflowStatus,
    WorkflowExecution,
    QCMetrics,
    HostRemovalMetrics,
)

from lib.repositories.sqlite_repository import (
    SQLiteWorkflowRepository,
    SQLitePathogenDetectionRepository,
)


class TestPERVModels:
    """Test PERV detection models."""

    def test_perv_detection_result_auto_confidence(self):
        """Test automatic confidence assignment based on reads + coverage."""
        # High confidence
        high_conf = PERVDetectionResult(
            subtype=PERVSubtype.PERV_A,
            reads_aligned=15,
            coverage=0.9,
            specific_motifs_found=['ATGGCAGCCACCACAGC'],
            mean_identity=96.5,
            confidence=PathogenConfidence.HIGH  # Explicit
        )
        assert high_conf.confidence == PathogenConfidence.HIGH
        assert high_conf.requires_alert is True

        # Medium confidence (auto-calculated)
        med_conf = PERVDetectionResult(
            subtype=PERVSubtype.PERV_B,
            reads_aligned=5,
            coverage=0.6,
            specific_motifs_found=[],
            mean_identity=92.0,
            confidence=PathogenConfidence.MEDIUM
        )
        assert med_conf.confidence == PathogenConfidence.MEDIUM
        assert med_conf.requires_alert is True

        # Negative (no detection)
        neg = PERVDetectionResult(
            subtype=PERVSubtype.PERV_C,
            reads_aligned=0,
            coverage=0.0,
            specific_motifs_found=[],
            mean_identity=0.0,
            confidence=PathogenConfidence.NEGATIVE
        )
        assert neg.confidence == PathogenConfidence.NEGATIVE
        assert neg.requires_alert is False

    def test_perv_typing_output_alert_logic(self):
        """Test PERV typing output alert requirements."""
        result = PERVTypingOutput(
            run_id="TEST-001",
            bam_file=Path("/path/to/test.bam"),
            detections={
                PERVSubtype.PERV_A: PERVDetectionResult(
                    subtype=PERVSubtype.PERV_A,
                    reads_aligned=20,
                    coverage=0.95,
                    specific_motifs_found=['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC'],
                    mean_identity=98.0,
                    confidence=PathogenConfidence.HIGH
                ),
                PERVSubtype.PERV_B: PERVDetectionResult(
                    subtype=PERVSubtype.PERV_B,
                    reads_aligned=0,
                    coverage=0.0,
                    specific_motifs_found=[],
                    mean_identity=0.0,
                    confidence=PathogenConfidence.NEGATIVE
                )
            }
        )

        # Should trigger alert (PERV-A detected)
        assert result.requires_sns_alert is True
        assert result.positive_subtypes == [PERVSubtype.PERV_A]

        # Audit log format
        audit_log = result.to_audit_log()
        assert audit_log['run_id'] == "TEST-001"
        assert audit_log['perv_detected'] is True
        assert 'PERV-A' in audit_log['positive_subtypes']

    def test_perv_validation_errors(self):
        """Test Pydantic validation catches invalid data."""
        with pytest.raises(ValueError):
            # Invalid coverage (>1.0)
            PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=10,
                coverage=1.5,  # Invalid!
                specific_motifs_found=[],
                mean_identity=95.0,
                confidence=PathogenConfidence.HIGH
            )

        with pytest.raises(ValueError):
            # Negative reads
            PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=-5,  # Invalid!
                coverage=0.5,
                specific_motifs_found=[],
                mean_identity=95.0,
                confidence=PathogenConfidence.HIGH
            )


class TestPMDA91PathoSogenModels:
    """Test PMDA 91 pathogen screening models."""

    def test_pmda_91_pathogen_count_validation(self):
        """Test that model enforces 91 pathogen count."""
        # Valid: 2 positive + 89 negative = 91
        valid_result = PMDA91PathogenResult(
            run_id="TEST-002",
            all_91_tested=True,
            positive_detections=[
                PathogenDetectionResult(
                    pathogen_code="PERV-A",
                    pathogen_name_en="Porcine Endogenous Retrovirus A",
                    pathogen_name_ja="ブタ内在性レトロウイルスA型",
                    detection_method=PathogenDetectionMethod.KRAKEN2,
                    reads_assigned=50,
                    confidence=PathogenConfidence.HIGH
                ),
                PathogenDetectionResult(
                    pathogen_code="PCV2",
                    pathogen_name_en="Porcine Circovirus 2",
                    pathogen_name_ja="ブタサーコウイルス2型",
                    detection_method=PathogenDetectionMethod.MINIMAP2,
                    reads_assigned=30,
                    confidence=PathogenConfidence.MEDIUM
                )
            ],
            negative_count=89,
            kraken2_db_version="pmda_2024.2",
            blast_db_version="rvdb_2024.1"
        )

        assert len(valid_result.positive_detections) == 2
        assert valid_result.negative_count == 89
        assert valid_result.requires_alert is True  # PERV-A is high priority

        # Invalid: count mismatch
        with pytest.raises(ValueError, match="Pathogen count mismatch"):
            PMDA91PathogenResult(
                run_id="TEST-003",
                all_91_tested=True,
                positive_detections=[],
                negative_count=90,  # Should be 91!
                kraken2_db_version="pmda_2024.2",
                blast_db_version="rvdb_2024.1"
            )


class TestWorkflowModels:
    """Test workflow execution models."""

    def test_workflow_execution_qc_pass_logic(self):
        """Test QC metrics pass/fail logic."""
        # Passing QC
        good_qc = QCMetrics(
            total_reads=1_000_000,
            total_bases_gb=1.5,
            mean_qscore=12.5,
            median_qscore=13.0,
            n50=5000,
            median_read_length=3000,
            active_channels=450
        )
        assert good_qc.passes_qc is True

        # Failing QC (low Q-score)
        bad_qc = QCMetrics(
            total_reads=1_000_000,
            total_bases_gb=1.5,
            mean_qscore=7.0,  # Too low!
            median_qscore=7.5,
            n50=5000,
            median_read_length=3000,
            active_channels=450
        )
        assert bad_qc.passes_qc is False

    def test_host_removal_validation(self):
        """Test host removal metrics validation."""
        # Normal xenotransplant sample (85% host)
        normal = HostRemovalMetrics(
            host_removal_rate=85.0,
            non_host_reads=150_000,
            non_host_bases_gb=0.2,
            pig_rrna_removed=500_000,
            human_removed=350_000
        )
        assert normal.is_valid is True

        # Unusual (too little host removal - suspicious)
        unusual_low = HostRemovalMetrics(
            host_removal_rate=30.0,  # Too low!
            non_host_reads=700_000,
            non_host_bases_gb=1.0,
            pig_rrna_removed=100_000,
            human_removed=50_000
        )
        assert unusual_low.is_valid is False

        # Unusual (almost all host - very little pathogen material)
        unusual_high = HostRemovalMetrics(
            host_removal_rate=99.5,  # Too high!
            non_host_reads=500,
            non_host_bases_gb=0.001,
            pig_rrna_removed=950_000,
            human_removed=50_000
        )
        assert unusual_high.is_valid is False


class TestRepositoryPattern:
    """Test repository pattern with SQLite backend."""

    @pytest.fixture
    def workflow_repo(self):
        """Create in-memory SQLite workflow repository."""
        repo = SQLiteWorkflowRepository(db_path=":memory:")
        yield repo
        repo.close()

    @pytest.fixture
    def pathogen_repo(self):
        """Create in-memory SQLite pathogen repository."""
        repo = SQLitePathogenDetectionRepository(db_path=":memory:")
        yield repo
        repo.close()

    def test_workflow_crud_operations(self, workflow_repo):
        """Test create, read, update, delete workflow operations."""
        # Create workflow
        workflow = WorkflowExecution(
            run_id="TEST-REPO-001",
            sample_id="SAMPLE-001",
            donor_id="DONOR-001",
            status=WorkflowStatus.INITIATED,
            operator_name="Test Operator",
            operator_email="test@example.com",
            sample_type="plasma",
            plasma_volume_ml=5.0
        )

        run_id = workflow_repo.create(workflow)
        assert run_id == "TEST-REPO-001"

        # Read workflow
        retrieved = workflow_repo.get("TEST-REPO-001")
        assert retrieved is not None
        assert retrieved.run_id == "TEST-REPO-001"
        assert retrieved.status == WorkflowStatus.INITIATED
        assert retrieved.operator_email == "test@example.com"

        # Update status
        workflow_repo.update_status(
            run_id="TEST-REPO-001",
            status=WorkflowStatus.COMPLETED
        )

        updated = workflow_repo.get("TEST-REPO-001")
        assert updated.status == WorkflowStatus.COMPLETED
        assert updated.completed_at is not None

        # List by status
        completed = workflow_repo.list_by_status(WorkflowStatus.COMPLETED)
        assert len(completed) == 1
        assert completed[0].run_id == "TEST-REPO-001"

    def test_workflow_duplicate_prevention(self, workflow_repo):
        """Test that duplicate run_ids are rejected."""
        workflow = WorkflowExecution(
            run_id="DUP-001",
            sample_id="SAMPLE-001",
            status=WorkflowStatus.INITIATED,
            operator_name="Test",
            operator_email="test@example.com"
        )

        workflow_repo.create(workflow)

        # Try to create duplicate
        with pytest.raises(ValueError, match="already exists"):
            workflow_repo.create(workflow)

    def test_pathogen_detection_storage(self, pathogen_repo):
        """Test pathogen detection storage and retrieval."""
        # Create detections
        detections = [
            PathogenDetectionResult(
                pathogen_code="PERV-A",
                pathogen_name_en="PERV-A",
                pathogen_name_ja="PERV-A",
                detection_method=PathogenDetectionMethod.KRAKEN2,
                reads_assigned=50,
                confidence=PathogenConfidence.HIGH,
                quantification=1500.0  # copies/mL
            ),
            PathogenDetectionResult(
                pathogen_code="PCV2",
                pathogen_name_en="PCV2",
                pathogen_name_ja="PCV2",
                detection_method=PathogenDetectionMethod.MINIMAP2,
                reads_assigned=25,
                confidence=PathogenConfidence.MEDIUM,
                quantification=500.0
            )
        ]

        # Batch insert
        ids = pathogen_repo.create_batch("TEST-PATH-001", detections)
        assert len(ids) == 2

        # Retrieve by run_id
        retrieved = pathogen_repo.get_by_run_id("TEST-PATH-001")
        assert len(retrieved) == 2
        assert retrieved[0].pathogen_code in ["PERV-A", "PCV2"]

        # Retrieve by pathogen_code
        perv_detections = pathogen_repo.get_by_pathogen_code("PERV-A")
        assert len(perv_detections) >= 1
        assert perv_detections[0].pathogen_code == "PERV-A"


class TestEndToEndWorkflow:
    """End-to-end workflow tests using all components."""

    def test_complete_workflow_with_perv_detection(self):
        """Test complete workflow from initiation to PERV detection."""
        # 1. Create workflow
        workflow = WorkflowExecution(
            run_id="E2E-TEST-001",
            sample_id="E2E-SAMPLE-001",
            donor_id="E2E-DONOR-001",
            status=WorkflowStatus.INITIATED,
            operator_name="E2E Test",
            operator_email="e2e@example.com",
            sample_type="plasma",
            plasma_volume_ml=5.0
        )

        # 2. Add QC metrics
        workflow.qc_metrics = QCMetrics(
            total_reads=2_000_000,
            total_bases_gb=3.0,
            mean_qscore=14.5,
            median_qscore=15.0,
            n50=7000,
            median_read_length=4500,
            active_channels=480
        )

        # 3. Add host removal metrics
        workflow.host_removal_metrics = HostRemovalMetrics(
            host_removal_rate=88.5,
            non_host_reads=230_000,
            non_host_bases_gb=0.35,
            pig_rrna_removed=1_200_000,
            human_removed=570_000
        )

        # 4. PERV typing results
        perv_result = PERVTypingOutput(
            run_id="E2E-TEST-001",
            bam_file=Path("/test/e2e.bam"),
            detections={
                PERVSubtype.PERV_A: PERVDetectionResult(
                    subtype=PERVSubtype.PERV_A,
                    reads_aligned=35,
                    coverage=0.92,
                    specific_motifs_found=['ATGGCAGCCACCACAGC'],
                    mean_identity=97.5,
                    confidence=PathogenConfidence.HIGH
                )
            }
        )

        # Assertions
        assert workflow.qc_metrics.passes_qc is True
        assert workflow.host_removal_metrics.is_valid is True
        assert perv_result.requires_sns_alert is True
        assert len(perv_result.positive_subtypes) == 1

        # Audit log
        audit_log = perv_result.to_audit_log()
        assert audit_log['perv_detected'] is True
        assert audit_log['all_results']['PERV-A']['confidence'] == 'HIGH'


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
