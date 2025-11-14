"""
Severity Classification Engine

Evaluates virus detections and assigns severity levels based on configured rules.
Integrates with notification router for alert distribution.
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import yaml


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SeverityEngine:
    """
    Engine for classifying virus detection severity

    Loads rules from YAML configuration and evaluates detections
    against multi-criteria decision matrix.
    """

    def __init__(self, rules_file: Optional[Path] = None):
        """
        Initialize severity engine

        Args:
            rules_file: Path to severity rules YAML file
        """
        if rules_file is None:
            # Default to config/severity_rules.yaml
            rules_file = Path(__file__).parent.parent / "config" / "severity_rules.yaml"

        self.rules_file = rules_file
        self.rules = self._load_rules()

        logger.info(f"Initialized Severity Engine with rules from: {rules_file}")

    def _load_rules(self) -> Dict[str, Any]:
        """
        Load severity rules from YAML configuration

        Returns:
            Rules dictionary
        """
        if not self.rules_file.exists():
            logger.error(f"Rules file not found: {self.rules_file}")
            raise FileNotFoundError(f"Severity rules file not found: {self.rules_file}")

        with open(self.rules_file, 'r', encoding='utf-8') as f:
            rules = yaml.safe_load(f)

        logger.info(f"Loaded severity rules: {len(rules.get('virus_rules', {}))} virus types")
        return rules

    def classify_detection(self, detection: Dict[str, Any]) -> Dict[str, Any]:
        """
        Classify a virus detection and assign severity level

        Args:
            detection: Detection record with virus_type, source, copies_per_ml, etc.

        Returns:
            Detection record with added severity, reason, and notification_config
        """
        virus_type = detection.get('virus_type')
        source = detection.get('source')
        copies_per_ml = detection.get('copies_per_ml', 0)

        if not virus_type:
            logger.warning("Detection missing virus_type, cannot classify")
            detection['severity'] = 'low'
            detection['reason'] = 'Unknown virus type'
            return detection

        # Get virus-specific rules
        virus_rules = self.rules.get('virus_rules', {}).get(virus_type)
        if not virus_rules:
            logger.warning(f"No rules found for virus: {virus_type}")
            detection['severity'] = 'medium'
            detection['reason'] = f'No specific rules for {virus_type}'
            return detection

        # Evaluate thresholds
        severity = 'low'
        reason = 'Default low severity'

        for threshold in virus_rules.get('severity_thresholds', []):
            condition = threshold.get('condition', {})

            if self._evaluate_condition(condition, detection):
                severity = threshold.get('severity', 'low')
                reason = threshold.get('reason', 'Threshold matched')
                logger.info(f"Matched severity rule: {virus_type} -> {severity} ({reason})")
                break

        # Check for compound rules (multiple viruses, external validation, etc.)
        severity, reason = self._apply_compound_rules(detection, severity, reason)

        # Add severity information to detection
        detection['severity'] = severity
        detection['reason'] = reason
        detection['severity_priority'] = self.rules.get('severity_levels', {}).get(severity, {}).get('priority', 4)
        detection['notification_config'] = self._get_notification_config(severity)
        detection['validation_required'] = self._get_validation_requirements(severity)

        logger.info(f"Classification complete: {virus_type} = {severity} (priority {detection['severity_priority']})")

        return detection

    def _evaluate_condition(self, condition: Dict[str, Any], detection: Dict[str, Any]) -> bool:
        """
        Evaluate if a condition matches the detection

        Args:
            condition: Rule condition
            detection: Detection record

        Returns:
            True if condition matches
        """
        # Check source
        if 'source' in condition:
            detection_source = detection.get('source', '')
            required_source = condition['source']

            # Handle internal_pipeline matching
            if required_source == 'internal_pipeline':
                if 'internal_pipeline' not in detection_source:
                    return False
            elif required_source not in detection_source:
                return False

        # Check copies_per_ml threshold
        if 'copies_per_ml' in condition:
            threshold_expr = condition['copies_per_ml']
            copies = detection.get('copies_per_ml', 0)

            if threshold_expr.startswith('>'):
                threshold = float(threshold_expr.strip('> '))
                if copies <= threshold:
                    return False
            elif threshold_expr.startswith('<'):
                threshold = float(threshold_expr.strip('< '))
                if copies >= threshold:
                    return False

        # Check detected flag
        if 'detected' in condition:
            if detection.get('detected', False) != condition['detected']:
                return False

        # Check keyword_match
        if 'keyword_match' in condition:
            if detection.get('metadata', {}).get('keyword_match', False) != condition['keyword_match']:
                return False

        # Check outbreak_reported (for external sources)
        if 'outbreak_reported' in condition:
            if detection.get('metadata', {}).get('outbreak_reported', False) != condition['outbreak_reported']:
                return False

        # Check new_case_report (for academic sources)
        if 'new_case_report' in condition:
            if detection.get('metadata', {}).get('new_case_report', False) != condition['new_case_report']:
                return False

        # All conditions matched
        return True

    def _apply_compound_rules(
        self,
        detection: Dict[str, Any],
        current_severity: str,
        current_reason: str
    ) -> tuple[str, str]:
        """
        Apply compound rules (multiple viruses, external validation, etc.)

        Args:
            detection: Detection record
            current_severity: Current severity level
            current_reason: Current reason

        Returns:
            Updated (severity, reason) tuple
        """
        compound_rules = self.rules.get('compound_rules', {})

        # Check external validation
        if detection.get('external_confirmation', False):
            validation_rules = compound_rules.get('external_validation', [])
            for rule in validation_rules:
                if self._evaluate_condition(rule.get('condition', {}), detection):
                    # Increase severity by one level
                    new_severity = self._increase_severity(current_severity)
                    return new_severity, f"{current_reason} + External validation"

        # Note: Multiple virus detection would need context of other detections
        # This would be handled by the notification router

        return current_severity, current_reason

    def _increase_severity(self, current: str) -> str:
        """
        Increase severity by one level

        Args:
            current: Current severity level

        Returns:
            Increased severity level
        """
        severity_order = ['low', 'medium', 'high', 'critical']
        try:
            current_index = severity_order.index(current)
            if current_index < len(severity_order) - 1:
                return severity_order[current_index + 1]
        except ValueError:
            pass

        return current

    def _get_notification_config(self, severity: str) -> Dict[str, Any]:
        """
        Get notification configuration for severity level

        Args:
            severity: Severity level

        Returns:
            Notification configuration
        """
        alert_routing = self.rules.get('alert_routing', {})
        return alert_routing.get(severity, {})

    def _get_validation_requirements(self, severity: str) -> List[str]:
        """
        Get validation requirements for severity level

        Args:
            severity: Severity level

        Returns:
            List of validation steps required
        """
        validation = self.rules.get('validation', {})
        requirements_key = f"{severity}_detections"
        return validation.get(requirements_key, [])

    def batch_classify(self, detections: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Classify multiple detections and check for compound rules

        Args:
            detections: List of detection records

        Returns:
            List of classified detections
        """
        logger.info(f"Batch classifying {len(detections)} detections")

        # First pass: individual classification
        classified = [self.classify_detection(d) for d in detections]

        # Second pass: compound rules (multiple viruses)
        critical_count = sum(1 for d in classified if d.get('severity') == 'critical')
        high_count = sum(1 for d in classified if d.get('severity') == 'high')
        total_viruses = len(set(d.get('virus_type') for d in classified))

        compound_rules = self.rules.get('compound_rules', {}).get('multiple_viruses', [])

        for rule in compound_rules:
            condition = rule.get('condition', {})
            action = rule.get('action', {})

            # Check condition
            viruses_detected_threshold = int(condition.get('viruses_detected', '0').strip('>= '))
            any_critical = condition.get('any_severity') == 'critical'
            all_high = condition.get('all_severity') == 'high'

            if total_viruses >= viruses_detected_threshold:
                if (any_critical and critical_count > 0) or \
                   (all_high and high_count == total_viruses):
                    # Apply action to all detections
                    for detection in classified:
                        detection['severity'] = action.get('severity', detection['severity'])
                        detection['reason'] += f" + {action.get('reason', '')}"
                        detection['compound_rule_applied'] = True

                    logger.warning(f"Compound rule applied: {action.get('reason')}")

        return classified

    def should_deduplicate(
        self,
        new_detection: Dict[str, Any],
        recent_detections: List[Dict[str, Any]]
    ) -> bool:
        """
        Check if detection should be deduplicated

        Args:
            new_detection: New detection to check
            recent_detections: Recent detections within time window

        Returns:
            True if should deduplicate (suppress alert)
        """
        dedup_rules = self.rules.get('deduplication', {})
        time_window = dedup_rules.get('time_window', 3600)  # seconds
        group_by = dedup_rules.get('group_by', [])

        current_time = datetime.now()

        for recent in recent_detections:
            # Check if within time window
            recent_time = datetime.fromisoformat(recent.get('timestamp', current_time.isoformat()))
            if (current_time - recent_time).total_seconds() > time_window:
                continue

            # Check if matches group_by criteria
            matches = all(
                new_detection.get(field) == recent.get(field)
                for field in group_by
            )

            if matches:
                logger.info(f"Deduplicating detection: {new_detection.get('virus_type')} "
                          f"matches recent detection within {time_window}s")
                return True

        return False


if __name__ == "__main__":
    # Example usage
    engine = SeverityEngine()

    # Test detection
    test_detection = {
        'virus_type': 'spumavirus',
        'source': 'internal_pipeline_kraken2',
        'copies_per_ml': 600,
        'detected': True,
        'timestamp': datetime.now().isoformat()
    }

    classified = engine.classify_detection(test_detection)
    print(f"Classification: {classified['severity']}")
    print(f"Reason: {classified['reason']}")
    print(f"Notification: {classified['notification_config']}")
