"""QC Review Agent — interprets QC metrics and decides on outlier removal."""

import json
import logging
from typing import Any, Dict, Optional

from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.models import QCDecision
from deseq2_agent.prompts import QC_REVIEW_PROMPT

logger = logging.getLogger(__name__)


class ConsoleIO:
    """Simple console interface for interactive mode. Swappable for testing."""

    def display(self, message: str) -> None:
        print(message)

    def ask(self, question: str) -> str:
        return input(question).strip()


class QCReviewAgent(BaseAgent):
    """Interprets QC metrics and produces a QCDecision (which samples to remove)."""

    prompt_template = QC_REVIEW_PROMPT
    output_model = QCDecision

    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        required = ["qc_metrics_summary", "metadata_summary"]
        for field in required:
            if field not in payload:
                raise ValueError(f"QCReviewAgent requires '{field}' in input")

        result = dict(payload)
        result.setdefault("outlier_flags_summary", "未检测到明显离群样本")
        result.setdefault("additional_context", "无")
        return result

    def invoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> QCDecision:
        return super().invoke(input, config)

    def invoke_interactive(
        self,
        input: Dict[str, Any],
        console_io: Optional[ConsoleIO] = None,
        config: Optional[Any] = None,
    ) -> QCDecision:
        """Run QC review in interactive mode: user confirms each outlier candidate.

        Args:
            input: Input payload (same as invoke())
            console_io: ConsoleIO instance (defaults to real console)
            config: Optional LangChain config

        Returns:
            QCDecision with user-confirmed samples_to_remove
        """
        if console_io is None:
            console_io = ConsoleIO()

        decision: QCDecision = super().invoke(input, config)

        if not decision.outlier_samples:
            console_io.display("\n[QC] No outlier samples detected by LLM.")
            return decision

        console_io.display(f"\n[QC] LLM identified {len(decision.outlier_samples)} potential outlier(s):")
        console_io.display(f"     Quality summary: {decision.quality_summary}\n")

        user_confirmed_remove = []
        for outlier in decision.outlier_samples:
            console_io.display(
                f"  Sample: {outlier.sample_id}\n"
                f"  Reason: {outlier.reason}\n"
                f"  Confidence: {outlier.confidence}"
            )
            answer = console_io.ask(f"  Remove sample '{outlier.sample_id}'? [y/N]: ")
            if answer.lower() in ("y", "yes"):
                user_confirmed_remove.append(outlier.sample_id)
                console_io.display(f"  -> Will REMOVE {outlier.sample_id}\n")
            else:
                console_io.display(f"  -> Keeping {outlier.sample_id}\n")

        # Return updated decision with user-confirmed removals
        return QCDecision(
            quality_summary=decision.quality_summary,
            outlier_samples=decision.outlier_samples,
            samples_to_remove=user_confirmed_remove,
            batch_effect_present=decision.batch_effect_present,
            batch_effect_severity=decision.batch_effect_severity,
            pca_interpretation=decision.pca_interpretation,
            data_usability=decision.data_usability,
            reasoning=decision.reasoning,
            qc_flags=decision.qc_flags,
        )
