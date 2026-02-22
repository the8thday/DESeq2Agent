"""Report Narrative Agent — generates structured narrative text for the HTML report."""

import logging
from typing import Any, Dict, Optional

from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.models import ReportNarrative
from deseq2_agent.prompts import REPORT_NARRATIVE_PROMPT

logger = logging.getLogger(__name__)


class ReportNarrativeAgent(BaseAgent):
    """Generates report narrative text from all prior analysis results."""

    prompt_template = REPORT_NARRATIVE_PROMPT
    output_model = ReportNarrative

    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        required = ["qc_narrative_input", "de_narrative_input", "pathway_narrative_input"]
        for field in required:
            if field not in payload:
                raise ValueError(f"ReportNarrativeAgent requires '{field}' in input")

        result = dict(payload)
        result.setdefault("study_context", "未提供研究背景")
        result.setdefault("analysis_params", "标准DESeq2参数")
        return result

    def invoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> ReportNarrative:
        return super().invoke(input, config)
