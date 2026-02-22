"""DE Review Agent — interprets differential expression results per contrast."""

import logging
from typing import Any, Dict, Optional

from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.models import DEReviewOutput
from deseq2_agent.prompts import DE_REVIEW_PROMPT

logger = logging.getLogger(__name__)


class DEReviewAgent(BaseAgent):
    """Interprets DE results for all contrasts and produces DEReviewOutput."""

    prompt_template = DE_REVIEW_PROMPT
    output_model = DEReviewOutput

    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        required = ["de_summary", "qc_context"]
        for field in required:
            if field not in payload:
                raise ValueError(f"DEReviewAgent requires '{field}' in input")

        result = dict(payload)
        result.setdefault("padj_threshold", "0.05")
        result.setdefault("lfc_threshold", "1.0")
        return result

    def invoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> DEReviewOutput:
        return super().invoke(input, config)
