"""Pathway Review Agent — interprets GSEA/ORA enrichment results."""

import logging
from typing import Any, Dict, Optional

from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.models import PathwayReviewOutput
from deseq2_agent.prompts import PATHWAY_REVIEW_PROMPT

logger = logging.getLogger(__name__)


class PathwayReviewAgent(BaseAgent):
    """Interprets pathway enrichment results for all contrasts."""

    prompt_template = PATHWAY_REVIEW_PROMPT
    output_model = PathwayReviewOutput

    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        required = ["enrichment_summary", "de_context"]
        for field in required:
            if field not in payload:
                raise ValueError(f"PathwayReviewAgent requires '{field}' in input")

        result = dict(payload)
        result.setdefault("biological_context", "未提供具体生物学背景")
        return result

    def invoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> PathwayReviewOutput:
        return super().invoke(input, config)
