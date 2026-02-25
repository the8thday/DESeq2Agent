"""DesignDetectionAgent — auto-detects experimental design from metadata."""

import json
from typing import Any, Dict

from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.models import DesignDecision
from deseq2_agent.prompts import DESIGN_DETECTION_PROMPT


class DesignDetectionAgent(BaseAgent):
    """Analyses sample metadata and recommends DESeq2 design + contrasts.

    Input keys:
        metadata_summary (str): JSON string built by pipeline._build_metadata_summary()

    Output:
        DesignDecision Pydantic model
    """

    prompt_template = DESIGN_DETECTION_PROMPT
    output_model = DesignDecision

    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        if "metadata_summary" not in payload:
            raise ValueError("DesignDetectionAgent requires 'metadata_summary'")
        # Ensure it's a plain string (not already a dict)
        if not isinstance(payload["metadata_summary"], str):
            payload["metadata_summary"] = json.dumps(
                payload["metadata_summary"], indent=2, ensure_ascii=False
            )
        return payload
