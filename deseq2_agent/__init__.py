"""DESeq2Agent — end-to-end RNA-seq differential expression analysis with LLM interpretation.

Combines R subprocesses (DESeq2, clusterProfiler) with LangChain LLM agents to:
1. Run DESeq2 QC, differential expression, and enrichment analysis
2. Interpret results with specialized AI agents
3. Generate a self-contained HTML report

Quick start::

    from deseq2_agent import DESeq2Pipeline
    from deseq2_agent.config import get_llm

    llm = get_llm()
    pipeline = DESeq2Pipeline(llm=llm, mode="auto")
    results = pipeline.run(
        counts_file="counts.csv",
        metadata_file="metadata.csv",
        contrasts=[{"name": "TreatvsCtrl", "variable": "condition",
                    "treatment": "Treatment", "control": "Control"}],
        species="human",
        output_dir="./results/",
    )
    print(f"Report: {results.report_path}")
"""

from deseq2_agent.config import LLMConfig, get_llm
from deseq2_agent.models import (
    DesignDecision,
    ContrastSuggestion,
    QCDecision,
    OutlierSample,
    DEReviewOutput,
    DEContrastInterpretation,
    PathwayReviewOutput,
    PathwayContrastInterpretation,
    PathwayTheme,
    ReportNarrative,
)
from deseq2_agent.pipeline import DESeq2Pipeline, PipelineConfig, PipelineResults, ContrastConfig
from deseq2_agent.runner import RScriptRunner, RScriptResult, RScriptError
from deseq2_agent.agents import (
    BaseAgent,
    DesignDetectionAgent,
    QCReviewAgent,
    DEReviewAgent,
    PathwayReviewAgent,
    ReportNarrativeAgent,
)
from deseq2_agent.report import HTMLReportBuilder

__version__ = "2.0.0"

__all__ = [
    # Config
    "LLMConfig",
    "get_llm",
    # Pipeline
    "DESeq2Pipeline",
    "PipelineConfig",
    "PipelineResults",
    "ContrastConfig",
    # Models
    "DesignDecision",
    "ContrastSuggestion",
    "QCDecision",
    "OutlierSample",
    "DEReviewOutput",
    "DEContrastInterpretation",
    "PathwayReviewOutput",
    "PathwayContrastInterpretation",
    "PathwayTheme",
    "ReportNarrative",
    # Runner
    "RScriptRunner",
    "RScriptResult",
    "RScriptError",
    # Agents
    "BaseAgent",
    "DesignDetectionAgent",
    "QCReviewAgent",
    "DEReviewAgent",
    "PathwayReviewAgent",
    "ReportNarrativeAgent",
    # Report
    "HTMLReportBuilder",
]
