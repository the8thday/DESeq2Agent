"""Pydantic v2 output models for DESeq2Agent pipeline."""

from typing import Dict, List, Optional
from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# QC Models
# ---------------------------------------------------------------------------

class OutlierSample(BaseModel):
    sample_id: str = Field(description="Sample identifier")
    reason: str = Field(description="Reason why this sample is flagged as outlier")
    confidence: str = Field(description="Confidence level: 'high', 'medium', or 'low'")


class QCDecision(BaseModel):
    quality_summary: str = Field(description="Overall data quality summary in 1-3 sentences")
    outlier_samples: List[OutlierSample] = Field(
        default_factory=list,
        description="List of samples flagged as potential outliers with reasons"
    )
    samples_to_remove: List[str] = Field(
        default_factory=list,
        description="Sample IDs that should be removed before DE analysis"
    )
    batch_effect_present: bool = Field(description="Whether a batch effect is detected")
    batch_effect_severity: str = Field(
        description="Severity of batch effect: 'none', 'mild', 'moderate', or 'severe'"
    )
    pca_interpretation: str = Field(
        description="Interpretation of PCA plot: sample clustering, separation by condition, outliers"
    )
    data_usability: str = Field(
        description="Overall usability assessment: 'proceed', 'proceed_with_caution', or 'stop'"
    )
    reasoning: str = Field(
        description="Detailed reasoning for the QC decisions made"
    )
    qc_flags: List[str] = Field(
        default_factory=list,
        description="List of specific QC concerns or warnings (e.g., 'low library size in sample X')"
    )


# ---------------------------------------------------------------------------
# DE Models
# ---------------------------------------------------------------------------

class DEContrastInterpretation(BaseModel):
    contrast_name: str = Field(description="Name of the contrast")
    overall_assessment: str = Field(
        description="Overall assessment of DE results for this contrast"
    )
    biological_coherence: str = Field(
        description="Biological coherence of results: 'strong', 'moderate', 'weak', or 'unclear'"
    )
    key_observations: List[str] = Field(
        description="Key biological observations from the DE results"
    )
    potential_concerns: List[str] = Field(
        default_factory=list,
        description="Potential technical or biological concerns about the results"
    )
    top_upregulated_commentary: str = Field(
        description="Commentary on the top upregulated genes and their biological significance"
    )
    top_downregulated_commentary: str = Field(
        description="Commentary on the top downregulated genes and their biological significance"
    )


class DEReviewOutput(BaseModel):
    contrasts: List[DEContrastInterpretation] = Field(
        description="Interpretation for each contrast"
    )
    cross_contrast_observations: Optional[str] = Field(
        default=None,
        description="Observations comparing results across multiple contrasts (if applicable)"
    )
    overall_data_quality_comment: str = Field(
        description="Comment on overall data quality based on DE results"
    )


# ---------------------------------------------------------------------------
# Pathway Models
# ---------------------------------------------------------------------------

class PathwayTheme(BaseModel):
    theme_name: str = Field(description="Name of the biological theme or process")
    supporting_pathways: List[str] = Field(
        description="Pathway names that support this theme"
    )
    direction: str = Field(
        description="Direction of enrichment: 'up', 'down', or 'mixed'"
    )


class PathwayContrastInterpretation(BaseModel):
    contrast_name: str = Field(description="Name of the contrast")
    mechanistic_summary: str = Field(
        description="Mechanistic summary of pathway findings for this contrast"
    )
    key_themes: List[PathwayTheme] = Field(
        description="Key biological themes identified from pathway analysis"
    )
    unexpected_findings: List[str] = Field(
        default_factory=list,
        description="Unexpected or surprising pathway findings"
    )
    missing_expected_pathways: List[str] = Field(
        default_factory=list,
        description="Pathways that would be expected biologically but are missing from results"
    )


class PathwayReviewOutput(BaseModel):
    contrasts: List[PathwayContrastInterpretation] = Field(
        description="Pathway interpretation for each contrast"
    )
    cross_contrast_pathway_observations: Optional[str] = Field(
        default=None,
        description="Cross-contrast observations about shared or divergent pathway themes"
    )


# ---------------------------------------------------------------------------
# Report Models
# ---------------------------------------------------------------------------

class ReportNarrative(BaseModel):
    executive_summary: str = Field(
        description="Executive summary of the entire analysis (3-5 sentences)"
    )
    methods_paragraph: str = Field(
        description="Methods paragraph describing the analysis pipeline"
    )
    qc_narrative: str = Field(
        description="Narrative description of QC results and decisions"
    )
    de_narrative: str = Field(
        description="Narrative description of differential expression results"
    )
    pathway_narrative: str = Field(
        description="Narrative description of pathway enrichment results"
    )
    key_findings: List[Dict[str, str]] = Field(
        description="List of key findings, each with 'finding' and 'priority' (high/medium/low) keys"
    )
    limitations: List[str] = Field(
        description="List of limitations and caveats of the analysis"
    )
    conclusions: str = Field(
        description="Conclusions section summarizing the biological interpretation"
    )
