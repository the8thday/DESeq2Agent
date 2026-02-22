"""DESeq2 Pipeline Orchestrator.

Coordinates R subprocess execution and LLM agents for end-to-end RNA-seq analysis:
counts.csv + metadata.csv + contrasts → HTML report
"""

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from langchain_core.language_models import BaseChatModel

from deseq2_agent.agents.qc_agent import ConsoleIO, QCReviewAgent
from deseq2_agent.agents.de_agent import DEReviewAgent
from deseq2_agent.agents.pathway_agent import PathwayReviewAgent
from deseq2_agent.agents.report_agent import ReportNarrativeAgent
from deseq2_agent.models import QCDecision, DEReviewOutput, PathwayReviewOutput, ReportNarrative
from deseq2_agent.report.html_builder import HTMLReportBuilder
from deseq2_agent.runner import RScriptRunner, RScriptError

logger = logging.getLogger(__name__)


@dataclass
class ContrastConfig:
    """Single contrast specification."""
    name: str
    variable: str
    treatment: str
    control: str

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "variable": self.variable,
            "treatment": self.treatment,
            "control": self.control,
        }


@dataclass
class PipelineConfig:
    """Full pipeline configuration."""
    counts_file: str
    metadata_file: str
    contrasts: List[ContrastConfig]
    species: str = "human"
    output_dir: str = "./results/"
    mode: str = "auto"                    # "auto" or "interactive"
    min_count_threshold: int = 10
    min_samples_fraction: float = 0.5
    padj_threshold: float = 0.05
    lfc_threshold: float = 1.0

    def to_r_config(self, samples_to_remove: Optional[List[str]] = None) -> dict:
        """Serialize to the config JSON that R scripts consume."""
        return {
            "counts_file": str(Path(self.counts_file).resolve()),
            "metadata_file": str(Path(self.metadata_file).resolve()),
            "contrasts": [c.to_dict() for c in self.contrasts],
            "species": self.species,
            "output_dir": str(Path(self.output_dir).resolve()) + "/",
            "min_count_threshold": self.min_count_threshold,
            "min_samples_fraction": self.min_samples_fraction,
            "padj_threshold": self.padj_threshold,
            "lfc_threshold": self.lfc_threshold,
            "samples_to_remove": samples_to_remove or [],
        }


@dataclass
class PipelineResults:
    """All outputs from a completed pipeline run."""
    qc_decision: Optional[QCDecision] = None
    de_review: Optional[DEReviewOutput] = None
    pathway_review: Optional[PathwayReviewOutput] = None
    narrative: Optional[ReportNarrative] = None
    report_path: Optional[str] = None
    output_dir: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "qc_decision": self.qc_decision.model_dump() if self.qc_decision else None,
            "de_review": self.de_review.model_dump() if self.de_review else None,
            "pathway_review": self.pathway_review.model_dump() if self.pathway_review else None,
            "narrative": self.narrative.model_dump() if self.narrative else None,
            "report_path": self.report_path,
            "output_dir": self.output_dir,
        }


class DESeq2Pipeline:
    """End-to-end RNA-seq analysis pipeline.

    Runs R scripts for computation and LLM agents for interpretation.
    Produces a self-contained HTML report.
    """

    def __init__(self, llm: BaseChatModel, mode: str = "auto"):
        """Initialize pipeline.

        Args:
            llm: LangChain chat model (ChatOpenAI or compatible)
            mode: "auto" (LLM decides outliers) or "interactive" (user confirms)
        """
        if mode not in ("auto", "interactive"):
            raise ValueError(f"mode must be 'auto' or 'interactive', got '{mode}'")
        self.llm = llm
        self.mode = mode

        self.qc_agent = QCReviewAgent(llm)
        self.de_agent = DEReviewAgent(llm)
        self.pathway_agent = PathwayReviewAgent(llm)
        self.report_agent = ReportNarrativeAgent(llm)
        self.html_builder = HTMLReportBuilder()

    def run(
        self,
        counts_file: str,
        metadata_file: str,
        contrasts: List[dict],
        species: str = "human",
        output_dir: str = "./results/",
    ) -> PipelineResults:
        """Execute the full pipeline.

        Args:
            counts_file: Path to counts CSV (rows=genes, cols=samples, first col=gene IDs)
            metadata_file: Path to metadata CSV (rows=samples, cols=covariates)
            contrasts: List of contrast dicts with keys: name, variable, treatment, control
            species: "human" or "mouse"
            output_dir: Directory to write all outputs

        Returns:
            PipelineResults with all agent outputs and path to HTML report
        """
        # --- Setup ---
        output_dir = str(Path(output_dir).resolve())
        os.makedirs(output_dir, exist_ok=True)

        contrast_objs = [ContrastConfig(**c) for c in contrasts]
        config = PipelineConfig(
            counts_file=counts_file,
            metadata_file=metadata_file,
            contrasts=contrast_objs,
            species=species,
            output_dir=output_dir,
            mode=self.mode,
        )

        runner = RScriptRunner(output_dir)

        if not runner.check_r_available():
            raise RuntimeError(
                "Rscript not found on PATH. Please install R (>=4.2) and ensure "
                "'Rscript' is accessible."
            )

        results = PipelineResults(output_dir=output_dir)

        # =========================================================
        # Step 1: Data Preparation
        # =========================================================
        logger.info("=== Step 1: Data Preparation ===")
        runner.write_config(config.to_r_config())
        runner.run("01_data_prep.R")

        missing = runner.validate_outputs(["data_prep_summary.json", "counts_ensembl.csv",
                                           "id_mapping.csv", "dds.rds"])
        if missing:
            raise RuntimeError(f"Data prep missing outputs: {missing}")

        with open(os.path.join(output_dir, "data_prep_summary.json"), "r") as f:
            data_prep_summary = json.load(f)
        logger.info(f"Data prep: {data_prep_summary['n_genes_filtered']} genes, "
                    f"{data_prep_summary['n_samples']} samples, "
                    f"ID type: {data_prep_summary['gene_id_type']}")

        # =========================================================
        # Step 2: QC Analysis (R)
        # =========================================================
        logger.info("=== Step 2: QC Analysis (R) ===")
        runner.run("02_qc_analysis.R")

        missing = runner.validate_outputs(["qc_metrics.json"])
        if missing:
            raise RuntimeError(f"QC analysis missing outputs: {missing}")

        with open(os.path.join(output_dir, "qc_metrics.json"), "r") as f:
            qc_metrics = json.load(f)

        # =========================================================
        # Step 3: QC Review (LLM)
        # =========================================================
        logger.info("=== Step 3: QC Review (LLM) ===")

        import pandas as pd
        metadata_df = pd.read_csv(metadata_file, index_col=0)
        metadata_summary = metadata_df.to_string()

        # R's c() on list-of-lists may serialize as a JSON object with numeric keys
        # rather than a JSON array — normalize to list either way
        outlier_flags_raw = qc_metrics.get("outlier_flags", [])
        if isinstance(outlier_flags_raw, dict):
            outlier_flags = list(outlier_flags_raw.values())
        else:
            outlier_flags = outlier_flags_raw

        outlier_summary = (
            "\n".join(
                f"- {f['sample']}: {f['metric']} = {f['value']:.2f} "
                f"(IQR bounds: [{f['threshold']['lower']:.2f}, {f['threshold']['upper']:.2f}])"
                for f in outlier_flags
                if isinstance(f, dict) and f.get("sample") and f.get("sample") != "pairwise"
            ) or "未检测到明显离群样本"
        )

        lib_summary = "\n".join(
            f"  {s}: {v/1e6:.2f}M reads"
            for s, v in qc_metrics.get("library_sizes", {}).items()
        )
        pca_var = qc_metrics.get("pca_variance_pct", {})
        pca_text = (f"PC1={list(pca_var.values())[0]}%, PC2={list(pca_var.values())[1]}%"
                    if len(pca_var) >= 2 else "N/A")

        qc_metrics_summary = (
            f"样本数: {qc_metrics.get('n_samples', 'N/A')}\n"
            f"过滤后基因数: {qc_metrics.get('n_genes_filtered', 'N/A')}\n"
            f"文库大小:\n{lib_summary}\n"
            f"PCA方差解释: {pca_text}\n"
            f"检测到的离群样本: {', '.join(qc_metrics.get('outlier_samples_flagged', [])) or '无'}"
        )

        qc_input = {
            "qc_metrics_summary": qc_metrics_summary,
            "metadata_summary": metadata_summary,
            "outlier_flags_summary": outlier_summary,
            "additional_context": f"物种: {species}, 对比组: {[c.name for c in contrast_objs]}",
        }

        if self.mode == "interactive":
            qc_decision = self.qc_agent.invoke_interactive(qc_input, console_io=ConsoleIO())
        else:
            qc_decision = self.qc_agent.invoke(qc_input)

        results.qc_decision = qc_decision
        logger.info(f"QC decision: usability={qc_decision.data_usability}, "
                    f"removing={qc_decision.samples_to_remove}")

        if qc_decision.data_usability == "stop":
            logger.warning("QC agent recommends stopping analysis. Data quality too poor.")

        # =========================================================
        # Step 4: Differential Expression (R)
        # =========================================================
        logger.info("=== Step 4: Differential Expression (R) ===")
        runner.write_config(config.to_r_config(
            samples_to_remove=qc_decision.samples_to_remove
        ))
        runner.run("03_de_analysis.R")

        missing = runner.validate_outputs(["de_summary.json"])
        if missing:
            raise RuntimeError(f"DE analysis missing outputs: {missing}")

        with open(os.path.join(output_dir, "de_summary.json"), "r") as f:
            de_summary = json.load(f)

        # =========================================================
        # Step 5: DE Review (LLM)
        # =========================================================
        logger.info("=== Step 5: DE Review (LLM) ===")

        de_summary_text = json.dumps(de_summary, indent=2, ensure_ascii=False)[:8000]
        qc_context = (
            f"QC摘要: {qc_decision.quality_summary}\n"
            f"移除样本: {', '.join(qc_decision.samples_to_remove) or '无'}\n"
            f"数据可用性: {qc_decision.data_usability}"
        )

        de_review = self.de_agent.invoke({
            "de_summary": de_summary_text,
            "qc_context": qc_context,
            "padj_threshold": str(config.padj_threshold),
            "lfc_threshold": str(config.lfc_threshold),
        })
        results.de_review = de_review
        logger.info(f"DE review complete for {len(de_review.contrasts)} contrast(s)")

        # =========================================================
        # Step 6: Enrichment Analysis (R)
        # =========================================================
        logger.info("=== Step 6: Enrichment Analysis (R) ===")
        try:
            runner.run("04_enrichment.R")
        except RScriptError as e:
            logger.warning(f"Enrichment analysis failed (non-fatal): {e}")

        enrichment_results = []
        for ct in contrast_objs:
            enrich_path = os.path.join(output_dir, f"enrichment_results_{ct.name}.json")
            if os.path.exists(enrich_path):
                with open(enrich_path, "r") as f:
                    enrichment_results.append(json.load(f))
            else:
                logger.warning(f"Enrichment results not found for {ct.name}")
                enrichment_results.append({
                    "contrast_name": ct.name,
                    "ora_skipped": True,
                    "gsea_go": None,
                    "gsea_kegg": None,
                    "ora_go": None,
                    "ora_kegg": None,
                    "plots": {},
                })

        # =========================================================
        # Step 7: Pathway Review (LLM)
        # =========================================================
        logger.info("=== Step 7: Pathway Review (LLM) ===")

        def summarize_enrichment(enrich: dict) -> str:
            lines = [f"对比度: {enrich['contrast_name']}"]
            if enrich.get("ora_skipped"):
                lines.append("  ORA已跳过（显著基因<10）")
            for key in ("gsea_go", "gsea_kegg", "ora_go", "ora_kegg"):
                rows = enrich.get(key)
                if rows:
                    top3 = [r.get("Description", r.get("ID", "")) for r in rows[:3]]
                    lines.append(f"  {key.upper()}: {', '.join(top3)}")
                else:
                    lines.append(f"  {key.upper()}: 无显著结果")
            return "\n".join(lines)

        enrichment_summary_text = "\n\n".join(
            summarize_enrichment(e) for e in enrichment_results
        )
        de_context_for_pathway = (
            f"DE结果概述:\n" +
            "\n".join(
                f"  {ct_name}: 上调{v.get('n_sig_up',0)}个，下调{v.get('n_sig_down',0)}个基因"
                for ct_name, v in de_summary.items()
            )
        )

        pathway_review = self.pathway_agent.invoke({
            "enrichment_summary": enrichment_summary_text,
            "de_context": de_context_for_pathway,
            "biological_context": f"物种: {species}",
        })
        results.pathway_review = pathway_review
        logger.info("Pathway review complete")

        # =========================================================
        # Step 8: Report Narrative (LLM)
        # =========================================================
        logger.info("=== Step 8: Report Narrative (LLM) ===")

        qc_narrative_input = (
            f"QC决策: {qc_decision.model_dump_json(indent=2)}"
        )
        de_narrative_input = (
            f"DE解读: {de_review.model_dump_json(indent=2)}"
        )
        pathway_narrative_input = (
            f"通路解读: {pathway_review.model_dump_json(indent=2)}"
        )
        analysis_params = (
            f"padj阈值={config.padj_threshold}, log2FC阈值={config.lfc_threshold}, "
            f"物种={species}, 对比组={'、'.join(c.name for c in contrast_objs)}, "
            f"DESeq2 + apeglm收缩, clusterProfiler富集分析"
        )

        narrative = self.report_agent.invoke({
            "study_context": f"RNA-seq差异表达分析，物种: {species}，对比组: {[c.name for c in contrast_objs]}",
            "qc_narrative_input": qc_narrative_input,
            "de_narrative_input": de_narrative_input,
            "pathway_narrative_input": pathway_narrative_input,
            "analysis_params": analysis_params,
        })
        results.narrative = narrative
        logger.info("Report narrative complete")

        # =========================================================
        # Step 9: HTML Report
        # =========================================================
        logger.info("=== Step 9: Generating HTML Report ===")
        report_path = os.path.join(output_dir, "report.html")

        results.report_path = self.html_builder.build(
            narrative=narrative,
            qc_decision=qc_decision,
            de_summary=de_summary,
            enrichment_results=enrichment_results,
            qc_metrics=qc_metrics,
            contrasts=[c.to_dict() for c in contrast_objs],
            species=species,
            output_path=report_path,
        )

        logger.info(f"Pipeline complete! Report: {results.report_path}")
        return results
