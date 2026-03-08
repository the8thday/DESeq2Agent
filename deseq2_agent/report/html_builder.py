"""HTML Report Builder — assembles self-contained HTML from plots + LLM narratives."""

import base64
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from jinja2 import Environment, FileSystemLoader

logger = logging.getLogger(__name__)

TEMPLATES_DIR = Path(__file__).parent / "templates"


class HTMLReportBuilder:
    """Builds a self-contained HTML report by embedding plots as base64 images."""

    def __init__(self):
        self.env = Environment(
            loader=FileSystemLoader(str(TEMPLATES_DIR)),
            autoescape=True,
        )

    @staticmethod
    def _pdf_to_png(pdf_path: Path) -> Optional[Path]:
        """Convert a single-page PDF to a temporary PNG using sips (macOS)."""
        png_path = pdf_path.with_suffix(".png")
        if png_path.exists() and png_path.stat().st_mtime >= pdf_path.stat().st_mtime:
            return png_path
        try:
            result = subprocess.run(
                ["sips", "-s", "format", "png",
                 "-s", "dpiWidth", "200", "-s", "dpiHeight", "200",
                 str(pdf_path), "--out", str(png_path)],
                capture_output=True, text=True, timeout=30,
            )
            if result.returncode == 0 and png_path.exists():
                return png_path
            logger.warning(f"sips PDF→PNG failed for {pdf_path}: {result.stderr}")
        except FileNotFoundError:
            logger.warning("sips not found (non-macOS?). PDF plots will not be embedded.")
        except Exception as e:
            logger.warning(f"PDF→PNG conversion failed for {pdf_path}: {e}")
        return None

    def _encode_plot_base64(self, plot_path: str) -> Optional[str]:
        """Return data URI for a plot file (PNG or PDF), or None if missing."""
        path = Path(plot_path)
        if not path.exists():
            logger.warning(f"Plot file not found: {plot_path}")
            return None
        try:
            # If PDF, convert to PNG first for HTML embedding
            if path.suffix.lower() == ".pdf":
                png_path = self._pdf_to_png(path)
                if png_path is None:
                    return None
                path = png_path

            with open(path, "rb") as f:
                data = base64.b64encode(f.read()).decode("ascii")
            return f"data:image/png;base64,{data}"
        except Exception as e:
            logger.warning(f"Failed to encode {plot_path}: {e}")
            return None

    def _encode_plots(self, plot_paths: Dict[str, str]) -> Dict[str, Optional[str]]:
        """Encode a dict of {name: path} plots to base64 data URIs."""
        return {name: self._encode_plot_base64(path) for name, path in plot_paths.items()}

    def build(
        self,
        narrative,           # ReportNarrative Pydantic model
        qc_decision,         # QCDecision Pydantic model
        de_summary: Dict,    # Parsed de_summary.json
        enrichment_results: List[Dict],  # List of parsed enrichment JSONs
        qc_metrics: Dict,    # Parsed qc_metrics.json
        contrasts: List[Dict],
        species: str,
        output_path: str,
        edger_sensitivity: Optional[Dict] = None,
    ) -> str:
        """Build and write the HTML report.

        Args:
            narrative: ReportNarrative from ReportNarrativeAgent
            qc_decision: QCDecision from QCReviewAgent
            de_summary: Parsed de_summary.json dict
            enrichment_results: List of parsed enrichment JSON dicts
            qc_metrics: Parsed qc_metrics.json dict
            contrasts: List of contrast config dicts
            species: Species string
            output_path: Where to write the HTML file

        Returns:
            Absolute path to the written HTML file
        """
        # Encode QC plots
        qc_plot_data = {}
        if "plots" in qc_metrics:
            qc_plot_data = self._encode_plots(qc_metrics["plots"])

        # Encode DE plots per contrast
        de_plot_data = {}
        for ct_name, ct_result in de_summary.items():
            if "plots" in ct_result:
                de_plot_data[ct_name] = self._encode_plots(ct_result["plots"])

        # Encode enrichment plots per contrast
        enrich_plot_data = {}
        for enrich in enrichment_results:
            ct_name = enrich.get("contrast_name", "unknown")
            if "plots" in enrich:
                enrich_plot_data[ct_name] = self._encode_plots(enrich["plots"])

        # Build enrichment index by contrast name for easy template lookup
        enrichment_by_contrast = {e["contrast_name"]: e for e in enrichment_results}

        template = self.env.get_template("report.html.jinja2")

        html_content = template.render(
            narrative=narrative,
            qc_decision=qc_decision,
            de_summary=de_summary,
            enrichment_by_contrast=enrichment_by_contrast,
            qc_plot_data=qc_plot_data,
            de_plot_data=de_plot_data,
            enrich_plot_data=enrich_plot_data,
            contrasts=contrasts,
            species=species,
            qc_metrics=qc_metrics,
            edger_sensitivity=edger_sensitivity or {},
        )

        output_path = os.path.abspath(output_path)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
        logger.info(f"HTML report written: {output_path} ({file_size_mb:.1f} MB)")
        return output_path
