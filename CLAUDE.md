# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DESeq2Agent is an end-to-end RNA-seq differential expression analysis pipeline (v2.0.0). It runs DESeq2, QC, and enrichment analysis via R subprocesses, with LLM agents interpreting results at key stages and producing a self-contained HTML report.

Input: raw counts CSV + metadata CSV + optional config JSON (contrasts auto-detected if omitted)
Output: `results/report.html` + all intermediate CSVs, JSONs, and PDFs

## Commands

```bash
# Development install (editable)
pip install -e ".[dev]"

# Minimal: auto-detect design and contrasts from metadata
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/ \
  --provider deepseek

# With config (contrasts defined, skips auto-detection)
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --config datademo/config.json \
  --output results/ \
  --provider deepseek

# With config (species only, contrasts auto-detected)
# config.json: {"species": "dog"}
python examples/run_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --config path/to/config.json \
  --output ./results/ \
  --provider deepseek

# Full example with all optional CLI overrides (CLI > config > defaults)
python examples/run_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --config path/to/config.json \
  --output ./results/ \
  --species human \       # human | mouse | rat | dog (overrides config)
  --provider deepseek \   # openai | deepseek
  --padj 0.05 \           # adjusted p-value threshold (overrides config default: 0.05)
  --lfc 1.0               # log2 fold-change threshold (overrides config default: 1.0)

# Interactive mode (LLM proposes outlier removal, user confirms)
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --config datademo/config.json \
  --interactive --provider deepseek

# Run tests
pytest tests/

# Run a single test file
pytest tests/smoke_test.py -v

# Lint and format
ruff check deseq2_agent/
black deseq2_agent/
```

## Configuration

Copy `.env.example` to `.env`:
- `OPENAI_API_KEY` + `OPENAI_MODEL` (default: `gpt-4o`)
- Or DeepSeek: `DEEPSEEK_API_KEY` + `DEEPSEEK_BASE_URL=https://api.deepseek.com/v1`
- Optional: `LLM_TEMPERATURE` (default: 0.2), `LLM_MAX_RETRIES` (default: 3)

Any OpenAI-compatible endpoint works via `base_url` in `LLMConfig`.

## Architecture

### 10-Step Pipeline (Step 0–9)

```
counts.csv + metadata.csv + optional config.json
  │
  ▼ [LLM] Step 0: DesignDetectionAgent
  │   Analyzes metadata structure (columns, unique values, sample counts)
  │   Detects design type: simple_two_group | longitudinal | paired | factorial | multi_group
  │   Suggests analysis_strategy: single_run | per_timepoint | per_condition
  │   Auto-generates contrasts (with subset_column/subset_value for complex designs)
  │   → DesignDecision; skipped if user provides contrasts in config
  │
  ▼ [R] Step 1: 01_data_prep.R
  │   Auto-detect gene IDs (ENSEMBL/ENTREZ/SYMBOL), normalize to ENSEMBL,
  │   filter low counts, write dds.rds + id_mapping.csv + counts_ensembl.csv
  │
  ▼ [R] Step 2: 02_qc_analysis.R
  │   Auto-detect informative metadata columns (2 ≤ n_unique ≤ 10 AND n_unique < n_samples)
  │   PCA (top 2000 variable genes) + MDS — one plot per qualifying column → pca_{col}.pdf / mds_{col}.pdf
  │   Correlation heatmap, hierarchical clustering (hclust.pdf), dispersion, library sizes
  │   IQR-based outlier flagging → qc_metrics.json + qc_plots/*.pdf
  │
  ▼ [LLM] Step 3: QCReviewAgent
  │   Decides outlier removal and data usability → QCDecision
  │   interactive mode: prompts user per candidate; auto mode: LLM decides
  │
  ▼ [R] Step 4: 03_de_analysis.R
  │   DESeq2 + lfcShrink(apeglm) excluding removed samples;
  │   Supports subset_column/subset_value for per-timepoint contrasts
  │   volcano + MA plots → {contrast}_results.csv + de_summary.json
  │
  ▼ [LLM] Step 5: DEReviewAgent → DEReviewOutput
  │
  ▼ [R] Step 6: 04_enrichment.R
  │   GSEA (GO+KEGG) + ORA (GO+KEGG) via clusterProfiler; set.seed(42)
  │   Skips ORA if <10 sig genes
  │   → enrichment_results_{contrast}.json + dotplots (PDF)
  │   → enrichment/tables/{contrast}/*.csv (full result tables)
  │
  ▼ [LLM] Step 7: PathwayReviewAgent → PathwayReviewOutput
  │
  ▼ [LLM] Step 8: ReportNarrativeAgent → ReportNarrative
  │
  ▼ [HTML] Step 9: HTMLReportBuilder
      Converts PDFs → PNGs (via sips on macOS), base64-encodes into
      single self-contained report.html (~1-5 MB)
```

### R-Python Communication

Before each R step: Python writes `{output_dir}/config.json` (atomic write via temp file + `os.replace`). R reads it, computes, writes outputs, exits 0 or non-zero. Python validates expected output files exist after each step.

**Critical**: All R scripts use `fromJSON(args[1], simplifyDataFrame = FALSE)` to keep JSON arrays-of-objects as lists-of-lists (not data frames).

### Key Modules

- `deseq2_agent/config.py`: `LLMConfig` dataclass + `get_llm()` factory
- `deseq2_agent/runner.py`: `RScriptRunner` — `write_config()`, `run()`, `check_r_available()`, `validate_outputs()`; `RScriptError` exception; `R_SCRIPTS_DIR = Path(__file__).parent.parent / "r_scripts"`
- `deseq2_agent/models.py`: Pydantic v2 output models — `DesignDecision`, `ContrastSuggestion`, `QCDecision`, `DEReviewOutput`, `PathwayReviewOutput`, `ReportNarrative` (and their nested types)
- `deseq2_agent/prompts.py`: 5 `ChatPromptTemplate` instances — **all in Simplified Chinese**
- `deseq2_agent/agents/base.py`: Abstract `BaseAgent`; builds LCEL chain as `prompt_template | llm.with_structured_output(output_model, method="function_calling")`
- `deseq2_agent/agents/design_agent.py`: `DesignDetectionAgent` — Step 0, auto-detects experimental design
- `deseq2_agent/agents/qc_agent.py`: `QCReviewAgent` + `ConsoleIO` for interactive mode
- `deseq2_agent/pipeline.py`: `DESeq2Pipeline.run()` — 10-step orchestrator; `PipelineConfig`, `PipelineResults`, `_build_metadata_summary()`
- `deseq2_agent/report/html_builder.py`: `HTMLReportBuilder.build()` — Jinja2 templating + PDF→PNG conversion + base64 embedding
- `deseq2_agent/report/templates/report.html.jinja2`: Self-contained HTML template

### QC Plots Naming Convention

`02_qc_analysis.R` auto-detects informative grouping columns from metadata using the heuristic:
`2 ≤ n_unique ≤ 10 AND n_unique < n_samples` — skips per-sample IDs and constant columns.

- PCA/MDS plots: `qc_plots/pca_{col}.pdf`, `qc_plots/mds_{col}.pdf` (one per qualifying column)
- Fixed plots: `qc_plots/heatmap.pdf`, `qc_plots/hclust.pdf`, `qc_plots/dispersion.pdf`, `qc_plots/library_sizes.pdf`, `qc_plots/detected_genes.pdf`
- `primary_col` = first qualifying column, used for heatmap/hclust/library_sizes/detected_genes coloring
- Falls back to first metadata column if no column qualifies

`qc_metrics.json → plots` dict is built dynamically with keys `pca_{col}` and `mds_{col}`.
`html_builder.py._encode_plots()` encodes all keys dynamically — no hardcoded plot names.
The Jinja2 template renders PCA/MDS sections by iterating keys with `pca_` / `mds_` prefix.

### BaseAgent Pattern

Each agent subclass declares two class attributes:
```python
class QCReviewAgent(BaseAgent):
    prompt_template = QC_REVIEW_PROMPT   # from prompts.py
    output_model = QCDecision            # from models.py
```

`validate_input()` must check required fields and serialize `dict`/`list` values to JSON strings via `json.dumps(..., ensure_ascii=False)`.

**Cross-provider compatibility**: `method="function_calling"` is required in `with_structured_output()` — DeepSeek does not support the default `json_schema` response format.

### DesignDetectionAgent (Step 0)

Runs before any R scripts. Analyzes metadata structure via `_build_metadata_summary()` which reads the raw CSV and produces per-column stats (n_unique, unique_values, value_counts).

**Contrast auto-detection flow**:
1. If user provides contrasts in config → use them directly, DesignDetectionAgent runs as advisory only
2. If no contrasts → DesignDetectionAgent generates `suggested_contrasts` → pipeline uses them
3. In interactive mode → user sees design decision and confirms before proceeding

**Supported design types**: `simple_two_group`, `longitudinal`, `paired`, `factorial`, `multi_group`

**Per-timepoint subsetting**: For longitudinal designs, contrasts include `subset_column` and `subset_value` fields. `03_de_analysis.R` filters samples to only those matching the subset before running DESeq2.

**Config priority**: CLI flags > config.json values > hardcoded defaults

### R Packages Required

```r
BiocManager::install(c("DESeq2", "DEGreport", "clusterProfiler",
                        "org.Hs.eg.db", "org.Mm.eg.db",
                        "org.Rn.eg.db",   # rat
                        "org.Cf.eg.db",   # dog (beagle)
                        "AnnotationDbi", "enrichplot"))
install.packages(c("jsonlite", "ggplot2", "pheatmap", "ggrepel",
                   "RColorBrewer", "dplyr", "tibble", "ggsci", "ggdendro"))
```

Requires: R ≥ 4.2, Bioconductor ≥ 3.16

### Known Serialization Gotchas

| Issue | Fix |
|-------|-----|
| R `c()` on list-of-lists → JSON object with integer keys | `pipeline.py` normalizes with `list(outlier_flags_raw.values())` if `isinstance(..., dict)` |
| `simplifyDataFrame = TRUE` (R default) collapses contrast arrays | All R scripts use `simplifyDataFrame = FALSE` |
| `samples_to_remove` from empty JSON list → R list not character vector | `as.character(unlist(config$samples_to_remove %||% list()))` |
| `ensembl_to_entrez()` must return same length as input (for building data frames) | NAs preserved; `entrez_vec()` is the filtered version for ORA/GSEA |
| GSEA results non-deterministic across runs | `set.seed(42)` in `04_enrichment.R` before any analysis |
| DE summary JSON truncated for many contrasts | `de_summary_text` truncated at 16000 chars (increased from 8000) |

## Prompts Language

All 5 prompts (`DESIGN_DETECTION_PROMPT`, `QC_REVIEW_PROMPT`, `DE_REVIEW_PROMPT`, `PATHWAY_REVIEW_PROMPT`, `REPORT_NARRATIVE_PROMPT`) are in **Simplified Chinese** with domain-specific bioinformatics roles. Do not modify without understanding Chinese.

## Design Constraints

Per the prompts, agents must not:
- Fabricate gene/pathway results not present in the provided data
- Overstate causal relationships or use language stronger than evidence supports
- Modify or reanalyze the original data
