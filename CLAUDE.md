# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DESeq2Agent is an end-to-end RNA-seq differential expression analysis pipeline. It runs DESeq2, QC, and enrichment analysis via R subprocesses, with LLM agents interpreting results at key stages and producing a self-contained HTML report.

Input: raw counts CSV + metadata CSV + contrast config
Output: `results/report.html` + all intermediate CSVs, JSONs, and PNGs

## Commands

```bash
# Development install (editable)
pip install -e ".[dev]"

# Run the pipeline with demo data
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/ \
  --provider deepseek   # or openai

# Interactive mode (LLM proposes outlier removal, user confirms)
python examples/run_pipeline.py --counts ... --interactive

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

### 9-Step Pipeline

```
counts.csv + metadata.csv + contrasts
  │
  ▼ [R] 01_data_prep.R
  │   Auto-detect gene IDs (ENSEMBL/ENTREZ/SYMBOL), normalize to ENSEMBL,
  │   filter low counts, write dds.rds + id_mapping.csv + counts_ensembl.csv
  │
  ▼ [R] 02_qc_analysis.R
  │   PCA, correlation heatmap, MDS, dispersion, library sizes;
  │   IQR-based outlier flagging → qc_metrics.json + qc_plots/*.png
  │
  ▼ [LLM] QCReviewAgent
  │   Decides outlier removal and data usability → QCDecision
  │   interactive mode: prompts user per candidate; auto mode: LLM decides
  │
  ▼ [R] 03_de_analysis.R
  │   DESeq2 + lfcShrink(apeglm) excluding removed samples;
  │   volcano + MA plots → {contrast}_results.csv + de_summary.json
  │
  ▼ [LLM] DEReviewAgent → DEReviewOutput
  │
  ▼ [R] 04_enrichment.R
  │   GSEA (GO+KEGG) + ORA (GO+KEGG) via clusterProfiler;
  │   skips ORA if <10 sig genes → enrichment_results_{contrast}.json + dotplots
  │
  ▼ [LLM] PathwayReviewAgent → PathwayReviewOutput
  │
  ▼ [LLM] ReportNarrativeAgent → ReportNarrative
  │
  ▼ [HTML] HTMLReportBuilder
      Base64-encodes all PNGs into single self-contained report.html (~1-5 MB)
```

### R-Python Communication

Before each R step: Python writes `{output_dir}/config.json` (atomic write via temp file + `os.replace`). R reads it, computes, writes outputs, exits 0 or non-zero. Python validates expected output files exist after each step.

**Critical**: All R scripts use `fromJSON(args[1], simplifyDataFrame = FALSE)` to keep JSON arrays-of-objects as lists-of-lists (not data frames).

### Key Modules

- `deseq2_agent/config.py`: `LLMConfig` dataclass + `get_llm()` factory
- `deseq2_agent/runner.py`: `RScriptRunner` — `write_config()`, `run()`, `check_r_available()`, `validate_outputs()`; `RScriptError` exception; `R_SCRIPTS_DIR = Path(__file__).parent.parent / "r_scripts"`
- `deseq2_agent/models.py`: Pydantic v2 output models — `QCDecision`, `DEReviewOutput`, `PathwayReviewOutput`, `ReportNarrative` (and their nested types)
- `deseq2_agent/prompts.py`: 4 `ChatPromptTemplate` instances — **all in Simplified Chinese**
- `deseq2_agent/agents/base.py`: Abstract `BaseAgent`; builds LCEL chain as `prompt_template | llm.with_structured_output(output_model, method="function_calling")`
- `deseq2_agent/agents/qc_agent.py`: `QCReviewAgent` + `ConsoleIO` for interactive mode
- `deseq2_agent/pipeline.py`: `DESeq2Pipeline.run()` — 9-step orchestrator; `PipelineConfig`, `PipelineResults`
- `deseq2_agent/report/html_builder.py`: `HTMLReportBuilder.build()` — Jinja2 templating + base64 PNG embedding
- `deseq2_agent/report/templates/report.html.jinja2`: Self-contained HTML template

### BaseAgent Pattern

Each agent subclass declares two class attributes:
```python
class QCReviewAgent(BaseAgent):
    prompt_template = QC_REVIEW_PROMPT   # from prompts.py
    output_model = QCDecision            # from models.py
```

`validate_input()` must check required fields and serialize `dict`/`list` values to JSON strings via `json.dumps(..., ensure_ascii=False)`.

**Cross-provider compatibility**: `method="function_calling"` is required in `with_structured_output()` — DeepSeek does not support the default `json_schema` response format.

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

## Prompts Language

All prompts are in **Simplified Chinese** with domain-specific bioinformatics roles. Do not modify without understanding Chinese.

## Design Constraints

Per the prompts, agents must not:
- Fabricate gene/pathway results not present in the provided data
- Overstate causal relationships or use language stronger than evidence supports
- Modify or reanalyze the original data
