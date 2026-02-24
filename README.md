# DESeq2Agent

**End-to-end RNA-seq differential expression pipeline powered by R + LLM agents**

DESeq2Agent takes raw count data and metadata, runs a full DESeq2 analysis via R subprocesses, and uses LLM agents to interpret results at each stage — delivering a self-contained HTML report with embedded plots, tables, and narrative.

## Features

- **True end-to-end pipeline**: raw counts CSV → HTML report, no pre-processing required
- **R subprocess integration**: DESeq2, clusterProfiler, and QC all run as R scripts; Python orchestrates and validates
- **4 LLM agents** with structured Pydantic outputs:
  - **QCReviewAgent**: interprets QC metrics, decides outlier removal
  - **DEReviewAgent**: interprets differential expression results
  - **PathwayReviewAgent**: interprets GSEA/ORA enrichment
  - **ReportNarrativeAgent**: synthesizes a full report narrative
- **Interactive mode**: LLM proposes outlier removal; user confirms before DE analysis
- **Multi-provider support**: OpenAI, DeepSeek, and any OpenAI-compatible endpoint
- **Self-contained HTML report**: all plots base64-encoded, ~1–5 MB, no external dependencies

## Requirements

**Python**
- Python ≥ 3.9
- langchain ≥ 0.3.0, langchain-openai ≥ 0.2.0
- pydantic ≥ 2.0, jinja2 ≥ 3.1.0, pandas ≥ 2.0.0, python-dotenv ≥ 1.0

**R** (≥ 4.2, Bioconductor ≥ 3.16)
```r
BiocManager::install(c("DESeq2", "DEGreport", "clusterProfiler",
                        "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "enrichplot"))
install.packages(c("jsonlite", "ggplot2", "pheatmap", "ggrepel",
                   "RColorBrewer", "dplyr", "tibble"))
```

## Installation

```bash
git clone https://github.com/your-org/DESeq2Agent
cd DESeq2Agent

# Install Python package (editable)
pip install -e ".[dev]"
```

## Configuration

```bash
cp .env.example .env
```

Edit `.env`:

```env
# OpenAI
OPENAI_API_KEY=sk-your-key
OPENAI_MODEL=gpt-4o

# DeepSeek (alternative)
# DEEPSEEK_API_KEY=your-key
# DEEPSEEK_BASE_URL=https://api.deepseek.com/v1
# DEEPSEEK_MODEL=deepseek-chat

# Optional
# LLM_TEMPERATURE=0.2
# LLM_MAX_RETRIES=3
```

## Quick Start

```bash
# Run with demo data (human, OpenAI)
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/

# Use DeepSeek
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/ \
  --provider deepseek

# Interactive mode — confirm outlier removal before DE
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/ \
  --interactive
```

Output: `results/report.html` plus all intermediate CSVs, JSONs, and PNGs.

## Input Format

**`counts.csv`** — raw integer counts, genes as rows, samples as columns:
```
gene_id,ctrl_1,ctrl_2,ctrl_3,treat_1,treat_2,treat_3
ENSG00000141510,120,134,118,89,92,95
ENSG00000012048,45,52,48,210,198,205
...
```
Gene IDs can be ENSEMBL, Entrez, or gene symbols — auto-detected and mapped.

**`metadata.csv`** — one row per sample:
```
sample,condition,batch
ctrl_1,Control,A
ctrl_2,Control,A
treat_1,Treatment,A
treat_2,Treatment,B
...
```

Contrasts are auto-generated from the `condition` column (all non-reference levels vs. the reference). Custom contrasts can be passed via the Python API.

## Pipeline Overview

```
counts.csv + metadata.csv
  │
  ▼ [R] 01_data_prep.R       Gene ID detection, ENSEMBL mapping, low-count filtering
  ▼ [R] 02_qc_analysis.R     PCA, heatmap, MDS, IQR-based outlier flagging
  ▼ [LLM] QCReviewAgent      Outlier removal decision → QCDecision
  ▼ [R] 03_de_analysis.R     DESeq2 + apeglm shrinkage, volcano/MA plots
  ▼ [LLM] DEReviewAgent      DE interpretation → DEReviewOutput
  ▼ [R] 04_enrichment.R      GSEA + ORA (GO/KEGG) via clusterProfiler
  ▼ [LLM] PathwayReviewAgent Pathway themes → PathwayReviewOutput
  ▼ [LLM] ReportNarrativeAgent Full narrative → ReportNarrative
  ▼ [HTML] HTMLReportBuilder  Self-contained report.html (~1–5 MB)
```

Python writes a `config.json` before each R step; R reads it, writes outputs, exits 0 or non-zero. Python validates expected files after each step.

## Project Structure

```
DESeq2Agent/
├── deseq2_agent/
│   ├── config.py              LLMConfig + get_llm() factory
│   ├── runner.py              RScriptRunner — subprocess + config I/O
│   ├── models.py              Pydantic v2 output models
│   ├── prompts.py             ChatPromptTemplates (Simplified Chinese)
│   ├── pipeline.py            DESeq2Pipeline 9-step orchestrator
│   ├── agents/
│   │   ├── base.py            Abstract BaseAgent (LCEL chain)
│   │   ├── qc_agent.py        QCReviewAgent + ConsoleIO (interactive)
│   │   ├── de_agent.py        DEReviewAgent
│   │   ├── pathway_agent.py   PathwayReviewAgent
│   │   └── report_agent.py    ReportNarrativeAgent
│   └── report/
│       ├── html_builder.py    HTMLReportBuilder (Jinja2 + base64 PNGs)
│       └── templates/
│           └── report.html.jinja2
├── r_scripts/
│   ├── 01_data_prep.R         Gene ID mapping, DDS construction
│   ├── 02_qc_analysis.R       QC plots + IQR outlier detection
│   ├── 03_de_analysis.R       DESeq2 + lfcShrink(apeglm)
│   └── 04_enrichment.R        clusterProfiler GSEA + ORA
├── datademo/
│   ├── counts.csv             100-gene × 6-sample synthetic demo
│   └── metadata.csv
├── examples/
│   └── run_pipeline.py        CLI entry point
└── tests/
    └── smoke_test.py
```

## Python API

```python
from deseq2_agent import DESeq2Pipeline, get_llm

llm = get_llm()  # reads from .env
pipeline = DESeq2Pipeline(llm, mode="auto")  # or mode="interactive"

results = pipeline.run(
    counts_file="datademo/counts.csv",
    metadata_file="datademo/metadata.csv",
    contrasts=None,   # auto-detect from metadata
    species="human",  # or "mouse"
    output_dir="results/"
)

print(results.report_path)           # path to report.html
print(results.qc_decision.samples_to_remove)
print(results.de_review.contrasts[0].overall_assessment)
print(results.narrative.executive_summary)
```

## CLI Options

```
python examples/run_pipeline.py --help

  --counts      Path to raw counts CSV (required)
  --metadata    Path to metadata CSV (required)
  --output      Output directory (default: results/)
  --species     human or mouse (default: human)
  --provider    openai or deepseek (default: openai)
  --interactive Prompt user to confirm outlier removal
  --padj        Adjusted p-value threshold (default: 0.05)
  --lfc         Log2 fold-change threshold (default: 1.0)
```

## Output Files

```
results/
├── report.html                       Self-contained HTML report
├── dds.rds                           DESeq2 object
├── counts_ensembl.csv                ENSEMBL-mapped, filtered counts
├── id_mapping.csv                    Gene ID cross-reference table
├── data_prep_summary.json
├── qc_metrics.json                   QC metrics + outlier flags
├── qc_plots/                         PCA, heatmap, MDS, dispersion, ...
├── de_summary.json
├── {contrast}_results.csv            Full DE results
├── {contrast}_sig_up.csv
├── {contrast}_sig_down.csv
├── de_plots/                         Volcano + MA plots per contrast
├── enrichment_results_{contrast}.json
└── enrichment/plots/                 GSEA/ORA dotplots
```

## License

MIT License
