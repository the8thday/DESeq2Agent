# DESeq2Agent

**End-to-end RNA-seq differential expression pipeline powered by R + LLM agents**

DESeq2Agent takes raw count data and metadata, automatically detects the experimental design, runs a full DESeq2 analysis via R subprocesses, and uses LLM agents to interpret results at each stage — delivering a self-contained HTML report with embedded plots, tables, and narrative.

## Features

- **Intelligent design detection**: `DesignDetectionAgent` auto-detects experimental design (simple two-group, longitudinal, paired, factorial, multi-group) and generates appropriate contrasts — no manual contrast specification required
- **True end-to-end pipeline**: raw counts CSV → HTML report, no pre-processing required
- **Reactome pathway analysis**: GSEA + ORA via ReactomePA for human samples, cross-validated with GO/KEGG
- **edgeR sensitivity analysis**: Independent cross-method validation (edgeR glmQLFTest) with concordance statistics
- **R subprocess integration**: DESeq2, edgeR, clusterProfiler, and QC all run as R scripts; Python orchestrates and validates
- **5 LLM agents** with structured Pydantic v2 outputs:
  - **DesignDetectionAgent**: analyzes metadata, detects design type, generates contrasts
  - **QCReviewAgent**: interprets QC metrics, decides outlier removal
  - **DEReviewAgent**: interprets differential expression results
  - **PathwayReviewAgent**: interprets GSEA/ORA enrichment
  - **ReportNarrativeAgent**: synthesizes a full report narrative
- **Flexible input**: works with or without a config file; contrasts auto-detected when omitted
- **Complex design support**: longitudinal/time-series data via per-timepoint subsetting (`subset_column`/`subset_value`)
- **Interactive mode**: LLM proposes outlier removal; user confirms before DE analysis
- **Multi-provider support**: OpenAI, DeepSeek, and any OpenAI-compatible endpoint
- **Multi-species**: human, mouse, rat, dog (beagle)
- **Publication-quality figures**: all plots saved as PDF with vector graphics
- **Self-contained HTML report**: PDFs converted to PNG and base64-encoded, ~1–5 MB, no external dependencies
- **Full enrichment tables**: GSEA/ORA results exported as CSV in addition to plots

## Requirements

**Python**
- Python ≥ 3.9
- langchain ≥ 0.3.0, langchain-openai ≥ 0.2.0
- pydantic ≥ 2.0, jinja2 ≥ 3.1.0, pandas ≥ 2.0.0, python-dotenv ≥ 1.0

**R** (≥ 4.2, Bioconductor ≥ 3.16)
```r
BiocManager::install(c("DESeq2", "DEGreport", "clusterProfiler",
                        "ReactomePA", "reactome.db",  # Reactome (human only)
                        "edgeR",                       # sensitivity analysis
                        "org.Hs.eg.db", "org.Mm.eg.db",
                        "org.Rn.eg.db",   # rat
                        "org.Cf.eg.db",   # dog (beagle)
                        "AnnotationDbi", "enrichplot"))
install.packages(c("jsonlite", "ggplot2", "pheatmap", "ggrepel",
                   "RColorBrewer", "dplyr", "tibble", "ggsci", "ggdendro"))
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

### Mode 1: Auto-detect everything (no config file)

The Agent analyzes your metadata and automatically generates appropriate contrasts:

```bash
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --output results/ \
  --provider deepseek
```

### Mode 2: Config file with species only (contrasts auto-detected)

```json
// config.json
{"species": "dog"}
```

```bash
python examples/run_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --config config.json \
  --output results/ \
  --provider deepseek
```

### Mode 3: Config file with explicit contrasts

```json
// config.json
{
  "contrasts": [
    {
      "name": "DrugvsControl",
      "variable": "group",
      "treatment": "drug",
      "control": "control"
    }
  ],
  "species": "dog"
}
```

```bash
python examples/run_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --config config.json \
  --output results/ \
  --provider deepseek
```

### Interactive mode

```bash
python examples/run_pipeline.py \
  --counts datademo/counts.csv \
  --metadata datademo/metadata.csv \
  --config datademo/config.json \
  --output results/ \
  --interactive --provider deepseek
```

Output: `results/report.html` plus all intermediate CSVs, JSONs, and PDFs.

### Parameter priority

`CLI flags` > `config.json values` > `hardcoded defaults`

```bash
# CLI --species overrides config species; --padj/--lfc override config thresholds
python examples/run_pipeline.py \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --config config.json \
  --species human \
  --padj 0.01 --lfc 1.5 \
  --output results/ --provider deepseek
```

## Config File Format

Minimal (contrasts auto-detected):
```json
{
  "species": "dog"
}
```

Full (all options):
```json
{
  "contrasts": [
    {
      "name": "DrugvsControl",
      "variable": "group",
      "treatment": "drug",
      "control": "control"
    },
    {
      "name": "Day7vsPreDose",
      "variable": "group",
      "treatment": "Day7",
      "control": "PreDose",
      "subset_column": "timepoint",
      "subset_value": "Day7"
    }
  ],
  "species": "dog",
  "padj_threshold": 0.05,
  "lfc_threshold": 1.0
}
```

Note: `counts_file`, `metadata_file`, and `output_dir` are always specified via CLI, not in config.

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

## Pipeline Overview

```
counts.csv + metadata.csv + optional config.json
  │
  ▼ [LLM] Step 0: DesignDetectionAgent
  │   Detect experimental design, generate contrasts
  │   (skipped if contrasts provided in config)
  │
  ▼ [R] Step 1: 01_data_prep.R
  │   Gene ID detection, ENSEMBL mapping, low-count filtering
  │
  ▼ [R] Step 2: 02_qc_analysis.R
  │   PCA, MDS, heatmap, hclust, IQR-based outlier flagging
  │
  ▼ [LLM] Step 3: QCReviewAgent
  │   Outlier removal decision → QCDecision
  │
  ▼ [R] Step 4: 03_de_analysis.R
  │   DESeq2 + apeglm shrinkage, volcano/MA plots
  │   Supports per-timepoint subsetting
  │
  ▼ [LLM] Step 5: DEReviewAgent → DEReviewOutput
  │
  ▼ [R] Step 6: 04_enrichment.R
  │   GSEA + ORA (GO/KEGG/Reactome) via clusterProfiler + ReactomePA
  │   Reactome only for human; full result tables exported as CSV
  │
  ▼ [R] Step 6b: 05_edger_sensitivity.R
  │   edgeR glmQLFTest cross-validation (non-fatal)
  │   Concordance stats: overlap, Spearman ρ, direction agreement
  │
  ▼ [LLM] Step 7: PathwayReviewAgent → PathwayReviewOutput
  │
  ▼ [LLM] Step 8: ReportNarrativeAgent → ReportNarrative
  │
  ▼ [HTML] Step 9: HTMLReportBuilder
      Self-contained report.html (~1–5 MB)
```

Python writes a `config.json` before each R step; R reads it, writes outputs, exits 0 or non-zero. Python validates expected files after each step.

## Project Structure

```
DESeq2Agent/
├── deseq2_agent/
│   ├── config.py              LLMConfig + get_llm() factory
│   ├── runner.py              RScriptRunner — subprocess + config I/O
│   ├── models.py              Pydantic v2 output models
│   ├── prompts.py             5 ChatPromptTemplates (Simplified Chinese)
│   ├── pipeline.py            DESeq2Pipeline 10-step orchestrator
│   ├── agents/
│   │   ├── base.py            Abstract BaseAgent (LCEL chain)
│   │   ├── design_agent.py    DesignDetectionAgent (Step 0)
│   │   ├── qc_agent.py        QCReviewAgent + ConsoleIO (interactive)
│   │   ├── de_agent.py        DEReviewAgent
│   │   ├── pathway_agent.py   PathwayReviewAgent
│   │   └── report_agent.py    ReportNarrativeAgent
│   └── report/
│       ├── html_builder.py    HTMLReportBuilder (Jinja2 + PDF→PNG + base64)
│       └── templates/
│           └── report.html.jinja2
├── r_scripts/
│   ├── 01_data_prep.R         Gene ID mapping, DDS construction
│   ├── 02_qc_analysis.R       QC plots (PDF) + IQR outlier detection
│   ├── 03_de_analysis.R       DESeq2 + lfcShrink(apeglm), subset support
│   ├── 04_enrichment.R        clusterProfiler + ReactomePA GSEA/ORA, CSV export
│   └── 05_edger_sensitivity.R edgeR glmQLFTest cross-validation
├── datademo/
│   ├── counts.csv             100-gene × 6-sample synthetic demo
│   ├── metadata.csv
│   └── config.json            Demo contrast config
├── examples/
│   └── run_pipeline.py        CLI entry point
└── tests/
    ├── smoke_test.py          Import + Pydantic verification
    └── test_runnable_protocol.py  API contract tests
```

## Python API

```python
from deseq2_agent import DESeq2Pipeline, get_llm

llm = get_llm()  # reads from .env
pipeline = DESeq2Pipeline(llm, mode="auto")  # or mode="interactive"

# Auto-detect contrasts from metadata
results = pipeline.run(
    counts_file="datademo/counts.csv",
    metadata_file="datademo/metadata.csv",
    contrasts=None,   # DesignDetectionAgent auto-generates
    species="human",
    output_dir="results/"
)

# Or provide explicit contrasts
results = pipeline.run(
    counts_file="datademo/counts.csv",
    metadata_file="datademo/metadata.csv",
    contrasts=[{
        "name": "TreatvsCtrl",
        "variable": "condition",
        "treatment": "Treatment",
        "control": "Control",
    }],
    species="human",
    output_dir="results/"
)

print(results.report_path)           # path to report.html
print(results.design_decision)       # DesignDecision from Step 0
print(results.qc_decision.samples_to_remove)
print(results.de_review.contrasts[0].overall_assessment)
print(results.narrative.executive_summary)
```

## Supported Experimental Designs

| Design Type | Example | Analysis Strategy |
|---|---|---|
| `simple_two_group` | Treatment vs Control | Single DESeq2 run |
| `multi_group` | 3+ treatment groups | Multiple pairwise contrasts |
| `longitudinal` | Time-series (Day0, Day3, Day7, ...) | Per-timepoint subsetting |
| `paired` | Before/after on same subjects | Paired design formula |
| `factorial` | Drug × Dose combinations | Interaction or per-condition |

For longitudinal designs, the agent generates contrasts with `subset_column`/`subset_value` to run DESeq2 on each timepoint independently.

## CLI Options

```
python examples/run_pipeline.py --help

  --counts      Path to raw counts CSV (required)
  --metadata    Path to metadata CSV (required)
  --config      Path to analysis config JSON (optional)
  --output      Output directory (default: results/)
  --species     human | mouse | rat | dog (default: human, overrides config)
  --provider    openai | deepseek (default: openai)
  --interactive Prompt user to confirm outlier removal
  --padj        Adjusted p-value threshold (default: 0.05, overrides config)
  --lfc         Log2 fold-change threshold (default: 1.0, overrides config)
```

## Output Files

```
results/
├── report.html                       Self-contained HTML report
├── dds.rds                           DESeq2 object
├── counts_ensembl.csv                ENSEMBL-mapped, filtered counts
├── id_mapping.csv                    Gene ID cross-reference table
├── config.json                       Pipeline config (written by Python)
├── data_prep_summary.json
├── qc_metrics.json                   QC metrics + outlier flags
├── qc_plots/                         PCA, MDS, heatmap, hclust, ... (PDF)
├── de_summary.json
├── {contrast}_results.csv            Full DE results
├── {contrast}_sig_up.csv
├── {contrast}_sig_down.csv
├── de_plots/                         Volcano + MA plots per contrast (PDF)
├── enrichment_results_{contrast}.json
├── enrichment/
│   ├── plots/{contrast}/             GSEA/ORA dotplots (PDF)
│   └── tables/{contrast}/            Full enrichment CSV tables
│       ├── gsea_go.csv
│       ├── gsea_kegg.csv
│       ├── ora_go.csv
│       ├── ora_kegg.csv
│       ├── gsea_reactome.csv         (human only)
│       └── ora_reactome.csv          (human only)
├── edger_{contrast}_results.csv      edgeR DE results (sensitivity)
├── edger_sensitivity.json            Cross-method concordance stats
```

## License

MIT License
