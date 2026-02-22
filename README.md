# DESeq2Agent

🧬 **LangChain-based AI Agent Framework for RNA-seq DESeq2 Analysis Interpretation**

DESeq2Agent is a multi-agent system that uses LLMs to interpret RNA-seq differential expression analysis results. It leverages LangChain's Expression Language (LCEL) for composable, type-safe chains.

## Features

- 🔗 **LangChain Integration**: Built on LangChain v0.3+ with LCEL for modern, composable chains
- 🎯 **Structured Outputs**: Pydantic models for type-safe, validated outputs
- 🤖 **5 Specialized Agents**:
  - **DataAgent**: Experimental design understanding and metadata validation
  - **QCAgent**: Quality control assessment and outlier detection
  - **DEAgent**: Differential expression interpretation
  - **PathwayAgent**: GSEA/ORA pathway analysis interpretation
  - **ReportAgent**: Structured report generation
- 🌐 **Multi-Provider Support**: Works with OpenAI, DeepSeek, and other OpenAI-compatible APIs
- 📊 **Pipeline Orchestration**: Complete workflow with automatic result passing between agents

## Installation

```bash
# Clone the repository
cd DESeq2Agent

# Install dependencies
pip install -r requirements.txt

# Or install as a package
pip install -e .
```

## Configuration

1. Copy the environment template:
```bash
cp .env.example .env
```

2. Edit `.env` with your API key:
```env
# For OpenAI
OPENAI_API_KEY=sk-your-api-key-here
OPENAI_MODEL=gpt-4o

# For DeepSeek (alternative)
# DEEPSEEK_API_KEY=your-deepseek-key
# DEEPSEEK_BASE_URL=https://api.deepseek.com/v1
# DEEPSEEK_MODEL=deepseek-chat
```

## Quick Start

### Using Individual Agents

```python
from deseq2_agent import get_llm, DEAgent

# Initialize LLM
llm = get_llm()  # Uses OPENAI_API_KEY from environment

# Create agent
de_agent = DEAgent(llm)

# Run analysis
result = de_agent.run(
    comparison="Treatment vs Control",
    biological_hypothesis="Treatment reduces inflammation",
    cutoff="padj < 0.05 & |log2FC| > 1",
    up_genes=[{"gene": "COL1A1", "log2FC": 2.5, "padj": 1e-10}],
    down_genes=[{"gene": "IL6", "log2FC": -1.5, "padj": 1e-5}]
)

print(result.overall_assessment)
```

### Using the Complete Pipeline

```python
from deseq2_agent import DESeq2Pipeline, get_llm

llm = get_llm()
pipeline = DESeq2Pipeline(llm)

results = pipeline.run(
    sample_metadata={"samples": [...], "design": "paired"},
    comparison_groups={"treatment": "T12", "control": "T0"},
    qc_metrics={"library_sizes": [...], "genes_detected": [...]},
    up_genes=[...],
    down_genes=[...],
    gsea_up=[...],
    gsea_down=[...],
    study_context="Phase II efficacy study"
)

# Access structured outputs
print(results.report.executive_summary)
print(results.report.key_findings)
```

### Using DeepSeek

```python
from deseq2_agent import get_llm, DESeq2Pipeline

# Use DeepSeek instead of OpenAI
llm = get_llm(provider="deepseek")
pipeline = DESeq2Pipeline(llm)
```

## Project Structure

```
DESeq2Agent/
├── deseq2_agent/
│   ├── __init__.py          # Package exports
│   ├── config.py             # LLM configuration
│   ├── models.py             # Pydantic output models
│   ├── prompts.py            # Agent prompt templates
│   ├── pipeline.py           # Pipeline orchestration
│   └── agents/
│       ├── base.py           # Base agent class
│       ├── data_agent.py     # Experimental design
│       ├── qc_agent.py       # Quality control
│       ├── de_agent.py       # Differential expression
│       ├── pathway_agent.py  # Pathway analysis
│       └── report_agent.py   # Report generation
└── examples/
    ├── basic_usage.py        # Individual agent examples
    └── full_pipeline.py      # Complete pipeline example
```

## Agent Outputs

Each agent produces a structured Pydantic model:

### DataAgentOutput
- `experimental_design_summary`: Experiment design overview
- `comparison_validity`: Comparison logic assessment
- `metadata_issues`: List of metadata problems found
- `confounding_factors`: Potential confounders identified
- `recommendations`: Suggestions for improvement

### QCAgentOutput
- `quality_summary`: Overall quality assessment
- `outlier_samples`: List of outliers with reasons
- `batch_effect_assessment`: Batch effect evaluation
- `data_usability`: Whether data supports analysis
- `qc_recommendations`: QC suggestions

### DEAgentOutput
- `overall_assessment`: Hypothesis support evaluation
- `supporting_evidence`: Evidence supporting hypothesis
- `potential_concerns`: Concerns or deviations
- `interpretation_boundary`: Interpretation limits

### PathwayAgentOutput
- `mechanistic_summary`: Core biological mechanisms
- `biological_context`: Context interpretation
- `key_pathways_up/down`: Key regulated pathways
- `missing_or_unexpected_signals`: Unexpected findings

### ReportAgentOutput
- `executive_summary`: 1-2 paragraph summary
- `key_findings`: Prioritized findings list
- `results_narrative`: Report-ready text
- `limitations`: Analysis limitations
- `conclusions`: Final conclusions

## Requirements

- Python >= 3.9
- langchain >= 0.3.0
- langchain-openai >= 0.2.0
- pydantic >= 2.0
- python-dotenv >= 1.0

## License

MIT License

## Related Projects

- [AutoBioInfoPipe](../AutoBioInfoPipe) - Original implementation without LangChain
