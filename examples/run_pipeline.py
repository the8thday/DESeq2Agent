"""Example: Run the full DESeq2Agent pipeline.

Usage:
    # Auto mode (LLM decides outlier removal):
    python examples/run_pipeline.py

    # Interactive mode (user confirms each outlier):
    python examples/run_pipeline.py --interactive

    # Custom inputs:
    python examples/run_pipeline.py \
        --counts path/to/counts.csv \
        --metadata path/to/metadata.csv \
        --output ./my_results/

Requirements:
    - .env file with OPENAI_API_KEY (or DEEPSEEK_API_KEY)
    - R ≥ 4.2 with required Bioconductor packages
    - pip install -e ".[dev]"
"""

import argparse
import logging
import os
import sys
from pathlib import Path

# Add project root to path when running as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv
load_dotenv()

from deseq2_agent import DESeq2Pipeline, get_llm
from deseq2_agent.config import LLMConfig

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="DESeq2Agent — end-to-end RNA-seq analysis pipeline"
    )
    parser.add_argument(
        "--counts",
        default="datademo/counts.csv",
        help="Path to counts CSV (rows=genes, cols=samples)",
    )
    parser.add_argument(
        "--metadata",
        default="datademo/metadata.csv",
        help="Path to metadata CSV (rows=samples)",
    )
    parser.add_argument(
        "--output",
        default="./results/",
        help="Output directory",
    )
    parser.add_argument(
        "--species",
        default="human",
        choices=["human", "mouse"],
        help="Species for annotation",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Use interactive mode: prompt user to confirm outlier removal",
    )
    parser.add_argument(
        "--provider",
        default="openai",
        choices=["openai", "deepseek"],
        help="LLM provider",
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold (default: 0.05)",
    )
    parser.add_argument(
        "--lfc",
        type=float,
        default=1.0,
        help="log2 fold change threshold (default: 1.0)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Validate inputs
    counts_path = Path(args.counts)
    metadata_path = Path(args.metadata)
    if not counts_path.exists():
        logger.error(f"Counts file not found: {counts_path}")
        sys.exit(1)
    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        sys.exit(1)

    # Load LLM
    logger.info(f"Loading LLM (provider={args.provider})...")
    try:
        llm = get_llm(provider=args.provider)
    except ValueError as e:
        logger.error(f"LLM initialization failed: {e}")
        logger.error("Please set OPENAI_API_KEY (or DEEPSEEK_API_KEY) in your .env file")
        sys.exit(1)

    # Define contrasts
    # This matches the demo metadata.csv which has a 'condition' column
    # with values 'Control' and 'Treatment'
    contrasts = [
        {
            "name": "TreatmentvsControl",
            "variable": "condition",
            "treatment": "Treatment",
            "control": "Control",
        }
    ]

    mode = "interactive" if args.interactive else "auto"
    logger.info(f"Starting pipeline: mode={mode}, species={args.species}")
    logger.info(f"Counts: {counts_path}")
    logger.info(f"Metadata: {metadata_path}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Contrasts: {[c['name'] for c in contrasts]}")

    pipeline = DESeq2Pipeline(llm=llm, mode=mode)

    try:
        results = pipeline.run(
            counts_file=str(counts_path),
            metadata_file=str(metadata_path),
            contrasts=contrasts,
            species=args.species,
            output_dir=args.output,
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    if results.qc_decision:
        print(f"QC: {results.qc_decision.data_usability} | "
              f"Batch effect: {results.qc_decision.batch_effect_severity}")
        if results.qc_decision.samples_to_remove:
            print(f"     Removed samples: {results.qc_decision.samples_to_remove}")
    if results.de_review:
        for ct in results.de_review.contrasts:
            print(f"DE [{ct.contrast_name}]: coherence={ct.biological_coherence}")
    if results.report_path:
        print(f"\nReport: {results.report_path}")
        size_mb = os.path.getsize(results.report_path) / (1024 * 1024)
        print(f"  Size: {size_mb:.1f} MB")
        print(f"  Open in browser: file://{os.path.abspath(results.report_path)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
