import sys
import os

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    print("Testing imports...")
    from deseq2_agent import (
        DESeq2Pipeline,
        DesignDetectionAgent,
        QCReviewAgent,
        DEReviewAgent,
        PathwayReviewAgent,
        ReportNarrativeAgent,
    )
    from deseq2_agent.config import LLMConfig
    print("✅ Imports successful")

    print("\nTesting DeepSeek dependency...")
    import langchain_deepseek
    print("✅ langchain-deepseek available")

    print("\nTesting Pydantic models...")
    from deseq2_agent.models import DesignDecision, QCDecision, DEReviewOutput
    # Verify Pydantic v2 usage
    assert hasattr(DesignDecision, "model_dump"), "Pydantic v2 .model_dump() not found"
    assert hasattr(QCDecision, "model_dump"), "Pydantic v2 .model_dump() not found"
    assert hasattr(DEReviewOutput, "model_dump"), "Pydantic v2 .model_dump() not found"
    print("✅ Pydantic v2 confirmed")

    print("\nTesting version requirements...")
    from importlib import metadata
    print(f"LangChain version: {metadata.version('langchain')}")
    print(f"LangChain Core version: {metadata.version('langchain-core')}")
    print(f"LangChain OpenAI version: {metadata.version('langchain-openai')}")

    print("\n✅ Verification complete!")

except Exception as e:
    print(f"\n❌ Verification Failed: {e}")
    sys.exit(1)
