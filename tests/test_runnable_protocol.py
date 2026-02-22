import unittest
import asyncio
import inspect
from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.pipeline import DESeq2Pipeline

class TestRunnableProtocol(unittest.TestCase):
    def test_base_agent_has_invoke_methods(self):
        """Verify BaseAgent has invoke/ainvoke methods."""
        self.assertTrue(hasattr(BaseAgent, "invoke"))
        self.assertTrue(hasattr(BaseAgent, "ainvoke"))
        
        # Verify legacy methods are removed
        self.assertFalse(hasattr(BaseAgent, "run"))
        self.assertFalse(hasattr(BaseAgent, "arun"))

    def test_pipeline_has_invoke_methods(self):
        """Verify DESeq2Pipeline has invoke/ainvoke methods."""
        self.assertTrue(hasattr(DESeq2Pipeline, "invoke"))
        self.assertTrue(hasattr(DESeq2Pipeline, "ainvoke"))
        
        # Verify legacy methods are removed
        self.assertFalse(hasattr(DESeq2Pipeline, "run"))
        self.assertFalse(hasattr(DESeq2Pipeline, "arun"))

    def test_agent_ainvoke_signatures(self):
        """Verify async signatures exist (static check via inspection)."""
        self.assertTrue(inspect.iscoroutinefunction(BaseAgent.ainvoke))
        self.assertTrue(inspect.iscoroutinefunction(DESeq2Pipeline.ainvoke))

if __name__ == "__main__":
    unittest.main()
