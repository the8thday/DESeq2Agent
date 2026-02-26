import unittest
import inspect
from deseq2_agent.agents.base import BaseAgent
from deseq2_agent.pipeline import DESeq2Pipeline


class TestAgentProtocol(unittest.TestCase):
    """Verify agent and pipeline method signatures match current API."""

    def test_base_agent_has_invoke_methods(self):
        """BaseAgent exposes invoke/ainvoke (LCEL Runnable protocol)."""
        self.assertTrue(hasattr(BaseAgent, "invoke"))
        self.assertTrue(hasattr(BaseAgent, "ainvoke"))

    def test_base_agent_no_legacy_methods(self):
        """Legacy run/arun removed from BaseAgent."""
        self.assertFalse(hasattr(BaseAgent, "run"))
        self.assertFalse(hasattr(BaseAgent, "arun"))

    def test_base_agent_ainvoke_is_async(self):
        """BaseAgent.ainvoke is a coroutine function."""
        self.assertTrue(inspect.iscoroutinefunction(BaseAgent.ainvoke))

    def test_pipeline_has_run_method(self):
        """DESeq2Pipeline uses run() as orchestrator (not Runnable protocol)."""
        self.assertTrue(hasattr(DESeq2Pipeline, "run"))

    def test_pipeline_run_signature(self):
        """DESeq2Pipeline.run() accepts the expected parameters."""
        sig = inspect.signature(DESeq2Pipeline.run)
        params = list(sig.parameters.keys())
        self.assertIn("counts_file", params)
        self.assertIn("metadata_file", params)
        self.assertIn("contrasts", params)
        self.assertIn("species", params)
        self.assertIn("output_dir", params)

    def test_pipeline_contrasts_optional(self):
        """DESeq2Pipeline.run() contrasts parameter defaults to None."""
        sig = inspect.signature(DESeq2Pipeline.run)
        contrasts_param = sig.parameters["contrasts"]
        self.assertEqual(contrasts_param.default, None)

    def test_pipeline_has_design_agent(self):
        """DESeq2Pipeline.__init__ creates a design_agent attribute."""
        self.assertIn("design_agent", DESeq2Pipeline.__init__.__code__.co_names
                       + DESeq2Pipeline.__init__.__code__.co_varnames)


if __name__ == "__main__":
    unittest.main()
