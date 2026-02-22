"""Base Agent Class for DESeq2Agent

Defines the abstract interface that all specialized agents must implement.
"""


import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Type, Optional

from langchain_core.language_models import BaseChatModel
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnableSequence
from pydantic import BaseModel

logger = logging.getLogger(__name__)


class BaseAgent(ABC):
    """Abstract base class for all DESeq2 analysis agents.
    
    Each agent is responsible for a specific aspect of RNA-seq analysis:
    - DataAgent: Experimental design understanding
    - QCAgent: Quality control assessment
    - DEAgent: Differential expression interpretation
    - PathwayAgent: Pathway analysis
    - ReportAgent: Report generation
    
    Agents use LangChain's LCEL for composable, type-safe chains.
    """
    
    # Subclasses must define these
    prompt_template: ChatPromptTemplate
    output_model: Type[BaseModel]
    
    def __init__(self, llm: BaseChatModel):
        """Initialize the agent with an LLM.
        
        Args:
            llm: LangChain chat model instance
        """
        self.llm = llm
        self._chain = self._build_chain()
    
    def _build_chain(self) -> RunnableSequence:
        """Build the LCEL chain for this agent.
        
        Returns:
            A runnable chain: prompt | llm | parser
        """
        # Use function_calling for cross-provider compatibility (DeepSeek, OpenAI, etc.)
        structured_llm = self.llm.with_structured_output(
            self.output_model, method="function_calling"
        )
        
        # Build the chain
        chain = self.prompt_template | structured_llm
        
        return chain
    
    @abstractmethod
    def validate_input(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Validate and preprocess input payload.
        
        Args:
            payload: Input data dictionary
            
        Returns:
            Validated/preprocessed payload
            
        Raises:
            ValueError: If required fields are missing
        """
        pass
    
    def invoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> BaseModel:
        """Execute the agent with the given input.
        
        Args:
            input: Input fields for the prompt template
            config: Optional configuration for the runnable
            
        Returns:
            Pydantic model instance with structured output
        """
        # Validate input
        validated_input = self.validate_input(input)
        
        # Fill in missing optional fields with empty strings
        for key in self.prompt_template.input_variables:
            if key not in validated_input:
                validated_input[key] = ""
        
        logger.info(f"Running {self.__class__.__name__}...")
        
        # Execute the chain
        result = self._chain.invoke(validated_input, config=config)
        
        logger.info(f"{self.__class__.__name__} completed successfully")
        
        return result
    
    async def ainvoke(self, input: Dict[str, Any], config: Optional[Any] = None) -> BaseModel:
        """Async version of invoke().
        
        Args:
            input: Input fields for the prompt template
            config: Optional configuration for the runnable
            
        Returns:
            Pydantic model instance with structured output
        """
        validated_input = self.validate_input(input)
        
        for key in self.prompt_template.input_variables:
            if key not in validated_input:
                validated_input[key] = ""
        
        logger.info(f"Running {self.__class__.__name__} (async)...")
        
        result = await self._chain.ainvoke(validated_input, config=config)
        
        logger.info(f"{self.__class__.__name__} completed successfully")
        
        return result
