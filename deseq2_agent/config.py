"""LLM Configuration and Initialization for DESeq2Agent"""

import os
from dataclasses import dataclass
from typing import Optional

from dotenv import load_dotenv
from langchain_openai import ChatOpenAI

# Load environment variables from .env file
load_dotenv()


@dataclass
class LLMConfig:
    """Configuration for LLM initialization."""
    
    api_key: Optional[str] = None
    model: str = "gpt-4o"
    temperature: float = 0.2
    max_retries: int = 3
    timeout: float = 60.0
    base_url: Optional[str] = None  # For DeepSeek or other providers
    
    @classmethod
    def from_env(cls, provider: str = "openai") -> "LLMConfig":
        """Create config from environment variables.
        
        Args:
            provider: LLM provider ("openai" or "deepseek")
            
        Returns:
            LLMConfig instance
        """
        if provider == "deepseek":
            return cls(
                api_key=os.getenv("DEEPSEEK_API_KEY"),
                model=os.getenv("DEEPSEEK_MODEL", "deepseek-chat"),
                temperature=float(os.getenv("LLM_TEMPERATURE", "0.2")),
                max_retries=int(os.getenv("LLM_MAX_RETRIES", "3")),
                base_url=os.getenv("DEEPSEEK_BASE_URL", "https://api.deepseek.com/v1"),
            )
        else:  # openai
            return cls(
                api_key=os.getenv("OPENAI_API_KEY"),
                model=os.getenv("OPENAI_MODEL", "gpt-4o"),
                temperature=float(os.getenv("LLM_TEMPERATURE", "0.2")),
                max_retries=int(os.getenv("LLM_MAX_RETRIES", "3")),
            )


def get_llm(
    config: Optional[LLMConfig] = None,
    provider: str = "openai",
    **kwargs
) -> ChatOpenAI:
    """Get a configured LangChain LLM instance.
    
    Args:
        config: Optional LLMConfig. If not provided, loads from environment.
        provider: LLM provider ("openai" or "deepseek")
        **kwargs: Additional arguments passed to ChatOpenAI
        
    Returns:
        Configured ChatOpenAI instance
        
    Example:
        >>> llm = get_llm()  # Uses OPENAI_API_KEY from environment
        >>> llm = get_llm(provider="deepseek")  # Uses DeepSeek
        >>> llm = get_llm(LLMConfig(api_key="sk-...", model="gpt-4o-mini"))
    """
    if config is None:
        config = LLMConfig.from_env(provider)
    
    if not config.api_key:
        raise ValueError(
            f"API key not found. Please set {'DEEPSEEK_API_KEY' if provider == 'deepseek' else 'OPENAI_API_KEY'} "
            "environment variable or pass api_key in LLMConfig."
        )
    
    llm_kwargs = {
        "api_key": config.api_key,
        "model": config.model,
        "temperature": config.temperature,
        "max_retries": config.max_retries,
        "timeout": config.timeout,
        **kwargs,
    }
    
    # Add base_url for non-OpenAI providers (e.g., DeepSeek)
    if config.base_url:
        llm_kwargs["base_url"] = config.base_url
    
    return ChatOpenAI(**llm_kwargs)
