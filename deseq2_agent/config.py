"""LLM Configuration and Initialization for DESeq2Agent.

Supported providers:
  - openai:             OpenAI API (gpt-4o, gpt-4o-mini, etc.)
  - deepseek:           DeepSeek chat models via ChatDeepSeek (deepseek-chat)
  - deepseek-reasoner:  DeepSeek reasoning model via ChatDeepSeek (deepseek-reasoner)
  - anthropic:          Anthropic Claude (requires langchain-anthropic)
  - google:             Google Gemini (requires langchain-google-genai)
  - <any>:              Any provider supported by langchain init_chat_model()

DeepSeek models use langchain-deepseek's ChatDeepSeek which handles:
  - reasoning_content extraction for reasoner models
  - disabled_params for tool_choice (reasoner doesn't support forced tool invocation)

OpenAI uses ChatOpenAI. Other providers use init_chat_model() for auto-routing.
"""

import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from dotenv import load_dotenv
from langchain_core.language_models import BaseChatModel
from langchain_openai import ChatOpenAI

# Load environment variables from .env file
load_dotenv()

logger = logging.getLogger(__name__)

# ── Provider presets ─────────────────────────────────────────────────────────
# Each preset defines: env var names, defaults, and any special kwargs.
# "model_provider" is used by init_chat_model for non-OpenAI-compatible providers.

PROVIDER_PRESETS: Dict[str, Dict[str, Any]] = {
    "openai": {
        "env_key": "OPENAI_API_KEY",
        "env_model": "OPENAI_MODEL",
        "default_model": "gpt-4o",
    },
    "deepseek": {
        "env_key": "DEEPSEEK_API_KEY",
        "env_model": "DEEPSEEK_MODEL",
        "default_model": "deepseek-chat",
    },
    "deepseek-reasoner": {
        "env_key": "DEEPSEEK_API_KEY",
        "env_model": "DEEPSEEK_MODEL",
        "default_model": "deepseek-reasoner",
        "temperature": None,  # reasoner does not accept temperature
    },
    "anthropic": {
        "env_key": "ANTHROPIC_API_KEY",
        "env_model": "ANTHROPIC_MODEL",
        "default_model": "claude-sonnet-4-20250514",
        "model_provider": "anthropic",  # signals: use init_chat_model
    },
    "google": {
        "env_key": "GOOGLE_API_KEY",
        "env_model": "GOOGLE_MODEL",
        "default_model": "gemini-2.0-flash",
        "model_provider": "google_genai",
    },
}

# Providers that use ChatDeepSeek
_DEEPSEEK_PROVIDERS = {"deepseek", "deepseek-reasoner"}


@dataclass
class LLMConfig:
    """Configuration for LLM initialization."""

    api_key: Optional[str] = None
    model: str = "gpt-4o"
    temperature: Optional[float] = 0.2     # None = omit (required for reasoner)
    max_retries: int = 3
    timeout: float = 60.0
    base_url: Optional[str] = None         # For OpenAI-compatible endpoints
    model_provider: Optional[str] = None   # For init_chat_model routing
    extra_kwargs: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_env(cls, provider: str = "openai") -> "LLMConfig":
        """Create config from environment variables.

        Args:
            provider: Provider name — one of PROVIDER_PRESETS keys,
                      or any provider supported by langchain init_chat_model.

        Returns:
            LLMConfig instance
        """
        preset = PROVIDER_PRESETS.get(provider, {})

        env_key = preset.get("env_key", f"{provider.upper()}_API_KEY")
        env_model = preset.get("env_model", f"{provider.upper()}_MODEL")
        default_model = preset.get("default_model", "gpt-4o")

        # Temperature: preset can override (e.g. None for reasoner)
        if "temperature" in preset:
            temperature = preset["temperature"]
        else:
            temperature = float(os.getenv("LLM_TEMPERATURE", "0.2"))

        return cls(
            api_key=os.getenv(env_key),
            model=os.getenv(env_model, default_model),
            temperature=temperature,
            max_retries=int(os.getenv("LLM_MAX_RETRIES", "3")),
            base_url=preset.get("base_url"),
            model_provider=preset.get("model_provider"),
        )


def _build_deepseek(config: LLMConfig) -> BaseChatModel:
    """Build a ChatDeepSeek instance with proper reasoner handling."""
    from langchain_deepseek import ChatDeepSeek

    ds_kwargs: Dict[str, Any] = {
        "model": config.model,
        "api_key": config.api_key,
        "max_retries": config.max_retries,
        **config.extra_kwargs,
    }
    if config.temperature is not None:
        ds_kwargs["temperature"] = config.temperature

    # deepseek-reasoner supports function calling (since R1-0528) but does NOT
    # support forced tool invocation (tool_choice). Disable it so that
    # with_structured_output(method="function_calling") works correctly.
    if "reasoner" in config.model.lower():
        ds_kwargs["disabled_params"] = {"tool_choice": None}
        logger.info(
            f"Reasoner model '{config.model}': disabled tool_choice for "
            "structured output compatibility"
        )

    logger.info(f"Initializing ChatDeepSeek: model={config.model}")
    return ChatDeepSeek(**ds_kwargs)


def get_llm(
    config: Optional[LLMConfig] = None,
    provider: str = "openai",
    **kwargs,
) -> BaseChatModel:
    """Get a configured LangChain chat model instance.

    Provider routing:
      - openai            → ChatOpenAI
      - deepseek*         → ChatDeepSeek (handles reasoning_content + disabled_params)
      - anthropic/google  → init_chat_model() auto-routing

    Args:
        config: Optional LLMConfig. If not provided, loads from environment.
        provider: Provider name (see PROVIDER_PRESETS or any langchain provider)
        **kwargs: Additional arguments passed to the chat model constructor

    Returns:
        Configured BaseChatModel instance

    Examples:
        >>> llm = get_llm()                              # OpenAI (default)
        >>> llm = get_llm(provider="deepseek")            # DeepSeek chat
        >>> llm = get_llm(provider="deepseek-reasoner")   # DeepSeek reasoner
        >>> llm = get_llm(provider="anthropic")           # Claude
        >>> llm = get_llm(provider="google")              # Gemini
    """
    if config is None:
        config = LLMConfig.from_env(provider)

    # Merge extra kwargs
    if kwargs:
        config.extra_kwargs.update(kwargs)

    if not config.api_key:
        preset = PROVIDER_PRESETS.get(provider, {})
        env_key = preset.get("env_key", f"{provider.upper()}_API_KEY")
        raise ValueError(
            f"API key not found for provider '{provider}'. "
            f"Please set {env_key} environment variable or pass api_key in LLMConfig."
        )

    # ── DeepSeek providers → ChatDeepSeek ─────────────────────────────────
    if provider in _DEEPSEEK_PROVIDERS:
        return _build_deepseek(config)

    # ── OpenAI → ChatOpenAI ───────────────────────────────────────────────
    if provider == "openai" or config.model_provider is None:
        llm_kwargs: Dict[str, Any] = {
            "api_key": config.api_key,
            "model": config.model,
            "max_retries": config.max_retries,
            "timeout": config.timeout,
            **config.extra_kwargs,
        }
        if config.temperature is not None:
            llm_kwargs["temperature"] = config.temperature
        if config.base_url:
            llm_kwargs["base_url"] = config.base_url

        logger.info(f"Initializing ChatOpenAI: model={config.model}")
        return ChatOpenAI(**llm_kwargs)

    # ── Other providers → init_chat_model ─────────────────────────────────
    try:
        from langchain.chat_models import init_chat_model
    except ImportError:
        raise ImportError(
            f"Provider '{provider}' requires the 'langchain' package and the "
            f"corresponding provider package. "
            f"Install with: pip install langchain langchain-{config.model_provider}"
        )

    init_kwargs: Dict[str, Any] = {
        "api_key": config.api_key,
        "max_retries": config.max_retries,
        **config.extra_kwargs,
    }
    if config.temperature is not None:
        init_kwargs["temperature"] = config.temperature

    logger.info(
        f"Initializing via init_chat_model: model={config.model}, "
        f"provider={config.model_provider}"
    )
    return init_chat_model(
        config.model,
        model_provider=config.model_provider,
        **init_kwargs,
    )
