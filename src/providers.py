"""
This module contains the provider configure args and endpoint mappings for benchmark_pipeline.
"""

import os
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class ProviderConfig(dict):
    """Configuration settings for a language model provider"""
    provider: str
    api_key: str
    requests_per_minute: Optional[int] = None
    extra_params: dict[str, str] = field(default_factory=dict)

    def __init__(self, provider: str, api_key: str, requests_per_minute: Optional[int] = None, 
                 extra_params: dict[str, str] = None):
        # Initialize the dict
        super().__init__()

        # Set dataclass fields
        self.provider = provider
        self.api_key = api_key
        self.requests_per_minute = requests_per_minute
        self.extra_params = extra_params or {}

        # Populate dictionary values
        self['provider'] = provider
        self['api_key'] = api_key
        if requests_per_minute is not None:
            self['requests_per_minute'] = requests_per_minute

        # Add any extra parameters
        for key, value in self.extra_params.items():
            self[key] = value

@dataclass
class ProviderEndpoints:
    """Mapping of provider names to their available models"""
    endpoints: dict[str, list[str]]

# Define providers
PROVIDERS = {
    "openai": ProviderConfig("openai", os.getenv("OPENAI_API_KEY"), 9950),
    "anthropic": ProviderConfig("anthropic", os.getenv("ANTHROPIC_API_KEY"), 300),
    "google": ProviderConfig("google", os.getenv("GOOGLE_API_KEY"), 200),
    "huggingface": ProviderConfig("huggingface", os.getenv("HUGGINGFACE_API_KEY")),
    "azureml_endpoint": ProviderConfig(
        "azureml_endpoint", os.getenv("AZUREML_API_KEY"), None,
        {"endpoint_name": "your-endpoint-name", "region": "your-region"}
    ),
    "bedrock": ProviderConfig(
        "bedrock", os.getenv("BEDROCK_API_KEY"), 30,
        {"region_name": "us-west-2", "aws_access_key_id": os.getenv("AWS_ACCESS_KEY_ID")}
    ),
}

# Define model endpoints for each provider
ENDPOINTS = ProviderEndpoints({
    "bedrock": [
        "meta.llama3-1-8b-instruct-v1:0",
        "meta.llama3-1-70b-instruct-v1:0",
        "meta.llama3-1-405b-instruct-v1:0",
        "cohere.command-r-plus-v1:0",
        "mistral.mistral-large-2407-v1:0",
        # 'amazon.titan-text-express-v1',   # This model doesn't work well enough to include in benchmarking analysis (failed preliminary testing)
        # 'amazon.titan-text-lite-v1',      # This model doesn't work well enough to include in benchmarking analysis (failed preliminary testing)
        # 'ai21.j2-ultra-v1'                # This model doesn't work well enough to include in benchmarking analysis (failed preliminary testing)
    ],
    "google": [
        "gemini-1.5-pro",
        "gemini-1.5-flash"
        ],
    "openai": [
        "gpt-4",
        "gpt-4o",
        "gpt-4o-mini"
        ],
    "anthropic": [
        "claude-3-5-sonnet-20240620",
        "claude-3-opus-20240229",
        "claude-3-haiku-20240307",
    ],
})

MODEL_TICK_LABELS = {'claude-3-5-sonnet-20240620': 'Claude 3.5 Sonnet',
 'gemini-1.5-pro': 'Gemini 1.5 Pro',
 'claude-3-opus-20240229': 'Claude 3 Opus',
 'gpt-4o': 'GPT-4o',
 'gpt-4': 'GPT-4',
 'meta.llama3-1-405b-instruct-v1:0': 'Llama 3.1 405B Instruct',
 'claude-3-haiku-20240307': 'Claude 3 Haiku',
 'meta.llama3-1-70b-instruct-v1:0': 'Llama 3.1 70B Instruct',
 'mistral.mistral-large-2407-v1:0': 'Mistral Large',
 'gpt-4o-mini': 'GPT-4o mini',
 'cohere.command-r-plus-v1:0': 'Command R Plus',
 'gemini-1.5-flash': 'Gemini 1.5 Flash',
 'meta.llama3-1-8b-instruct-v1:0': 'Llama 3.1 8B Instruct'}

REMOVE_TICK_LABELS = {
    "claude-3-5-sonnet-20240620": "",
    "gemini-1.5-pro": "",
    "claude-3-opus-20240229": "",
    "gpt-4o": "",
    "gpt-4": "",
    "meta.llama3-1-405b-instruct-v1:0": "",
    "claude-3-haiku-20240307": "",
    "meta.llama3-1-70b-instruct-v1:0": "",
    "mistral.mistral-large-2407-v1:0": "",
    "gpt-4o-mini": "",
    "cohere.command-r-plus-v1:0": "",
    "gemini-1.5-flash": "",
    "meta.llama3-1-8b-instruct-v1:0": "",
}
