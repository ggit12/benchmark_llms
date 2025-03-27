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
    requests_per_minute: Optional[int] = None
    extra_params: dict[str, str] = field(default_factory=dict)

    def __init__(self, provider: str, requests_per_minute: Optional[int] = None,
                 extra_params: dict[str, str] = None):
        # Initialize the dict
        super().__init__()

        # Set dataclass fields
        self.provider = provider
        self.requests_per_minute = requests_per_minute
        self.extra_params = extra_params or {}

        # Populate dictionary values
        self['provider'] = provider
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
    "openai": ProviderConfig("openai", 8000, {"api_key": os.getenv("OPENAI_API_KEY")}),
    "anthropic": ProviderConfig("anthropic", 100, {"api_key": os.getenv("ANTHROPIC_API_KEY")}),
    "google": ProviderConfig("google", 200, {"api_key": os.getenv("GOOGLE_API_KEY")}),
    "huggingface": ProviderConfig("huggingface", {"api_key": os.getenv("HUGGINGFACE_API_KEY")}),
    "azureml_endpoint": ProviderConfig(
        "azureml_endpoint", None,
        {
            "api_key": os.getenv("AZUREML_API_KEY"),
            "endpoint_name": "your-endpoint-name",
            "region": "your-region"
        }
    ),
    "bedrock": ProviderConfig(
        "bedrock", 10,
        {
            "region_name": "us-west-2",
            "aws_access_key_id": os.getenv("AWS_ACCESS_KEY_ID"),
            "aws_secret_access_key": os.getenv("AWS_SECRET_ACCESS_KEY")
        }
    ),
}

# Define model endpoints for each provider
ENDPOINTS = ProviderEndpoints({
    "bedrock": [
        # "meta.llama3-1-8b-instruct-v1:0", # These models disabled temporarily due to rate limit issues with the Bedrock API (too low)
        # "meta.llama3-1-70b-instruct-v1:0",
        # "meta.llama3-1-405b-instruct-v1:0",
        # "cohere.command-r-plus-v1:0",
        # "mistral.mistral-large-2407-v1:0", # End of disabled models
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
        "gpt-4o-mini",
        # "o1-mini-2024-09-12", # o-series models not yet supported in AnnDictionary
        # "o3-mini-2025-01-31",
        ],
    "anthropic": [
        "claude-3-7-sonnet-20250219",
        "claude-3-5-sonnet-20240620",
        "claude-3-5-haiku-20241022",
        "claude-3-opus-20240229",
        "claude-3-haiku-20240307",
    ],
})

MODEL_TICK_LABELS = {
 # Anthropic
 'claude-3-7-sonnet-20250219': 'Claude 3.7 Sonnet',
 'claude-3-5-sonnet-20240620': 'Claude 3.5 Sonnet',
 'claude-3-5-haiku-20241022': 'Claude 3.5 Haiku',
 'claude-3-opus-20240229': 'Claude 3 Opus',
 'claude-3-haiku-20240307': 'Claude 3 Haiku',

 # OpenAI
 'gpt-4o': 'GPT-4o',
 'gpt-4': 'GPT-4',
 'gpt-4o-mini': 'GPT-4o mini',
 'o1-mini-2024-09-12': 'o1 mini',
 'o3-mini-2025-01-31': 'o3 mini',

 # Google
 'gemini-1.5-flash': 'Gemini 1.5 Flash',
 'gemini-1.5-pro': 'Gemini 1.5 Pro',

 # Bedrock
 'meta.llama3-1-405b-instruct-v1:0': 'Llama 3.1 405B Instruct',
 'meta.llama3-1-70b-instruct-v1:0': 'Llama 3.1 70B Instruct',
 'meta.llama3-1-8b-instruct-v1:0': 'Llama 3.1 8B Instruct',
 'mistral.mistral-large-2407-v1:0': 'Mistral Large',
 'cohere.command-r-plus-v1:0': 'Command R Plus',
 }

REMOVE_TICK_LABELS = {key: "" for key in MODEL_TICK_LABELS}
# REMOVE_TICK_LABELS = {
#     # Anthropic
#     "claude-3-7-sonnet-20250219": "",
#     "claude-3-5-sonnet-20240620": "",
#     "claude-3-5-haiku-20241022": "",
#     "claude-3-opus-20240229": "",
#     "claude-3-haiku-20240307": "",

#     # OpenAI
#     "gpt-4o": "",
#     "gpt-4": "",
#     "gpt-4o-mini": "",
#     "o1-mini-2024-09-12": "",
#     "o3-mini-2025-01-31": "",


#     # Google
#     "gemini-1.5-pro": "",
#     "gemini-1.5-flash": "",

#     # Bedrock
#     "meta.llama3-1-405b-instruct-v1:0": "",
#     "meta.llama3-1-70b-instruct-v1:0": "",
#     "meta.llama3-1-8b-instruct-v1:0": "",
#     "mistral.mistral-large-2407-v1:0": "",
#     "cohere.command-r-plus-v1:0": "",
# }
