"""
This script takes model, input, and output. This script runs analysis on input using the specified model and saves to output.
"""

import os
import sys
import pickle
import argparse
import anndict as adt

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

from src import PROVIDERS, ENDPOINTS, run_multiple_providers_models # pylint: disable=wrong-import-position

parser = argparse.ArgumentParser()
parser.add_argument("--provider", required=True, help="Name of the provider")
parser.add_argument("--input", required=True, help="Path to input adata_dict directory")
parser.add_argument("--outdir", required=True, help="Output directory")
args = parser.parse_args()

adata_dict = adt.read_adata_dict(args.input)

# Get the list of models for this provider
all_models = ENDPOINTS.endpoints[args.provider]

# Check outdir to see which models have already been run
completed_models = [
    f.replace(".pkl", "").replace(f"{args.provider}_", "")
    for f in os.listdir(args.outdir)
    if args.provider in f
]
model_list = [m for m in all_models if m not in completed_models]

# Print information for debugging
print(f"All models: {all_models}")
print(f"Files in directory: {os.listdir(args.outdir)}")
print(f"Completed models: {completed_models}")
print(f"Models to run: {model_list}")

#End if all models have been run
if not model_list:
    print(f"All models for {args.provider} have already been run.")
    sys.exit(0)

# Create a dict mapping this provider to its models
provider_endpoint_dict = {args.provider: model_list}

# Create a dict for provider config using the global PROVIDERS dictionary
provider_config = {args.provider: PROVIDERS[args.provider]}

results = run_multiple_providers_models(
    adata_dict, provider_endpoint_dict, provider_config
)

os.makedirs(args.outdir, exist_ok=True)

# Write out separate pickle file for each model
for model in model_list:
    model_key = f"{args.provider}_{model}"
    model_results = results.get(model_key, {})
    pickle_path = os.path.join(args.outdir, f"{model_key}.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(model_results, f)
