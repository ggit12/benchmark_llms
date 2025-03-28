"""
This script takes model, input, and output. This script runs analysis on input using the specified model and saves to output.
"""
# pylint: disable=invalid-name

import os
import sys
import pickle
import argparse
import shutil

import anndict as adt
import pandas as pd

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

# Get the list of models for this provider
all_models = ENDPOINTS.endpoints[args.provider]

# Check if cache directory exists and move files back to outdir
cache_dir = os.path.join(args.outdir, "cached_outputs")
if os.path.exists(cache_dir):
    cache_files = [f for f in os.listdir(cache_dir) if f.startswith(args.provider)]
    for file in cache_files:
        source = os.path.join(cache_dir, file)
        dest = os.path.join(args.outdir, file)
        shutil.copy2(source, dest)  # preserves metadata
        os.remove(source)  # remove from cache

# Check outdir to see which models have already been run
os.makedirs(args.outdir, exist_ok=True)
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
    #Write a done file to indicate that the script has completed
    # with open(os.path.join(args.outdir, f"{args.provider}.done"), "w") as f:
    #     f.write("done")
    print(f"All models for {args.provider} have already been run.", flush=True)
    sys.exit(0)

# Create a dict mapping this provider to its models
provider_endpoint_dict = {args.provider: model_list}

# Create a dict for provider config using the global PROVIDERS dictionary
provider_config = {args.provider: PROVIDERS[args.provider]}

adata_dict = adt.read_adata_dict(args.input)

results = run_multiple_providers_models(
    adata_dict, provider_endpoint_dict, provider_config
)

# Write out separate pickle file for each model
# Before writing, the output is validated.
# If the output is not as expected, writing is skipped for that model's output.
all_models_completed = True
for model in model_list:
    model_key = f"{args.provider}_{model}"
    model_results = results.get(model_key, {})

    # Check each DataFrame in the dictionary
    this_model_completed = True
    for label, df in model_results['label_results'].items():
        if not isinstance(df, pd.DataFrame):
            print(f"Result for label '{label}' in {model_key} is not a pd.DataFrame.", flush=True)
            this_model_completed = False

        cell_type_col = [col for col in df.columns if col.endswith('_ai_cell_type')]
        if cell_type_col and df[cell_type_col].isna().any().any():
            print(f"Results for label '{label}' in {model_key} has NaN values in the '_ai_cell_type' column.", flush=True)
            this_model_completed = False

    if not this_model_completed:
        all_models_completed = False
        continue

    # Save results to pickle file
    pickle_path = os.path.join(args.outdir, f"{model_key}.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(model_results, f)


if not all_models_completed:
    # Exit with error if not all models have been run successfully
    raise ValueError("Some models did not complete successfully. Check logs for details and rerun this rule.")
# else:
#     #Write a done file to indicate that the script has completed
#     with open(os.path.join(args.outdir, f"{args.provider}.done"), "w") as f:
#         f.write("done")
