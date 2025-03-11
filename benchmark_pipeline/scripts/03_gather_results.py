"""
This script concatenates the results of cell type annotation. 
Then, annotation columns are made that map to a shared set of 
labels for comparison. Comparisons are then also run.
"""
# pylint: disable=line-too-long
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

import pickle
import os
import glob
import gc

import anndict as adt

# Load the results
print("Loading results", flush=True)
results_paths = [f for f in glob.glob('./res/02_run_provider/*.pkl')]

results = {}
for path in results_paths:
    with open(path, 'rb') as f:
        model_results = pickle.load(f)
        model_name = os.path.basename(path).replace('.pkl', '')
        results[model_name] = model_results
print("Loaded results", flush=True)

# Load the adata_dict
print("Loading adata_dict", flush=True)
adata_dict = adt.read_adata_dict('./dat/preprocessed_tissue_adt_ts2')
print("Loaded adata_dict", flush=True)

def merge_obs_to_adata(adata_dict, results):
    """Merges obs columns that have annotations into adata_dict"""
    for organ in adata_dict.keys():
        for model_key in results.keys():
            if 'obs_dict' in results[model_key] and organ in results[model_key]['obs_dict']:
                obs_df = results[model_key]['obs_dict'][organ]

                # Identify new columns
                new_cols = [col for col in obs_df.columns if col not in adata_dict[organ].obs.columns]

                if new_cols:
                    # Assume the index of obs_df matches the index of adata_dict[organ].obs
                    # If not, you might need to specify a column to use as the index
                    merged_obs = adata_dict[organ].obs.join(obs_df[new_cols], how='left')

                    # Update the obs DataFrame
                    adata_dict[organ].obs = merged_obs

# Assuming adata_dict and results are already defined and populated
# adata_dict = merge_obs_to_adata(adata_dict, results)
print("Merging obs columns into adata_dict", flush=True)
merge_obs_to_adata(adata_dict, results)
print("Merged obs columns into adata_dict", flush=True)

# Free up memory
del results
gc.collect()
print("deleted ``results`` to free memory", flush=True)

# Merge the adata_dict
adata = adt.concatenate_adata_dict(adata_dict, new_col_name=None) # To achieve future default behavior, setting new_col_name to None
print("Concatenated adata_dict", flush=True)

# Free up memory
del adata_dict
gc.collect()
print("deleted ``adata_dict`` to free memory", flush=True)

# Write out all generated objects
base_path = './res/03_gather_results/'

#adata_dict
# pickle.dump(adata_dict, open(base_path + "ts2_de_novo_llm_annotated_adt", 'wb'))

#adata
# del adata.obs["adt_key"] # can't write a tuple column in obs of anndata
adata.write_h5ad('./res/03_gather_results/ts2_de_novo_llm_annotated.h5ad')
print("Wrote adata", flush=True)
