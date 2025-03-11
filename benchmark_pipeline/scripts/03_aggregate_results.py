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
import sys

import anndict as adt

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    cell_type_by_plurality,
    PROVIDERS,
)


# Load the results
results_paths = [f for f in glob.glob('./res/02_run_provider/*.pkl')]

results = {}
for path in results_paths:
    with open(path, 'rb') as f:
        model_results = pickle.load(f)
        model_name = os.path.basename(path).replace('.pkl', '')
        results[model_name] = model_results

# Load the adata_dict
adata_dict = adt.read_adata_dict('./dat/preprocessed_tissue_adt_ts2')

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

    # return adata_dict

# Assuming adata_dict and results are already defined and populated
# adata_dict = merge_obs_to_adata(adata_dict, results)
merge_obs_to_adata(adata_dict, results)
print("Merged obs columns into adata_dict")

# Merge the adata_dict
adata = adt.concatenate_adata_dict(adata_dict)
print("Concatenated adata_dict")

# After merging the AnnData, the labels are then processed so that various metrics can be calculated using them.

# Process the labels
#Configure the backend to work with a specific provider and model
provider = os.environ['PROVIDER_FOR_DOWNSTREAM_ANALYSIS']
llm_config = {provider: PROVIDERS[provider]}
llm_config['model'] = os.environ['MODEL_FOR_DOWNSTREAM_ANALYSIS']

adt.configure_llm_backend(**llm_config)

# adt.configure_llm_backend(provider='anthropic',
#                           model='claude-3-5-sonnet-20240620',
#                           api_key=os.environ['ANTHROPIC_API_KEY'],
#                           requests_per_minute=400
#                           )
# adt.configure_llm_backend(provider='openai',
#                           model='gpt-4o',
#                           api_key=os.environ['OPENAI_API_KEY'],
#                           requests_per_minute=9950)


#get cell type columns
# cell_type_cols = adt.get_adata_columns(adata, col_endswith=['ai_cell_sub_type', 'simplified_ai_cell_type'], not_col_startswith=['raw', 'agreement'])
cell_type_cols = adt.get_adata_columns(adata, ends_with=['simplified_ai_cell_type'], not_starts_with=['raw', 'agreement'], not_contains=['consistent'])

# Define manual cell type column
manual_cell_type_col = 'cell_ontology_class'

#unify category labels across all ai annotations and manual annotation.
label_map_with_manual = adt.ensure_label_consistency_adata(adata, cell_type_cols + [manual_cell_type_col], simplification_level='unified', new_col_prefix='consistent_including_manual')

#get unified cols
unified_cell_types_with_manual = adt.get_adata_columns(adata, contains=['consistent_including_manual'])



#calculate a cell type by majority vote of all the LLMs
llm_celltype_cols = adt.get_adata_columns(adata, contains=['consistent_including_manual'], not_contains=[manual_cell_type_col])
cell_type_by_plurality(adata, rater_cols=llm_celltype_cols, new_col_name='cell_type_by_plurality')
print("Calculated cell type by plurality")


#assess the agreement between the 'consistified' manual annotations and the ai-generated annotations
label_agreement_binary = adt.ai_compare_cell_type_labels_pairwise(adata, ['consistent_including_manual_' + manual_cell_type_col], llm_celltype_cols, new_col_prefix='binary_agreement', comparison_level='binary')
print("Calculated binary agreement")

#get these column names
binary_agreement_cols = adt.get_adata_columns(adata, contains = ['binary_agreement_consistent_including_manual'])

#also assess at the partial agreement level
label_agreement_categorical = adt.ai_compare_cell_type_labels_pairwise(adata, ['consistent_including_manual_' + manual_cell_type_col], llm_celltype_cols, new_col_prefix='categorical_agreement', comparison_level='categorical')
print("Calculated categorical agreement")

#get these column names
categorical_agreement_cols = adt.get_adata_columns(adata, contains = ['categorical_agreement_consistent_including_manual'])

#change scale of categorical label agreement cols from 0 to 1 (raw has values 0, 1, and 2)
adata.obs[categorical_agreement_cols] = adata.obs[categorical_agreement_cols]/2

#use categorical labels to count only perfect matches, and set partial matches to 0
perfect_only_categorical_agreement_cols = ["perfect_only_" + col for col in categorical_agreement_cols]
adata.obs[perfect_only_categorical_agreement_cols] = adata.obs[categorical_agreement_cols].replace(0.5, 0)
print("Calculated perfect only categorical agreement")

# Write out all generated objects
base_path = './res/03_gather_outputs/'

#adata_dict
# pickle.dump(adata_dict, open(base_path + "ts2_de_novo_llm_annotated_adt", 'wb'))

#adata
del adata.obs["adt_key"] # can't write a tuple column in obs of anndata
adata.write_h5ad('./res/03_gather_outputs/ts2_de_novo_llm_annotated.h5ad')
print("Wrote adata")

#label_map_with_manual
pickle.dump(label_map_with_manual, open(base_path + 'label_map_with_manual.pkl', 'wb'))

#label_agreement_binary
pickle.dump(label_agreement_binary, open(base_path + 'label_agreement_binary.pkl', 'wb'))

#label_agreement_categorical
pickle.dump(label_agreement_categorical, open(base_path + 'label_agreement_categorical.pkl', 'wb'))

#Write all these
pickle.dump(manual_cell_type_col, open(base_path + 'manual_cell_type_col.pkl', 'wb'))
pickle.dump(llm_celltype_cols, open(base_path + 'llm_celltype_cols.pkl', 'wb'))
pickle.dump(binary_agreement_cols, open(base_path + 'binary_agreement_cols.pkl', 'wb'))
pickle.dump(categorical_agreement_cols, open(base_path + 'categorical_agreement_cols.pkl', 'wb'))
pickle.dump(perfect_only_categorical_agreement_cols, open(base_path + 'perfect_only_categorical_agreement_cols.pkl', 'wb'))
print("Wrote all outputs")
