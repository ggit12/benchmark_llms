"""
Calculates the direct string agreement between the cell type labels 
and saves the results into the 04_postprocess_results directory.
"""
# pylint: disable=line-too-long
# pylint: disable=invalid-name

import pickle
import os
import sys

import anndict as adt
import scanpy as sc

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    direct_compare_cell_type_labels_pairwise,
)

# Read the merged AnnData
adata = sc.read_h5ad('./res/04_postprocess_results/adt_de_novo_llm_annotated.h5ad')

# Load the cell type columns
base_path = './res/04_postprocess_results/'
with open(base_path + 'llm_celltype_cols.pkl', 'rb') as f:
    llm_celltype_cols = pickle.load(f)

# Read manual cell type column
with open("../../dat/manual_cell_type_col.pkl", 'rb') as f:
    manual_cell_type_col = pickle.load(f)

#assess the agreement using direct string comparison
consistent_manual_cell_type_col = 'consistent_including_manual_' + manual_cell_type_col
direct_compare_cell_type_labels_pairwise(adata, [consistent_manual_cell_type_col], llm_celltype_cols, new_col_prefix='direct_string_agreement')
print("Calculated direct string agreement", flush=True)

# get these column names
direct_string_agreement_cols = adt.get_adata_columns(adata, contains=['direct_string_agreement_consistent_including_manual'])

# Write out all generated objects
adata.write_h5ad(base_path + 'adt_de_novo_llm_annotated.h5ad')
pickle.dump(direct_string_agreement_cols, open(base_path + 'direct_string_agreement_cols.pkl', 'wb'))
