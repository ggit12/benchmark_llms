"""
This script preprocesses the TSPv2 dataset for annotation benchmarking.

Runs from the root of the benchmark_pipeline directory.
"""
# pylint: disable=line-too-long

import os
import gc
import csv
import sys
import pickle

import pkg_resources
import anndict as adt
import scanpy as sc
import psutil


from dotenv import load_dotenv
load_dotenv()

#read data
adata = sc.read_h5ad(os.environ["INPUT_DATA"]) 

#remove extra obsm and layers for memory purposes
# Delete all the extra data
for key in ['X_pca', 'X_scvi', 'X_umap', 'X_umap_scvi_full_donorassay', 'X_uncorrected_alltissues_umap', 'X_uncorrected_umap']:
    if key in adata.obsm:
        del adata.obsm[key]

for key in ['decontXcounts', 'log_normalized', 'scale_data']:
    if key in adata.layers:
        del adata.layers[key]

#set X to be raw_counts
adata.X = adata.layers['raw_counts'].copy()

# Take only protein coding genes
source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

#load list of protein-coding genes
protein_coding_path = pkg_resources.resource_filename('src', 'dat/protein_coding_genes_list.csv')
with open(protein_coding_path, 'r', encoding='utf-8') as file:
    protein_coding = list(csv.reader(file, delimiter=','))
protein_coding = [i[6] for i in protein_coding[1:] if i[6]] #extra if at the end removes empty genes (i.e. "") from the list

#filter adata to only include protein coding genes
adata = adata[:, adata.var.index.isin(protein_coding)]

#build adata_dict
adata_dict = adt.build_adata_dict(adata, ['tissue'])


#remove a standard list of uninformative genes
abundant_rnas = [
    "MALAT1",
    "NEAT1",
    "XIST",
    "KCNQ1OT1",
    "RPPH1",
    "RN7SL1",
    "RMRP",
    "SNHG1",
    "MIAT",
    "H19"
]

adt.wrappers.remove_genes_adata_dict(adata_dict, abundant_rnas)


#free memory
del adata
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

#Need to set use_multithreading to False to avoid overloading the memory

#Run leiden clustering on each adata independently
#adata.X is raw counts, so run standard preprocessing
# Normalize each AnnData in the dictionary
adt.wrappers.normalize_adata_dict(adata_dict, use_multithreading=False)
print("Normalized")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Log transform each AnnData in the dictionary
adt.wrappers.log_transform_adata_dict(adata_dict, use_multithreading=False)
print("Log transformed")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Optionally, you might subset the data to only high-variance genes
adt.wrappers.set_high_variance_genes_adata_dict(adata_dict, n_top_genes=2000, subset=False, use_multithreading=False)
print("Set high variance genes")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Scale each AnnData in the dictionary
adt.wrappers.scale_adata_dict(adata_dict, use_multithreading=False)
print("Scaled")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Perform PCA on each AnnData in the dictionary
adt.wrappers.pca_adata_dict(adata_dict, n_comps=50, mask_var='highly_variable', use_multithreading=False)
print("Performed PCA")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

#Calculate the neighborhood graph
adt.wrappers.neighbors_adata_dict(adata_dict) # Already disables multithreading because sc.pp.neighbors is multithreaded
print("Calculated neighborhood graph")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

#Calculate the UMAP
adt.wrappers.calculate_umap_adata_dict(adata_dict) # Same as above
print("Calculated UMAP")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

#get leiden clusters
adt.wrappers.leiden_adata_dict(adata_dict, resolution=0.5, use_multithreading=False)
print("Ran leiden")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Set the var index to the gene name because the index is what will be used by the LLM
# adata_dict.set_var_index('feature_name')

# Remove clusters with fewer than 5 cells so that DEGs can be calculated
# adata_dict = adt.sample_and_drop_adata_dict(adata_dict, strata_keys=['leiden'], min_num_cells=5)

# Run diffexp analysis
adt.wrappers.rank_genes_groups_adata_dict(adata_dict, groupby='leiden', use_multithreading=False)
print("Ran diffexp")
print("Memory before deletion:", psutil.virtual_memory())
gc.collect()
print("Memory after deletion:", psutil.virtual_memory())

# Write the preprocessed AdataDict
adt.write_adata_dict(adata_dict, "./dat/preprocessed_tissue_adt_ts2")
print("Wrote preprocessed AdataDict")

# Write the manual cell type column as pickle
manual_cell_type_col = 'cell_ontology_class' # pylint: disable=invalid-name
with open("./dat/manual_cell_type_col.pkl", 'wb') as f:
    pickle.dump(manual_cell_type_col, f)
print("Wrote manual cell type column")
