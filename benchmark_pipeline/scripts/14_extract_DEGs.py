
"""
This script extracts the differentially expressed genes for each cluster of each tissue
"""
# pylint: disable=line-too-long
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

import pickle

import anndict as adt
import scanpy as sc
import pandas as pd

# Load the adata_dict
print("Loading adata_dict", flush=True)
adata_dict = adt.read_adata_dict('./dat/preprocessed_tissue_adt_ts2')
print("Loaded adata_dict", flush=True)

#Extract the differentially expressed genes for each cluster of each tissue
def extract_DEGs(adata, n_genes=10, adt_key=None):
    """
    Extract the differentially expressed genes for each cluster of each tissue
    """
    DEG_df = sc.get.rank_genes_groups_df(adata, None).groupby('group', observed=False).head(n_genes).reset_index(drop=True)
    if adt_key is not None:
        DEG_df['adt_key'] = [adt_key for _ in range(len(DEG_df))]
    return DEG_df

top_n_genes = adata_dict.fapply(extract_DEGs)
top_n_genes = pd.concat(top_n_genes.values(), ignore_index=True)

# Save the DEGs
print("Saving DEGs as csv", flush=True)
top_n_genes.to_csv('res/14_extract_DEGs/DEGs.csv', index=False)
print("Saved DEGs as csv", flush=True)

# Save the DEGs as a pickle file
print("Saving DEGs as pickle", flush=True)
with open('res/14_extract_DEGs/DEGs.pkl', 'wb') as f:
    pickle.dump(top_n_genes, f)
print("Saved DEGs as pickle", flush=True)

# Create output file to record DEG analysis
with open('res/14_extract_DEGs/DEG_frequencies.txt', 'w', encoding="utf-8") as out_file:
    # Output number of clusters analyzed (number of unique rows in DEG[["groups","adt_key"]]):
    cluster_count = len(top_n_genes[["group", "adt_key"]].drop_duplicates())
    out_file.write(f"Number of clusters analyzed: {cluster_count}\n\n")

    # Output number of occurrences of these genes
    genes = ["IGHM", "IGHD", "YBX3", "TNFRSF13B", "CD27", "ATXN1"]
    gene_counts = top_n_genes['names'].value_counts().reindex(genes, fill_value=0).astype(int)
    out_file.write("Number of occurrences for each gene:\n")
    for gene, count in gene_counts.items():
        out_file.write(f"{gene}: {count}\n")

print("DEG frequency analysis written to res/14_extract_DEGs/DEG_frequencies.txt")
