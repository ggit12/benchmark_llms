"""
A function to plot genes on a UMAP.
"""

import matplotlib.pyplot as plt
import scanpy as sc


def plot_genes_umap(adata, genes, n_cols=4, fig_size=(5, 5), layer=None, title=None):
    """Plots the genes on a UMAP using scanpy's sc.pl.umap function."""

    # Filter genes to those present in the data
    genes_to_plot = [gene for gene in genes if gene in adata.var_names]

    if not genes_to_plot:
        print("None of the specified genes are present in the data.")
        return

    # Calculate number of rows based on number of genes and desired columns
    n_rows = -(-len(genes_to_plot) // n_cols)  # Ceiling division

    # Set up the plot
    fig, axes = plt.subplots(n_rows, n_cols, figsize=fig_size)
    fig.suptitle(f"{title} Marker Genes", fontsize=12)

    # Flatten axes array for easier indexing
    axes = axes.flatten()

    for i, gene in enumerate(genes_to_plot):
        ax = axes[i]
        sc.pl.umap(
            adata, color=gene, ax=ax, show=False, title=gene, layer=layer, vmax=3
        )

    # Remove extra subplots
    for j in range(i + 1, len(axes)): # pylint: disable=undefined-loop-variable
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()

    # Print genes not found in the data
    missing_genes = set(genes) - set(genes_to_plot)
    if missing_genes:
        print(
            f"The following genes were not found in the data: {', '.join(missing_genes)}"
        )

    # Print which layer was used for coloring
    if layer:
        print(f"Used '{layer}' layer for gene expression values.")
    else:
        print("Used default expression values (usually from adata.X).")

    return fig, axes
