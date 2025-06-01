"""
Code to make figure 3, the investigation on the top 10 largest cell types
"""

# pylint: disable=line-too-long

import sys
import os
import pickle

import anndict as adt
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt

from anndict.utils.anndata_ import filter_gene_list

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_clustermap,
    customize_figure,
    customize_scatterplot,
    compute_agreement_df,
    find_indices_closest_to_4_corners,
    plot_agreement_simple,
    extract_table_from_fig,
)

from src import REMOVE_TICK_LABELS as remove_tick_labels

# Use a non-GUI backend
matplotlib.use('Agg')

# Read the results
adata = sc.read_h5ad("./res/04_postprocess_results/adt_de_novo_llm_annotated.h5ad")

# And the various column names
manual_cell_type_col = pickle.load(
    open("../../dat/manual_cell_type_col.pkl", "rb")
)

binary_agreement_cols_top_models = pickle.load(
    open("./res/05_figure_2_and_table_2/binary_agreement_cols_top_models.pkl", "rb")
)

perfect_only_categorical_agreement_cols_top_models = pickle.load(
    open(
        "./res/05_figure_2_and_table_2/perfect_only_categorical_agreement_cols_top_models.pkl",
        "rb",
    )
)

llm_celltype_cols_top_models = pickle.load(
    open("./res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl", "rb")
)

# get the top 10 cell types by count
consistent_manual_cell_type_col = "consistent_including_manual_" + manual_cell_type_col
adata_large_celltypes = adt.sample_and_drop(adata, strata_keys=[manual_cell_type_col], n_largest_groups=10)

# Panel A

# Top of A
# plot the binary agreement
agreement_plot_top_celltypes = adt.plot_model_agreement(
    adata_large_celltypes,
    group_by=consistent_manual_cell_type_col,
    sub_group_by="tissue",
    agreement_cols=binary_agreement_cols_top_models,
    granularity=1,
)
customize_figure(
    agreement_plot_top_celltypes,
    remove_legend=True,
    x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
    new_ylabel="Agreement with Manual Annotation (yes/no)",
    fig_width=2.4,
    fig_height=3,
)

# Save the plot as an SVG file
agreement_plot_top_celltypes[0].savefig(
    "./res/07_figure_3_and_s2/agreement_plot_largest_celltypes.svg", format="svg"
)

# Extract values from the plot
agreement_table_top_celltypes = extract_table_from_fig(agreement_plot_top_celltypes, value_col_name="Overall Binary (% of Cells)", x_tick_label_name="Cell Type")

# Write for later aggregation
agreement_table_top_celltypes.to_pickle(
    "./res/07_figure_3_and_s2/agreement_table_largest_celltypes.pkl"
)

# Bottom of A
# plot the perfect match only agreement
agreement_plot_overall_categorical_perfect_top_celltypes = adt.plot_model_agreement(
    adata_large_celltypes,
    group_by=consistent_manual_cell_type_col,
    sub_group_by="tissue",
    agreement_cols=perfect_only_categorical_agreement_cols_top_models,
    granularity=1,
)
customize_figure(
    agreement_plot_overall_categorical_perfect_top_celltypes,
    remove_legend=True,
    x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
    new_ylabel="Agreement with Manual Annotation (% perfect match)",
    fig_width=2.4,
    fig_height=3,
)
agreement_plot_overall_categorical_perfect_top_celltypes[0].savefig(
    "./res/07_figure_3_and_s2/agreement_plot_largest_celltypes_perfect_only.svg",
    format="svg",
)

# Extract values from the plot
agreement_table_overall_categorical_perfect_top_celltypes = extract_table_from_fig(
    agreement_plot_overall_categorical_perfect_top_celltypes,
    value_col_name="Perfect Match (% of Cells)",
    x_tick_label_name="Cell Type",)

# Write for later aggregation
agreement_table_overall_categorical_perfect_top_celltypes.to_pickle(
    "./res/07_figure_3_and_s2/agreement_table_largest_celltypes_perfect_only.pkl"
)

# Figure S2

# Loop over the binary agreement columns and generate a plot for each
for col in binary_agreement_cols_top_models:
    agreement_plot_temp = adt.plot_model_agreement(
        adata_large_celltypes,
        group_by=consistent_manual_cell_type_col,
        sub_group_by="tissue",
        agreement_cols=[col],  # Use the current binary agreement column
        granularity=2,
    )

    # Write a version of the plot with the legend
    model_used = col.replace(
        "binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_",
        "",
    )

    agreement_plot_temp[0].savefig(
        f"./res/07_figure_3_and_s2/agreement_plot_tissue_celltype_{model_used}_withlegend.svg",
        format="svg",
    )

    agreement_plot_custom = customize_clustermap(
        agreement_plot_temp,
        remove_legend=True,
        x_tick_substrings=[
            "binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_",
            "_simplified_ai_cell_type",
        ]
        + [
            "('" + i + "', '"
            for i in adata_large_celltypes.obs["tissue"].unique().tolist()
        ]
        + ["')"],
        new_ylabel="",
        fig_width=2,
        fig_height=3,
        new_tick_labels=remove_tick_labels,
    )

    # write the plot
    agreement_plot_custom[0].savefig(
        f"./res/07_figure_3_and_s2/agreement_plot_tissue_celltype_{model_used}_top_celltypes.svg",
        format="svg",
    )

# Write a done file for the above agreement plots
with open(
    "./res/07_figure_3_and_s2/agreement_plots_by_tissue_celltype_top_celltypes_done",
    "w", encoding="utf-8",
) as f:
    f.write("done")


# To determine which cell types to look at in Panels B and D
# compute interrater agreement
# Consistency itself is calculated on the processed labeled columns that have been unified (i.e. made to share a common set of labels)
# For further investigation into cell types, we use the original labels themselves (after honing in on which labels to investigate based on the agreement scores)
agreement_df_large_celltypes = compute_agreement_df(
    adata_large_celltypes,
    rater_cols=llm_celltype_cols_top_models,
    manual_col=consistent_manual_cell_type_col,
    agreement_type="plurality",
    normalize_values=False,
)
celltypes_of_interest_large_celltypes = find_indices_closest_to_4_corners(
    agreement_df_large_celltypes, n=3
)


agreement_scatterplot_large_celltypes = plot_agreement_simple(
    agreement_df=agreement_df_large_celltypes
)

customize_scatterplot(
    agreement_scatterplot_large_celltypes,
    remove_legend=True,
    new_ylabel="Inter-rater Agreement",
    new_xlabel="Agreement with Manual Annotation",
)

# Save the plot as an SVG file
agreement_scatterplot_large_celltypes[0].savefig(
    "./res/07_figure_3_and_s2/agreement_scatterplot_largest_celltypes_top_llms.svg",
    format="svg",
)


# TODO: Code below here won't run unless on full object
# record the cell types that are closest to the top-left corner as a txt
with open(
    "./res/07_figure_3_and_s2/celltypes_of_interest_large_celltypes.txt",
    "w", encoding="utf-8",
) as f:
    f.write(
        "Large cell types of interest in the top left corner of the agreement plot:\n"
    )
    for celltype in celltypes_of_interest_large_celltypes["top_left"]:
        f.write(f"{celltype}\n")


# get adata_dict of celltypes closest to top-left (of the largest cell types)
adata_top_left_cells = adt.build_adata_dict(
    adata_large_celltypes,
    strata_keys=[manual_cell_type_col],
    desired_strata=[("basal cell",), ("stromal cell of ovary",)],
)

# Panel B and D

# Set the AI cell type column to the best performing model by overall binary agreement (the first element in './res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl')
with open('./res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl', 'rb') as f:
    ai_cell_type_col = pickle.load(f)[0]

# Or to manually specify a model (included for flexibility)
# ai_cell_type_col = 'consistent_including_manual_' + os.environ['MODEL_FOR_DETAILED_ANALYSIS'] + '_simplified_ai_cell_type'


# make sankey plots among the top models and the manual annotation
sankey_ai_to_man = adt.wrappers.plot_sankey_adata_dict(
    adata_top_left_cells,
    cols=[manual_cell_type_col]
    + [
        ai_cell_type_col
    ],
    # params={'edge_color': "grey"}
)
adt.wrappers.save_sankey_adata_dict(
    sankey_ai_to_man,
    filename="./res/07_figure_3_and_s2/ai_to_manual_top_left_cells.svg",  # Dynamic filename for each plot
)

# Write a done file for the above sankey plots
with open(
    "./res/07_figure_3_and_s2/ai_to_manual_top_left_cells_done",
    "w", encoding="utf-8",
) as f:
    f.write("done")


# Above might fail for basal cells, so do:
sankey_ai_to_man_basal_cells = adt.plot_sankey(
    adata_top_left_cells[("basal cell",)],
    cols=[manual_cell_type_col]
    + [
        ai_cell_type_col
    ],
    # params={'edge_color': "grey"}
)
adt.save_sankey(
    sankey_ai_to_man_basal_cells,
    filename="./res/07_figure_3_and_s2/ai_to_manual_top_left_cells_basal_cells.svg",  # Dynamic filename for each plot
)


# Panel C


adata_basal = adata_top_left_cells[("basal cell",)].copy()
sc.pp.neighbors(adata_basal, n_neighbors=15)

# Step 2: Recalculate the UMAP embedding
sc.tl.umap(adata_basal)

# Panel C

# Defining gene sets for epithelial and basal cell markers
epithelial_genes = ["CDH1", "EPCAM", "KRT8"]  # Epithelial Cell Markers
basal_genes = ["KRT5", "KRT14", "TP63"]  # Basal Cell Markers
# keratinocyte_genes = ['KRT1', 'KRT10', 'IVL'] #Keratinocyte Markers


# # Function to filter genes present in the dataset
# # pylint: disable=redefined-outer-name
# def filter_genes(genes, adata):
#     """Removes genes from ``genes`` that are not present in ``adata``."""
#     return [gene for gene in genes if gene in adata.var_names]


# Filter gene lists to include only genes present in the dataset for 'epithelial' and 'basal' cells
epithelial_genes_filtered = filter_gene_list(adata_basal, epithelial_genes,)
basal_genes_filtered = filter_gene_list(adata_basal, basal_genes)
# keratinocyte_genes_filtered = adt.filter_gene_list(keratinocyte_genes, adata_basal)

# Print filtered gene lists (optional)
print("Epithelial genes used:", epithelial_genes_filtered)
print("Basal genes used:", basal_genes_filtered)
# print('Keratinocyte genes used:', keratinocyte_genes_filtered)

# Calculate module scores for each gene set
sc.tl.score_genes(
    adata_basal, gene_list=epithelial_genes_filtered, score_name="epithelial_score"
)
sc.tl.score_genes(adata_basal, gene_list=basal_genes_filtered, score_name="basal_score")
# sc.tl.score_genes(adata_basal, gene_list=keratinocyte_genes_filtered, score_name='keratinocyte_score')

# Remove outliers by capping scores at a reasonable quantile (e.g., 99th percentile)
# def remove_outliers(data, score_name):
#     upper_limit = np.percentile(data[score_name], 99)
#     data[score_name] = np.clip(data[score_name], None, upper_limit)

# remove_outliers(adata_basal.obs, 'epithelial_score')
# remove_outliers(adata_basal.obs, 'basal_score')
# remove_outliers(adata_basal.obs, 'keratinocyte_score')

# Create a DataFrame with the module scores
module_scores = adata_basal.obs[["epithelial_score", "basal_score"]]

# Calculate the mean module score for each gene set
mean_scores = module_scores.mean()

# Create figure and axes using plt.subplots
basal_module_fig, basal_module_ax = plt.subplots(figsize=(8, 6))

# Plot the mean module scores as a bar plot
mean_scores.plot(kind="bar", color=["skyblue", "salmon", "lightgreen"])
plt.ylabel("Mean Module Score")
plt.title("Mean Module Scores for Epithelial and Basal Cells")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# UMAP colored by epithelial module score
epi_fig = sc.pl.umap(
    adata_basal,
    color="epithelial_score",
    title="Epithelial Module Score",
    vmax=1.5,
    vmin=0,
    return_fig=True,
)
epi_ax = epi_fig.axes[0]


# UMAP colored by basal module score
basal_fig = sc.pl.umap(
    adata_basal,
    color="basal_score",
    title="Basal Module Score",
    vmax=1.5,
    vmin=0,
    return_fig=True,
)
basal_ax = basal_fig.axes[0]

# UMAP colored by basal module score
# ker_fig = sc.pl.umap(adata_basal, color='keratinocyte_score', title='Keratinocyte Module Score', vmax=1.5, vmin=0, return_fig=True)

customize_figure(
    (basal_module_fig, basal_module_ax),
    fig_width=1.6,
    fig_height=3,
    new_tick_labels={
        "epithelial_score": "Epithelial Cell Score",
        "basal_score": "Basal Cell Score",
    },
)
basal_module_fig.savefig(
    "./res/07_figure_3_and_s2/gene_module_scores_in_basal_cells.svg", format="svg"
)

customize_scatterplot((epi_fig, epi_ax))
epi_fig.savefig(
    "./res/07_figure_3_and_s2/epithelial_module_umap_in_basal_cells.svg", format="svg"
)

customize_scatterplot((basal_fig, basal_ax))
basal_fig.savefig(
    "./res/07_figure_3_and_s2/basal_module_umap_in_basal_cells.svg", format="svg"
)


# Panel E

adata_stromal_ov = adata_top_left_cells[("stromal cell of ovary",)].copy()

sc.pp.neighbors(adata_stromal_ov, n_neighbors=15)

# Step 2: Recalculate the UMAP embedding
sc.tl.umap(adata_stromal_ov)

# Importing necessary libraries

# Defining gene sets for granulosa and stromal cells
granulosa_genes = ["AMH", "HSD17B1", "SERPINE2", "GSTA1"]
stromal_genes = ["DCN", "LUM"]


# Filter gene lists to include only genes present in the dataset for 'stromal' and 'granulosa'
granulosa_genes_filtered = filter_gene_list(adata_stromal_ov, granulosa_genes,)
stromal_genes_filtered = filter_gene_list(adata_stromal_ov, stromal_genes)

# Print filtered gene lists (optional)
print("Granulosa genes used:", granulosa_genes_filtered)
print("Stromal genes used:", stromal_genes_filtered)

# Calculate module scores for each gene set
sc.tl.score_genes(
    adata_stromal_ov, gene_list=granulosa_genes_filtered, score_name="granulosa_score"
)
sc.tl.score_genes(adata_stromal_ov, gene_list=stromal_genes_filtered, score_name="stromal_score")

# Create a DataFrame with the module scores
module_scores = adata_stromal_ov.obs[["granulosa_score", "stromal_score"]]

# Calculate the mean module score for each gene set
mean_scores = module_scores.mean()

# Create figure and axes using plt.subplots
ov_module_fig, ov_module_ax = plt.subplots(figsize=(8, 6))

# Plot the mean module scores as a bar plot
mean_scores.plot(kind="bar", color=["skyblue", "salmon"])
plt.ylabel("Mean Module Score")
plt.title("Mean Module Scores for Granulosa and Stromal Cells")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# UMAP colored by epithelial module score
gran_fig = sc.pl.umap(
    adata_stromal_ov,
    color="granulosa_score",
    title="Granulosa Module Score",
    vmax=1.5,
    vmin=0,
    return_fig=True,
)
gran_ax = gran_fig.axes[0]


# UMAP colored by basal module score
stromal_fig = sc.pl.umap(
    adata_stromal_ov,
    color="stromal_score",
    title="Stromal Module Score",
    vmax=1.5,
    vmin=0,
    return_fig=True,
)
stromal_ax = stromal_fig.axes[0]

customize_figure(
    (ov_module_fig, ov_module_ax),
    fig_width=1.6,
    fig_height=3,
    new_tick_labels={
        "granulosa_score": "Granulosa Cell Score",
        "stromal_score": "Stromal Cell Score",
    },
)

ov_module_ax.set_ylim(-0.05, 0.2)

# Add a horizontal dashed line at y=1
ov_module_ax.axhline(y=0, color="black", linestyle="--", linewidth=0.5)

ov_module_fig.savefig(
    "./res/07_figure_3_and_s2/gene_module_scores_in_stromal_cells.svg", format="svg"
)

customize_scatterplot((gran_fig, gran_ax))
gran_fig.savefig(
    "./res/07_figure_3_and_s2/granulosa_module_umap_in_stromal_cells.svg", format="svg"
)

customize_scatterplot((stromal_fig, stromal_ax))
stromal_fig.savefig(
    "./res/07_figure_3_and_s2/stromal_module_umap_in_stromal_cells.svg", format="svg"
)
