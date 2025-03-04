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


# from anndict.utils.anndata_ import filter_gene_list

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
)

from src import REMOVE_TICK_LABELS as remove_tick_labels

# Use a non-GUI backend
matplotlib.use('Agg')

# Read the results
adata = sc.read_h5ad("./res/03_gather_outputs/ts2_de_novo_llm_annotated.h5ad")

# And the various column names
manual_cell_type_col = pickle.load(
    open("./res/03_gather_outputs/manual_cell_type_col.pkl", "rb")
)

binary_agreement_cols_top_models = pickle.load(
    open("./res/04_figure_2_and_table_2/binary_agreement_cols_top_models.pkl", "rb")
)

perfect_only_categorical_agreement_cols_top_models = pickle.load(
    open(
        "./res/04_figure_2_and_table_2/perfect_only_categorical_agreement_cols_top_models.pkl",
        "rb",
    )
)

llm_celltype_cols_top_models = pickle.load(
    open("./res/03_gather_outputs/llm_celltype_cols.pkl", "rb")
)

# get the top 10 cell types by count
adata_large_celltypes = adt.wrappers.sample_and_drop_adata_dict(
    {"adata": adata}, strata_keys=["cell_ontology_class"], n_largest_groups=10
)["adata"]


# Panel A

# Top of A
# plot the binary agreement
agreement_plot_top_celltypes = adt.plot_model_agreement(
    adata_large_celltypes,
    group_by="consistent_including_manual_" + manual_cell_type_col,
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
    "./res/06_figure_3_and_s2/agreement_plot_largest_celltypes.svg", format="svg"
)

# Bottom of A
# plot the perfect match only agreement
agreement_plot_overall_categorical_perfect_celltypes = adt.plot_model_agreement(
    adata_large_celltypes,
    group_by="consistent_including_manual_" + manual_cell_type_col,
    sub_group_by="tissue",
    agreement_cols=perfect_only_categorical_agreement_cols_top_models,
    granularity=1,
)
customize_figure(
    agreement_plot_overall_categorical_perfect_celltypes,
    remove_legend=True,
    x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
    new_ylabel="Agreement with Manual Annotation (% perfect match)",
    fig_width=2.4,
    fig_height=3,
)
agreement_plot_overall_categorical_perfect_celltypes[0].savefig(
    "./res/06_figure_3_and_s2/agreement_plot_largest_celltypes_perfect_only.svg",
    format="svg",
)


# Figure S2
# Initialize an empty list to store the generated plots
agreement_plots_by_tissue_celltype_top_celltypes = []

# Loop over the binary agreement columns and generate a plot for each
for col in binary_agreement_cols_top_models:
    agreement_plot_temp = adt.plot_model_agreement(
        adata_large_celltypes,
        group_by="consistent_including_manual_" + manual_cell_type_col,
        sub_group_by="tissue",
        agreement_cols=[col],  # Use the current binary agreement column
        granularity=2,
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

    # Append the customized plot to the list
    agreement_plots_by_tissue_celltype_top_celltypes.append(agreement_plot_custom)

    # write the plot
    model_used = col.replace(
        "binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_",
        "",
    )
    agreement_plot_custom[0].savefig(
        f"./res/06_figure_3_and_s2/agreement_plot_tissue_celltype_{model_used}_top_celltypes.svg",
        format="svg",
    )

# Write a done file for the above agreement plots
with open(
    "./res/06_figure_3_and_s2/agreement_plots_by_tissue_celltype_top_celltypes_done",
    "w", encoding="utf-8",
) as f:
    f.write("done")


# To determine which cell types to look at in Panels B and D
# compute interrater agreement
agreement_df_large_celltypes = compute_agreement_df(
    adata_large_celltypes,
    rater_cols=llm_celltype_cols_top_models,
    manual_col="consistent_including_manual_" + manual_cell_type_col,
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
    "./res/06_figure_3_and_s2/agreement_scatterplot_largest_celltypes_top_llms.svg",
    format="svg",
)


# TODO: Code below here won't run unless on bigger object


# # get adata_dict of celltypes in top-left
# adata_top_left_cells = adt.build_adata_dict(
#     adata_large_celltypes,
#     strata_keys=["consistent_including_manual_cell_ontology_class"],
#     desired_strata=[("Basal Cell",), ("Stromal Cell",)],
# )

# # Panel B and D
# # make sankey plots among the top models and the manual annotation
# sankey_ai_to_man = adt.wrappers.plot_sankey_adata_dict(
#     adata_top_left_cells,
#     cols=["cell_ontology_class"]
#     + [
#         "consistent_including_manual_claude-3-5-sonnet-20240620_simplified_ai_cell_type"
#     ],
#     # params={'edge_color': "grey"}
# )
# adt.wrappers.save_sankey_adata_dict(
#     sankey_ai_to_man,
#     filename="./res/06_figure_3_and_s2/ai_to_manual_top_left_cells.svg",  # Dynamic filename for each plot
# )

# # Above might fail for basal cells, so do:
# sankey_ai_to_man_basal_cells = adt.plot_sankey(
#     adata_top_left_cells[("Basal Cell",)],
#     cols=["cell_ontology_class"]
#     + [
#         "consistent_including_manual_claude-3-5-sonnet-20240620_simplified_ai_cell_type"
#     ],
#     # params={'edge_color': "grey"}
# )
# adt.save_sankey(
#     sankey_ai_to_man_basal_cells,
#     filename="./res/06_figure_3_and_s2/ai_to_manual_top_left_cells_basal_cells.svg",  # Dynamic filename for each plot
# )

# # Get adata for basal and stromal cells
# # make a sankey plot for Basal Cells from manual to claude 3.5 sonnet (best performing)
# adata_basal_cells = adt.build_adata_dict(
#     adata_large_celltypes,
#     strata_keys=["consistent_including_manual_cell_ontology_class"],
#     desired_strata=[("Basal Cell",), ("Stromal Cell",)],
# )

# # # Panel C
# # test = adt.build_adata_dict(adata_top_left_cells['Basal Cell'].copy(), strata_keys=['tissue'])

# # adt.neighbors_adata_dict(test, n_neighbors=15)

# # # Step 2: Recalculate the UMAP embedding
# # adt.calculate_umap_adata_dict(test)

# test = adata_top_left_cells[("Basal Cell",)].copy()
# sc.pp.neighbors(test, n_neighbors=15)

# # Step 2: Recalculate the UMAP embedding
# sc.tl.umap(test)

# # Panel C

# # Defining gene sets for epithelial and basal cell markers
# epithelial_genes = ["CDH1", "EPCAM", "KRT8"]  # Epithelial Cell Markers
# basal_genes = ["KRT5", "KRT14", "TP63"]  # Basal Cell Markers
# # keratinocyte_genes = ['KRT1', 'KRT10', 'IVL'] #Keratinocyte Markers


# # # Function to filter genes present in the dataset
# # # pylint: disable=redefined-outer-name
# # def filter_genes(genes, adata):
# #     """Removes genes from ``genes`` that are not present in ``adata``."""
# #     return [gene for gene in genes if gene in adata.var_names]


# # Filter gene lists to include only genes present in the dataset for 'epithelial' and 'basal' cells
# epithelial_genes_filtered = adt.filter_gene_list(epithelial_genes, test)
# basal_genes_filtered = adt.filter_gene_list(basal_genes, test)
# # keratinocyte_genes_filtered = adt.filter_gene_list(keratinocyte_genes, test)

# # Print filtered gene lists (optional)
# print("Epithelial genes used:", epithelial_genes_filtered)
# print("Basal genes used:", basal_genes_filtered)
# # print('Keratinocyte genes used:', keratinocyte_genes_filtered)

# # Calculate module scores for each gene set
# sc.tl.score_genes(
#     test, gene_list=epithelial_genes_filtered, score_name="epithelial_score"
# )
# sc.tl.score_genes(test, gene_list=basal_genes_filtered, score_name="basal_score")
# # sc.tl.score_genes(test, gene_list=keratinocyte_genes_filtered, score_name='keratinocyte_score')

# # Remove outliers by capping scores at a reasonable quantile (e.g., 99th percentile)
# # def remove_outliers(data, score_name):
# #     upper_limit = np.percentile(data[score_name], 99)
# #     data[score_name] = np.clip(data[score_name], None, upper_limit)

# # remove_outliers(test.obs, 'epithelial_score')
# # remove_outliers(test.obs, 'basal_score')
# # remove_outliers(test.obs, 'keratinocyte_score')

# # Create a DataFrame with the module scores
# module_scores = test.obs[["epithelial_score", "basal_score"]]

# # Calculate the mean module score for each gene set
# mean_scores = module_scores.mean()

# # Create figure and axes using plt.subplots
# basal_module_fig, basal_module_ax = plt.subplots(figsize=(8, 6))

# # Plot the mean module scores as a bar plot
# mean_scores.plot(kind="bar", color=["skyblue", "salmon", "lightgreen"])
# plt.ylabel("Mean Module Score")
# plt.title("Mean Module Scores for Epithelial and Basal Cells")
# plt.xticks(rotation=45)
# plt.tight_layout()
# plt.show()

# # UMAP colored by epithelial module score
# epi_fig = sc.pl.umap(
#     test,
#     color="epithelial_score",
#     title="Epithelial Module Score",
#     vmax=1.5,
#     vmin=0,
#     return_fig=True,
# )
# epi_ax = epi_fig.axes[0]


# # UMAP colored by basal module score
# basal_fig = sc.pl.umap(
#     test,
#     color="basal_score",
#     title="Basal Module Score",
#     vmax=1.5,
#     vmin=0,
#     return_fig=True,
# )
# basal_ax = basal_fig.axes[0]

# # UMAP colored by basal module score
# # ker_fig = sc.pl.umap(test, color='keratinocyte_score', title='Keratinocyte Module Score', vmax=1.5, vmin=0, return_fig=True)

# customize_figure(
#     (basal_module_fig, basal_module_ax),
#     fig_width=1.6,
#     fig_height=3,
#     new_tick_labels={
#         "epithelial_score": "Epithelial Cell Score",
#         "basal_score": "Basal Cell Score",
#     },
# )
# basal_module_fig.savefig(
#     "./res/06_figure_3_and_s2/gene_module_scores_in_basal_cells.svg", format="svg"
# )

# customize_scatterplot((epi_fig, epi_ax))
# epi_fig.savefig(
#     "./res/06_figure_3_and_s2/epithelial_module_umap_in_basal_cells.svg", format="svg"
# )

# customize_scatterplot((basal_fig, basal_ax))
# basal_fig.savefig(
#     "./res/06_figure_3_and_s2/basal_module_umap_in_basal_cells.svg", format="svg"
# )


# # Panel E

# test_ov = adata_top_left_cells[("Stromal Cell",)].copy()

# sc.pp.neighbors(test_ov, n_neighbors=15)

# # Step 2: Recalculate the UMAP embedding
# sc.tl.umap(test_ov)

# # Importing necessary libraries

# # Defining gene sets for granulosa and stromal cells
# granulosa_genes = ["AMH", "HSD17B1", "SERPINE2", "GSTA1"]
# stromal_genes = ["DCN", "LUM"]


# # Filter gene lists to include only genes present in the dataset for 'stromal' and 'granulosa'
# granulosa_genes_filtered = adt.filter_gene_list(granulosa_genes, test_ov)
# stromal_genes_filtered = adt.filter_gene_list(stromal_genes, test_ov)

# # Print filtered gene lists (optional)
# print("Granulosa genes used:", granulosa_genes_filtered)
# print("Stromal genes used:", stromal_genes_filtered)

# # Calculate module scores for each gene set
# sc.tl.score_genes(
#     test_ov, gene_list=granulosa_genes_filtered, score_name="granulosa_score"
# )
# sc.tl.score_genes(test_ov, gene_list=stromal_genes_filtered, score_name="stromal_score")

# # Create a DataFrame with the module scores
# module_scores = test_ov.obs[["granulosa_score", "stromal_score"]]

# # Calculate the mean module score for each gene set
# mean_scores = module_scores.mean()

# # Create figure and axes using plt.subplots
# ov_module_fig, ov_module_ax = plt.subplots(figsize=(8, 6))

# # Plot the mean module scores as a bar plot
# mean_scores.plot(kind="bar", color=["skyblue", "salmon"])
# plt.ylabel("Mean Module Score")
# plt.title("Mean Module Scores for Granulosa and Stromal Cells")
# plt.xticks(rotation=45)
# plt.tight_layout()
# plt.show()

# # UMAP colored by epithelial module score
# gran_fig = sc.pl.umap(
#     test_ov,
#     color="granulosa_score",
#     title="Granulosa Module Score",
#     vmax=1.5,
#     vmin=0,
#     return_fig=True,
# )
# gran_ax = gran_fig.axes[0]


# # UMAP colored by basal module score
# stromal_fig = sc.pl.umap(
#     test_ov,
#     color="stromal_score",
#     title="Stromal Module Score",
#     vmax=1.5,
#     vmin=0,
#     return_fig=True,
# )
# stromal_ax = stromal_fig.axes[0]

# customize_figure(
#     (ov_module_fig, ov_module_ax),
#     fig_width=1.6,
#     fig_height=3,
#     new_tick_labels={
#         "granulosa_score": "Granulosa Cell Score",
#         "stromal_score": "Stromal Cell Score",
#     },
# )

# ov_module_ax.set_ylim(-0.05, 0.2)

# # Add a horizontal dashed line at y=1
# ov_module_ax.axhline(y=0, color="black", linestyle="--", linewidth=0.5)

# ov_module_fig.savefig(
#     "./res/06_figure_3_and_s2/gene_module_scores_in_stromal_cells.svg", format="svg"
# )

# customize_scatterplot((gran_fig, gran_ax))
# gran_fig.savefig(
#     "./res/06_figure_3_and_s2/granulosa_module_umap_in_stromal_cells.svg", format="svg"
# )

# customize_scatterplot((stromal_fig, stromal_ax))
# stromal_fig.savefig(
#     "./res/06_figure_3_and_s2/stromal_module_umap_in_stromal_cells.svg", format="svg"
# )
