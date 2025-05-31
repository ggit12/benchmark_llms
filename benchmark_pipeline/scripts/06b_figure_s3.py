"""
Figure S3.
"""
# pylint: disable=line-too-long

import sys
import os
import pickle

import anndict as adt
import scanpy as sc
import matplotlib


from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_figure,
    extract_table_from_fig,
)

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

# # Set the AI cell type column to the best performing model by overall binary agreement (the first element in './res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl')
ai_cell_type_col = llm_celltype_cols_top_models[0]

# Make bar plots by tissue as well

# plot the binary agreement
agreement_plot_by_tissue = adt.plot_model_agreement(
    adata,
    group_by="tissue",
    sub_group_by="consistent_including_manual_" + manual_cell_type_col, #this arg is ignored when granularity=1
    agreement_cols=binary_agreement_cols_top_models,
    granularity=1,
)
customize_figure(
    agreement_plot_by_tissue,
    remove_legend=True,
    x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
    new_ylabel="Agreement with Manual Annotation (yes/no)",
    fig_width=2.4,
    fig_height=3,
)

# Save the plot as an SVG file
agreement_plot_by_tissue[0].savefig(
    "./res/06b_figure_s3/agreement_plot_by_tissue.svg", format="svg"
)

# Extract values from the plot
agreement_table_by_tissue = extract_table_from_fig(agreement_plot_by_tissue, value_col_name="Overall Binary (% of Cells)", x_tick_label_name="Tissue")

# Write for later aggregation
agreement_table_by_tissue.to_pickle(
    "./res/06b_figure_s3/agreement_table_by_tissue.pkl"
)

# Bottom of A
# plot the perfect match only agreement
agreement_plot_by_tissue_categorical_perfect_only = adt.plot_model_agreement(
    adata,
    group_by="tissue",
    sub_group_by="consistent_including_manual_" + manual_cell_type_col,  # this arg is ignored when granularity=1
    agreement_cols=perfect_only_categorical_agreement_cols_top_models,
    granularity=1,
)
customize_figure(
    agreement_plot_by_tissue_categorical_perfect_only,
    remove_legend=True,
    x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
    new_ylabel="Agreement with Manual Annotation (% perfect match)",
    fig_width=2.4,
    fig_height=3,
)
agreement_plot_by_tissue_categorical_perfect_only[0].savefig(
    "./res/06b_figure_s3/agreement_plot_by_tissue_perfect_only.svg",
    format="svg",
)

# Extract values from the plot
agreement_table_by_tissue_categorical_perfect_only = extract_table_from_fig(
    agreement_plot_by_tissue_categorical_perfect_only,
    value_col_name="Perfect Match (% of Cells)",
    x_tick_label_name="Tissue",)

# Write for later aggregation
agreement_table_by_tissue_categorical_perfect_only.to_pickle(
    "./res/06b_figure_s3/agreement_table_by_tissue_perfect_only.pkl"
)

# For the tissues with lowest agreement, how were the cells annotated?
sankey_ai_to_man_basal_cells = adt.plot_sankey(
    adata[("Ear",)],
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