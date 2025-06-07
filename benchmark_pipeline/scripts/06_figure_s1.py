"""
Figure S1.
"""
# pylint: disable=line-too-long

import os
import sys
import pickle

import anndict as adt
import scanpy as sc
import matplotlib

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import customize_clustermap
from src import REMOVE_TICK_LABELS as remove_tick_labels

# Use a non-GUI backend
matplotlib.use('Agg')


# Read the results
adata = sc.read_h5ad("./res/04_postprocess_results/adt_de_novo_llm_annotated.h5ad")

# And the various column names
manual_cell_type_col = pickle.load(
    open("../../dat/manual_cell_type_col.pkl", "rb")
)
binary_agreement_cols = pickle.load(
    open("./res/04_postprocess_results/binary_agreement_cols.pkl", "rb")
)
categorical_agreement_cols = pickle.load(
    open("./res/04_postprocess_results/categorical_agreement_cols.pkl", "rb")
)
perfect_only_categorical_agreement_cols = pickle.load(
    open("./res/04_postprocess_results/perfect_only_categorical_agreement_cols.pkl", "rb")
)
binary_agreement_cols_top_models = pickle.load(
    open("./res/05_figure_2_and_table_2/binary_agreement_cols_top_models.pkl", "rb")
)


# Loop over the binary agreement columns and generate a plot for each
for col in binary_agreement_cols_top_models:
    agreement_plot = adt.plot_model_agreement(
        adata,
        group_by="consistent_including_manual_" + manual_cell_type_col,
        sub_group_by="tissue",
        agreement_cols=[col],  # Use the current binary agreement column
        granularity=2,
        legend=True,
    )

    # Write a version of the plot with the legend
    model_used = col.replace(
        "binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_",
        "",
    )

    agreement_plot.fig.savefig(
        f"res/06_figure_s1/agreement_plot_tissue_celltype_{model_used}_withlegend.svg",
        format="svg",
    )

    # Remove the tissue legend from the plot
    agreement_plot.ax_col_dendrogram.get_legend().remove() 

    agreement_plot_custom = customize_clustermap(
        agreement_plot,
        remove_legend=True,
        x_tick_substrings=[
            "binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_",
            "_simplified_ai_cell_type",
        ]
        + ["('" + i + "', '" for i in adata.obs["tissue"].unique().tolist()]
        + ["')"],
        new_ylabel="",
        fig_width=2.5,
        fig_height=11,
        new_tick_labels=remove_tick_labels,
    )

    # write the plot without the legend
    agreement_plot_custom[0].savefig(
        f"res/06_figure_s1/agreement_plot_tissue_celltype_{model_used}.svg",
        format="svg",
    )

# Write a done file to indicate that the script has completed
with open("res/06_figure_s1/done", "w", encoding="utf-8") as f:
    f.write("done")
