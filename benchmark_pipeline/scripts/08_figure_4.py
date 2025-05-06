"""
Code for plotting Figure 4, the scatterplots of inter-rater agreement vs manual agreement.
"""
# pylint: disable=line-too-long

import sys
import os
import pickle

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_scatterplot,
    compute_agreement_df,
    plot_agreement,
    calculate_weights,
)

# Use a non-GUI backend
matplotlib.use('Agg')

# Read the results
adata = sc.read_h5ad("./res/04_postprocess_results/ts2_de_novo_llm_annotated.h5ad")

# And the various column names
manual_cell_type_col = pickle.load(
    open("../../dat/manual_cell_type_col.pkl", "rb")
)

binary_agreement_cols_top_models = pickle.load(
    open("./res/05_figure_2_and_table_2/binary_agreement_cols_top_models.pkl", "rb")
)

categorical_agreement_cols_top_models = pickle.load(
    open(
        "./res/05_figure_2_and_table_2/categorical_agreement_cols_top_models.pkl",
        "rb",
    )
)

perfect_only_categorical_agreement_cols_top_models = pickle.load(
    open(
        "./res/05_figure_2_and_table_2/perfect_only_categorical_agreement_cols_top_models.pkl",
        "rb",
    )
)

llm_celltype_cols_top_models = pickle.load(
    open("./res/04_postprocess_results/llm_celltype_cols.pkl", "rb")
)

ten_largest_cell_types = pickle.load(
    open('./res/07_figure_3_and_s2/ten_largest_cell_types.pkl', 'rb')
    )

#compute interrater agreement across whole atlas
agreement_df = compute_agreement_df(adata, rater_cols=llm_celltype_cols_top_models, manual_col='consistent_including_manual_' + manual_cell_type_col, agreement_type='plurality', normalize_values=False)


agreement_weights_df = calculate_weights(adata, 'consistent_including_manual_' + manual_cell_type_col, equal_weights=False)

plt.rcParams['svg.fonttype'] = 'none'


# Dots unweighted, marginal densities unormalized (panel A)
agreement_scatterplot_unweighted_unnormalized = plot_agreement(agreement_df=agreement_df, weights=agreement_weights_df, show_labels=False, show_legend=False, normalize_kdes=False)

customize_scatterplot(agreement_scatterplot_unweighted_unnormalized, remove_legend=True, new_ylabel='Inter-rater Agreement', new_xlabel='Agreement with Manual Annotation')
agreement_scatterplot_unweighted_unnormalized[0].savefig('./res/08_figure_4/agreement_scatterplot_overall_unweighted_unnormalized.svg', format='svg')


# Dot size wieghted by cell type size, marginal densities scaled to cell type sizes (panel B)
agreement_scatterplot_weighted_normalized = plot_agreement(agreement_df=agreement_df, weights=agreement_weights_df, show_labels=False, show_legend=False, normalize_kdes=True)

customize_scatterplot(agreement_scatterplot_weighted_normalized, remove_legend=True, new_ylabel='Inter-rater Agreement', new_xlabel='Agreement with Manual Annotation')
agreement_scatterplot_weighted_normalized[0].savefig('./res/08_figure_4/agreement_scatterplot_overall_weighted_normalized.svg', format='svg')

# Write agreement_df and agreement_weights_df
pickle.dump(agreement_df, open('./res/08_figure_4/agreement_df.pkl', 'wb'))
pickle.dump(agreement_weights_df, open('./res/08_figure_4/agreement_weights_df.pkl', 'wb'))


# Recompute agreement_df and agreement_weights_df without the ten largest cell types (used in downstream analysis 09_figure_5.py)
# This is because the largest cell types are analyzed separately
adata = adata[~adata.obs[manual_cell_type_col].isin(ten_largest_cell_types)].copy()

agreement_df_without_largest = compute_agreement_df(adata, rater_cols=llm_celltype_cols_top_models, manual_col='consistent_including_manual_' + manual_cell_type_col, agreement_type='plurality', normalize_values=False)
agreement_weights_df_without_largest = calculate_weights(adata, 'consistent_including_manual_' + manual_cell_type_col, equal_weights=False)

pickle.dump(agreement_df_without_largest, open('./res/08_figure_4/agreement_df_without_largest.pkl', 'wb'))
pickle.dump(agreement_weights_df_without_largest, open('./res/08_figure_4/agreement_weights_df_without_largest.pkl', 'wb'))
