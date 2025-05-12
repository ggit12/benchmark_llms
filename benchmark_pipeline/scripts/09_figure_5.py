"""
This script contains the code to generate Figure 5 in the manuscript. 
This plot is the 10 celltypes that have the highest inter-rater agreement 
and lowest agreement with manual annotation.

Writes out the following files:

1. 'res/08_figure_5/confusion_matrix_for_cells_topleft_of_agreement.svg'
# 2. 'res/08_figure_5/sankey_topleft_of_agreement.svg'
3. 'res/08_figure_5/gene_module_scores_in_phagocytes.svg'
4. 'res/08_figure_5/macrophage_module_umap_in_phagocytes.svg'
5. 'res/08_figure_5/monocyte_module_umap_in_phagocytes.svg'
6. 'res/08_figure_5/dendritic_module_umap_in_phagocytes.svg'

"""
# pylint: disable=line-too-long
# pylint: disable=invalid-name

import sys
import os
import pickle

import anndict as adt
import scanpy as sc
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from anndict.utils.anndata_ import filter_gene_list

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    find_indices_closest_to_4_corners,
    customize_figure,
    customize_scatterplot,
    PROVIDERS,
)

# Use a non-GUI backend
matplotlib.use('Agg')

# Configure LLM backend
provider = os.environ['PROVIDER_FOR_POST_PROCESSING']
llm_config = PROVIDERS[provider].copy()
llm_config['model'] = os.environ['MODEL_FOR_POST_PROCESSING']

adt.configure_llm_backend(**llm_config)

# Read the results
adata = sc.read_h5ad("./res/04_postprocess_results/ts2_de_novo_llm_annotated.h5ad")

# Read the agreement_df and agreement_weights_df
agreement_df = pickle.load(open("./res/08_figure_4/agreement_df.pkl", "rb"))
agreement_weights_df = pickle.load(open("./res/08_figure_4/agreement_weights_df.pkl", "rb"))

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

cm_colorbar_df = agreement_df.join(agreement_weights_df.rename('agreement_weight'))
cm_colorbar_df.columns = ['Agreement with Manual', 'Inter-rater Agreement', 'Cell Type Abundance (% of Atlas)']

#get the celltypes with lowest agreement and make confusion matrix
celltypes_topleft = find_indices_closest_to_4_corners(agreement_df, n=10)['top_left']
celltypes_topleft = [(i,) for i in celltypes_topleft] # Make compatible with adt.build_adata_dict

#get adata with these celltypes
# same as in 07_figure_3_and_s2.py, we use the consistent_including_manual cell type column to hone in on the cell types, and directly use their manual annotation for further investigation
consistent_manual_cell_type_col = "consistent_including_manual_" + manual_cell_type_col
adata_topleft = adt.concatenate_adata_dict(adt.build_adata_dict(adata, strata_keys=[consistent_manual_cell_type_col], desired_strata=celltypes_topleft))

# Panel A
# ai_cell_type_col = 'consistent_including_manual_claude-3-5-sonnet-20240620_simplified_ai_cell_type'
# ai_cell_type_col = 'consistent_including_manual_gpt-4o_simplified_ai_cell_type'

# Set the AI cell type column to the best performing model by overall binary agreement (the first element in './res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl')
with open('./res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl', 'rb') as f:
    ai_cell_type_col = pickle.load(f)[0]

# Or to manually specific a specific model (included for flexibility)
# ai_cell_type_col = 'consistent_including_manual_' + os.environ['MODEL_FOR_DETAILED_ANALYSIS'] + '_simplified_ai_cell_type'

cm_topleft = adt.plot_confusion_matrix_from_adata(adata_topleft,
    true_label_key=consistent_manual_cell_type_col,
    predicted_label_key=ai_cell_type_col,
    diagonalize=True,
    figsize=(5,5),
    # true_label_color_df=cm_colorbar_df, # This arg is defined in an alternate version of this function (check orignal benchmarking notebook for def)
    true_ticklabels=True,
    predicted_ticklabels=True,
    annot=False)

cm_topleft.savefig('./res/09_figure_5/confusion_matrix_for_cells_topleft_of_agreement.svg', format='svg')

# Not in the figure
# sankey_topleft = adt.plot_sankey(adata_topleft, cols=[ai_cell_type_col, manual_cell_type_col,])
# adt.save_sankey(sankey_topleft, filename='./res/08_figure_5/sankey_topleft_of_agreement.svg')

#look at the cell types closest to the top left corner of the agreement plot, (5 of the largest of these by abundance)
adata_topleft_large = adt.sample_and_drop(adata_topleft, strata_keys=[consistent_manual_cell_type_col], n_largest_groups=5)
adata_dict = adt.build_adata_dict(adata_topleft_large, strata_keys=[consistent_manual_cell_type_col])


def plot_llm_markers_adata_dict(adata, adt_key=None): # pylint: disable=redefined-outer-name
    """
    Retrieve a set of marker genes and plot its score for a specific cell type in an AnnData object.
    """

    cell_type = adt_key[0]

    # Get up to the top 3 most frequent LLM labels for this cell type
    top_3_labels = adata.obs[ai_cell_type_col].value_counts().index[:3].tolist()

    #look up marker genes and calculate scores
    adt.cell_type_marker_gene_score(adata, cell_types=[cell_type] + top_3_labels, species='Human', list_length=None, score_name="_marker_gene_score")

    # Get the score column names
    score_cols = adt.get_adata_columns(adata, ends_with='_marker_gene_score')

    # Get the marker genes used for each score
    marker_genes = {col: adata.var.index[adata.var[col]].tolist() for col in score_cols}
    marker_genes_df = (
    pd.Series(marker_genes)
      .reset_index()
      .rename(columns={
          'index': 'score_col',
          0: 'genes_used'
      })
    )

    marker_genes_df.to_html(f'./res/09_figure_5/genes_used_in_scores_{cell_type}.html', index=False)

    # Plot the scores
    fig, axes = adt.module_score_umap(adata, score_cols=score_cols + [consistent_manual_cell_type_col] + [ai_cell_type_col])

    # Save a version of the plot with the legend
    fig.savefig(f'./res/09_figure_5/marker_gene_scores_withlegend_{cell_type}.svg', format='svg')

    for ax in np.ravel(axes):
        customize_scatterplot((fig, ax), fig_width=10.5, fig_height=7)

    # Save a version of the plot without the legend
    fig.savefig(f'./res/09_figure_5/marker_gene_scores_{cell_type}.svg', format='svg')

    #TODO: add bar plots


adata_dict.fapply(plot_llm_markers_adata_dict, use_multithreading=False, catch_errors=False)

#save a .done file for the marker gene score plots
with open('./res/09_figure_5/marker_gene_scores.done', 'w', encoding='utf-8') as f:
    f.write('done')


# # Get up to the top 3 most frequent LLM labels for this cell type
# top_3_labels = adata_single.obs[ai_cell_type_col].value_counts().index[:3].tolist()

# #look up marker genes and calculate scores
# adt.cell_type_marker_gene_score(adata, cell_types=[largest_cell_type] + top_3_labels, species='Human', list_length=None, score_name="_marker_gene_score")

# # Get the score column names
# score_cols = adt.get_adata_columns(adata, ends_with='_marker_gene_score')

# # Get the marker genes used for each score
# marker_genes = {col: adata.var.index[adata.var[col]].tolist() for col in score_cols}
# pd.DataFrame(marker_genes).to_html('./res/09_figure_5/genes_used_in_scores.html')

# # Plot the scores
# fig, axes = adt.module_score_umap(adata, score_cols=score_cols + [manual_cell_type_col] + [ai_cell_type_col])

# # Save a version of the plot with the legend
# fig.savefig('./res/09_figure_5/marker_gene_scores_withlegend.svg', format='svg')

# for ax in np.ravel(axes):
#     customize_scatterplot((fig, ax))

# # Save a version of the plot without the legend
# fig.savefig('./res/09_figure_5/marker_gene_scores.svg', format='svg')

# Old code for generating gene lists and plotting just looking at mononuclear phagocytes
# Code below here won't run unless on bigger object


# Looking into mononuclear phagocytes
# adt.display_html_summary(adt.summarize_metadata(adata_dict[('Stomach',)], cols=['donor*cell_ontology_class', 'donor']))

#Also run the manually curated mononuclear phagocyte analysis if present
# Check if phagocytes are present
phagocyte_present = 'Phagocyte' in adata_topleft.obs[consistent_manual_cell_type_col].unique() or \
                    'mononuclear phagocyte' in adata_topleft.obs[manual_cell_type_col].unique()

if phagocyte_present:
    # Filter the adata object to get phagocytes
    adata_mp = adata_topleft[adata_topleft.obs[manual_cell_type_col].isin(['mononuclear phagocyte'])]

    gene_lists = {
        "General_Mononuclear_Phagocyte": [
            "CD68",
            "CD14",
            "CSF1R",
            "ITGAM"  # CD11b
        ],

        "Macrophage_Specific": [
            "ADGRE1",  # EMR1, human equivalent of F4/80
            "CD163",
            "MERTK",
            "MSR1",   # CD204
            "MRC1"    # CD206
        ],

        "Stomach_Tissue_Resident_Macrophage": [
            "LYVE1",
            "TIMD4",
            "GATA6",
            "SIGLEC1"  # CD169, often on tissue-resident macrophages
        ],

        "Monocyte_Markers": [
            "FCN1",
            "S100A8",
            "S100A9",
            "CCR2"
        ],

        "Dendritic_Cell_Markers": [
            "ZBTB46",
            "FLT3",
            "CD1C",
            "CLEC10A"  # CD301, often on human DCs
        ],

        "Macrophage_Function": [
            "MARCO",
            "MMP9",
            "TIMP3",
            "FTL",
            "FTH1",
            "CD64"     # FCGR1A, important for phagocytosis
        ],

        "Stomach_Specific_Factors": [
            "SLC9A3",  # Sodium/hydrogen exchanger (acid resistance)
            "ATP4A",   # Proton pump (acid production)
            "TFF2",    # Trefoil factor 2 (mucosal protection)
            "MUC5AC",  # Mucin 5AC (mucosal protection)
            "CXCL8"    # IL-8, often produced by human stomach macrophages
        ],

        "Human_Specific_Macrophage_Markers": [
            "CD16",   # FCGR3A, expressed on some human macrophage subsets
            "CD32",   # FCGR2A, another Fc receptor on human macrophages
            "CD64",   # FCGR1A, high-affinity Fc receptor
            "CD11c",  # ITGAX, expressed on human macrophages and DCs
            "HLA-DRA" # MHC Class II, for antigen presentation
        ]
    }


    # Define the gene sets
    # macrophage_genes = ['CD68', 'CD163', 'MRC1', 'MERTK', 'CSF1R', 'EMR1']
    # monocyte_genes = ['CD14', 'CCR2', 'S100A8', 'S100A9']
    # dendritic_genes = ['ITGAX', 'ZBTB46', 'CLEC9A']

    macrophage_genes = gene_lists['Macrophage_Specific']
    monocyte_genes = gene_lists['Monocyte_Markers']  + ['CD14']
    dendritic_genes = gene_lists['Dendritic_Cell_Markers']

    # Filter gene lists to include only genes present in the dataset
    macrophage_genes_filtered = filter_gene_list(adata_mp, macrophage_genes)
    monocyte_genes_filtered = filter_gene_list(adata_mp, monocyte_genes)
    dendritic_genes_filtered = filter_gene_list(adata_mp, dendritic_genes)

    # Print filtered gene lists (optional)
    print('Macrophage genes used:', macrophage_genes_filtered)
    print('Monocyte genes used:', monocyte_genes_filtered)
    print('Dendritic genes used:', dendritic_genes_filtered)

    # Calculate module scores for each gene set
    sc.tl.score_genes(adata_mp, gene_list=macrophage_genes_filtered, score_name='macrophage_score')
    sc.tl.score_genes(adata_mp, gene_list=monocyte_genes_filtered, score_name='monocyte_score')
    sc.tl.score_genes(adata_mp, gene_list=dendritic_genes_filtered, score_name='dendritic_score')

    # Create a DataFrame with the module scores
    module_scores = adata_mp.obs[['macrophage_score', 'monocyte_score', 'dendritic_score']]

    # Calculate the mean module score for each gene set
    mean_scores = module_scores.mean()

    # Create figure and axes using plt.subplots()


    # Plot the mean module scores as a bar plot
    module_fig, module_ax = plt.subplots(figsize=(8, 6))
    mean_scores.plot(kind='bar', color=['skyblue', 'salmon', 'lightgreen'])
    plt.ylabel('Mean Module Score')
    plt.title('Mean Module Scores for Mononuclear Phagocyte Subtypes')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    #Panel B
    customize_figure((module_fig, module_ax), fig_width=2.4, fig_height=3, new_tick_labels={"macrophage_score":"Macrophage Score",
                                                                                            "monocyte_score":"Monocyte Score",
                                                                                            "dendritic_score":"Dendritic Cell Score"})

    module_fig.savefig('./res/09_figure_5/gene_module_scores_in_phagocytes.svg', format='svg')


    # Panel C
    # Plot and customize the first UMAP (Macrophage Score)
    fig1 = sc.pl.umap(adata_mp, color='macrophage_score', title='Macrophage Module Score', vmax='p99', return_fig=True)
    ax1 = fig1.axes[0]
    customize_scatterplot((fig1, ax1))
    fig1.savefig('./res/09_figure_5/macrophage_module_umap_in_phagocytes.svg', format='svg')

    # Plot and customize the second UMAP (Monocyte Score)
    fig2 = sc.pl.umap(adata_mp, color='monocyte_score', title='Monocyte Module Score', vmax='p99', return_fig=True)
    ax2 = fig2.axes[0]
    customize_scatterplot((fig2, ax2))
    fig2.savefig('./res/09_figure_5/monocyte_module_umap_in_phagocytes.svg', format='svg')

    # Plot and customize the third UMAP (Dendritic Cell Score)
    fig3 = sc.pl.umap(adata_mp, color='dendritic_score', title='Dendritic Cell Module Score', vmax='p99', return_fig=True)
    ax3 = fig3.axes[0]
    customize_scatterplot((fig3, ax3))
    fig3.savefig('./res/09_figure_5/dendritic_module_umap_in_phagocytes.svg', format='svg')


    #To confirm which UMAP is which score
    sc.pl.umap(adata_mp, color='macrophage_score', title='Macrophage Module Score', vmax='p99')
    sc.pl.umap(adata_mp, color='monocyte_score', title='Monocyte Module Score', vmax='p99')
    sc.pl.umap(adata_mp, color='dendritic_score', title='Dendritic Cell Module Score', vmax='p99')
