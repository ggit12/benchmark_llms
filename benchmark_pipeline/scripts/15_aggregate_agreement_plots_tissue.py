"""
Aggregates agreement plots of the top models for the largest cell types from multiple runs.
"""
# pylint: disable=line-too-long

import argparse
import pickle
import os
import sys

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_figure,
    customize_legend_labels,
)

from src import MODEL_TICK_LABELS as model_tick_labels

#define a global variable to store sorted cell types for consistent ordering across plots
sorted_cell_types = None  # pylint: disable=invalid-name

def create_cell_type_model_plot(performance_df, perfect_match=False):
    """
    Create a grouped bar plot with cell types as major groups and models as subgroups.

    Parameters:
    -----------
    performance_df : DataFrame
        Agreement data for plotting with cell types as index.
    perfect_match : bool
        If True, use the perfect match data, otherwise use overall binary data.

    Returns:
    --------
    tuple
        (fig, ax) for the plot.
    """
    # Reset index to make cell types a column
    plot_df = performance_df.reset_index()

    # Column prefix based on perfect_match parameter
    column_prefix = "Perfect Match (% of Cells)_perfect_only_categorical_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_" \
        if perfect_match else \
            "Overall Binary (% of Cells)_binary_agreement_consistent_including_manual_cell_ontology_class_consistent_including_manual_"
    column_suffix = "_simplified_ai_cell_type"

    # Get the model names from the columns
    model_columns = [col for col in performance_df.columns if col.startswith(column_prefix)]
    model_names = [col.replace(column_prefix, "").replace(column_suffix, "") for col in model_columns]

    # Create a DataFrame for plotting
    plot_data = []
    for _, row in plot_df.iterrows():
        cell_type = row['Cell Type']  # The index was named 'Model' in the original DataFrame

        for col, model_name in zip(model_columns, model_names):
            plot_data.append({
                'cell_type': cell_type,
                'model_name': model_name,
                'agreement': row[col]
            })

    plot_df = pd.DataFrame(plot_data)

    # Sort cell types based on their mean agreement value across all models
    # Keep the order the same as the first plot for ease of comparison on plots
    first_plot = 'sorted_cell_types' not in globals()  # Check if this is the first plot
    if first_plot:
        global sorted_cell_types # pylint: disable=global-statement
        sorted_cell_types = plot_df.groupby('cell_type')['agreement'].mean().sort_values(ascending=False).index.tolist()
    else:
        pass


    # Convert cell_type to categorical with specific order
    plot_df['cell_type'] = pd.Categorical(
        plot_df['cell_type'], categories=sorted_cell_types, ordered=True # pylint: disable=possibly-used-before-assignment
    )

    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 8))

    # Set the color palette
    palette = plt.get_cmap('Paired')(np.linspace(0, 1, plot_df['model_name'].nunique()))

    # Plot the aggregated bars
    sns.barplot(
        data=plot_df,
        x='cell_type',
        y='agreement',
        hue='model_name',
        palette=palette,
        ax=ax,
        order=sorted_cell_types,
        errorbar='sd',
        capsize=0.2,
        err_kws={"linewidth": 0.5, "color": "black"},
    )

    # Add individual data points
    sns.stripplot(
        data=plot_df,
        x='cell_type',
        y='agreement',
        hue='model_name',
        dodge=True,
        ax=ax,
        order=sorted_cell_types,
        alpha=1.0,
        size=1.5,
        color="black",
        legend=False,
    )

    # # Add proportion labels on top of each bar
    # for p in ax.patches:
    #     height = p.get_height()
    #     if height > 0:
    #         ax.text(
    #             p.get_x() + p.get_width() / 2.0,
    #             height + 0.01,
    #             f"{height * 100:.0f}%",
    #             ha="center",
    #         )

    # Rotate x-axis tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # Set labels and title
    plot_type = "Perfect Match" if perfect_match else "Binary"
    ax.set_xlabel("Cell Type")
    ax.set_ylabel(f"{plot_type}\n(Proportion)")
    ax.set_title("Agreemet\n(with manual annotation)")
    ax.set_ylim(0, 1.05)
    ax.legend(title="Models", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    return fig, ax


def main():
    parser = argparse.ArgumentParser(
        description="Create agreement plots for largest cell types."
    )
    parser.add_argument(
        "--input-tables",
        nargs="+",
        required=True,
        help="Input agreement table pickle files for largest cell types",
    )
    parser.add_argument(
        "--input-tables-perfect-only",
        nargs="+",
        required=True,
        help="Input agreement table pickle files for perfect match only",
    )

    args = parser.parse_args()

    print(f"Input tables: {args.input_tables}")
    # print(f"Input tables perfect only: {args.input_tables_perfect_only}")

    # Load agreement tables
    agreement_tables = [pickle.load(open(f, 'rb')) for f in args.input_tables]
    print(f"Loaded agreement tables[0]:\n{agreement_tables[0]}")
    agreement_tables_perfect = [pickle.load(open(f, 'rb')) for f in args.input_tables_perfect_only]

    # Reset index for each table for concatenaton
    [table.reset_index(inplace=True) for table in agreement_tables] # pylint: disable=expression-not-assigned
    print(f"Reset index for agreement tables[0]:\n{agreement_tables[0]}")
    [table.reset_index(inplace=True) for table in agreement_tables_perfect] # pylint: disable=expression-not-assigned

    # Concatenate tables from different runs
    agreement_tables = pd.concat(agreement_tables)
    print(f"Concatenated agreement tables:\n{agreement_tables}")
    agreement_tables_perfect = pd.concat(agreement_tables_perfect)

    agreement_tables.reset_index(drop=True, inplace=True)
    agreement_tables_perfect.reset_index(drop=True, inplace=True)


    # Create overall binary agreement plot
    # Create agreement plot for largest cell types
    agreement_plot = create_cell_type_model_plot(
        agreement_tables, perfect_match=False
    )

    # Fix the legend labels as well
    customize_legend_labels(
        agreement_plot,
        label_map=model_tick_labels,
    )

    # Fix legend position
    agreement_plot[1].legend(title="Models", bbox_to_anchor=(1.05, 1), loc="upper right")

    # Save a version of the plot with the legend
    agreement_plot[0].savefig(
        "res/15_aggregated_plots_tissue/agreement_plot_by_tissue_withlegend.svg", format="svg"
    )

    # Customize figure and remove the legend
    customize_figure(
        agreement_plot,
        remove_legend=True,
        x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
        fig_width=4.8,
        fig_height=3,
    )

    # Save a version of the plot with the legend
    agreement_plot[0].savefig(
        "res/15_aggregated_plots_tissue/agreement_plot_by_tissue.svg", format="svg"
    )


    # Create perfect match agreement plot
    agreement_plot_perfect_only = create_cell_type_model_plot(
        agreement_tables_perfect, perfect_match=True
    )

    # Fix the legend labels as well
    customize_legend_labels(
        agreement_plot_perfect_only,
        label_map=model_tick_labels,
    )

    # Fix legend position
    agreement_plot_perfect_only[1].legend(title="Models", bbox_to_anchor=(1.05, 1), loc="upper right")

    # Save a version of the plot with the legend
    agreement_plot_perfect_only[0].savefig(
        "res/15_aggregated_plots_tissue/agreement_plot_by_tissue_perfect_only_withlegend.svg", format="svg"
    )

    # Customize figure and remove the legend
    customize_figure(
        agreement_plot_perfect_only,
        remove_legend=True,
        x_tick_substrings=["agreement_cell_ontology_class_", "_simplified_ai_cell_type"],
        fig_width=4.8,
        fig_height=3,
    )

    # Save a version of the plot without the legend
    agreement_plot_perfect_only[0].savefig(
        "res/15_aggregated_plots_tissue/agreement_plot_by_tissue_perfect_only.svg", format="svg"
    )



if __name__ == "__main__":
    main()
