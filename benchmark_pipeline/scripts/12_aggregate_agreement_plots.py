#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pickle
import os
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_figure,
)

from src import MODEL_TICK_LABELS as model_tick_labels


def create_agreement_plot(performance_tables, by_cell_types=False):
    """
    Create a grouped bar plot of agreement categories with individual data points.

    Parameters:
    -----------
    plot_df : DataFrame
        Aggregated data for plotting.
    individual_df : DataFrame
        Individual data points for scatter overlay.
    by_cell_types : bool
        If True, title includes "Cell Types", otherwise "Cells".

    Returns:
    --------
    tuple
        (fig, ax) for the plot.
    """
    column_type = "Cell Types" if by_cell_types else "Cells"
    column_prefix = f"Categorical Agreement (% of {column_type})"

    # Filter columns for categorical agreement
    agreement_cols = [col for col in performance_tables.columns if col.startswith(column_prefix)]

    # Get the columns
    plot_df = performance_tables[["Model"] + agreement_cols].copy()

    # Create a DataFrame for plotting
    plot_data = []
    for col in agreement_cols:

        # Extract agreement level (1.0, 0.5, 0.0)
        agreement_level = col.split('_')[-1]

        for row_idx in plot_df.index:
            plot_data.append({
                'model_name': plot_df.loc[row_idx, 'Model'],
                'agreement': agreement_level,
                'proportion': plot_df.loc[row_idx, col]
            })

    plot_df = pd.DataFrame(plot_data)

    # Sort models based on their mean value in the highest agreement category (1.0)
    highest_agreement = plot_df[plot_df["agreement"] == "1.0"].copy()
    highest_agreement = highest_agreement.groupby("model_name")["proportion"].mean().reset_index()
    sorted_models = highest_agreement.sort_values("proportion", ascending=False)[
        "model_name"
    ]

    # Convert model_name to categorical with specific order
    plot_df["model_name"] = pd.Categorical(
        plot_df["model_name"], categories=sorted_models, ordered=True
    )

    # Define the order of agreement categories (1.0, 0.5, 0.0 from left to right)
    reversed_categories = ["1.0", "0.5", "0.0"]

    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 8))

    # Plot the aggregated bars
    sns.barplot(
        data=plot_df,
        x="model_name",
        y="proportion",
        hue="agreement",
        hue_order=reversed_categories,
        ax=ax,
        order=sorted_models,
        errorbar='sd',
        capsize=0.2,
        err_kws={"linewidth": 0.5, "color": "black"},
    )

    # Add individual data points (can be done with stripplot)
    sns.stripplot(
        data=plot_df,
        x="model_name",
        y="proportion",
        hue="agreement",
        hue_order=reversed_categories,
        dodge=True,
        ax=ax,
        order=sorted_models,
        alpha=1.0,
        size=1.5,
        color="black",
    )

    # Add proportion labels on top of each bar
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.text(
                p.get_x() + p.get_width() / 2.0,
                height + 0.01,
                f"{height * 100:.0f}%",
                ha="center",
            )

    # Rotate x-axis tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # Set labels and title
    axis_type = "Cell Types" if by_cell_types else "Cells"
    ax.set_xlabel("Model")
    ax.set_ylabel("Proportion")
    ax.set_title(f"Proportion of Agreement Categories by Model (% of {axis_type})")
    ax.set_ylim(0, 1.05)
    ax.legend(title="Agreement Categories", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    return fig, ax


def main():
    parser = argparse.ArgumentParser(
        description="Create agreement plots from performance tables."
    )
    parser.add_argument(
        "--input-tables",
        nargs="+",
        required=True,
        help="Input performance table pickle files",
    )

    args = parser.parse_args()

    # Load manual cell type column
    with open("./dat/manual_cell_type_col.pkl", "rb") as f:
        manual_cell_type_col = pickle.load(f)

    # Load performance tables
    performance_tables = [pickle.load(open(f, 'rb')) for f in args.input_tables]

    # Reset the index of each table
    [table.reset_index(inplace=True) for table in performance_tables] # pylint: disable=expression-not-assigned

    # Concatenate the tables to create a single DataFrame
    performance_tables = pd.concat(performance_tables)
    performance_tables.reset_index(drop=True, inplace=True)

    # Create agreement plot for % of Cells
    agreement_plot_categorical = create_agreement_plot(
        performance_tables, by_cell_types=False
    )

    # Customize and save the % of Cells plot
    customize_figure(
        agreement_plot_categorical,
        remove_legend=True,
        x_tick_substrings=[
            "categorical_agreement_consistent_including_manual_"
            + manual_cell_type_col
            + "_consistent_including_manual_",
            "_ai_cell_type",
        ],
        new_ylabel="Agreement with Manual Annotation (by level of agreement)",
        remove_bar_labels=True,
        fig_width=2.4,
        fig_height=3,
        new_tick_labels=model_tick_labels,
    )

    agreement_plot_categorical[0].savefig(
        "res/12_aggregated_agreement_plots/agreement_plot_categorical.svg", format="svg"
    )

    # Create agreement plot for % of Cell Types
    agreement_plot_categorical_unweighted = create_agreement_plot(
        performance_tables, by_cell_types=True
    )

    # Customize and save the % of Cell Types plot
    customize_figure(
        agreement_plot_categorical_unweighted,
        remove_legend=True,
        x_tick_substrings=[
            "categorical_agreement_consistent_including_manual_"
            + manual_cell_type_col
            + "_consistent_including_manual_",
            "_ai_cell_type",
        ],
        new_ylabel="Agreement with Manual Annotation (by level of agreement)",
        remove_bar_labels=True,
        fig_width=2.4,
        fig_height=3,
        new_tick_labels=model_tick_labels,
    )
    agreement_plot_categorical_unweighted[0].savefig(
        "res/12_aggregated_agreement_plots/agreement_plot_categorical_unweighted.svg",
        format="svg",
    )


if __name__ == "__main__":
    main()
