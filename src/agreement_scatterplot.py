"""
This module contains functions to plot agreement scatterplots.
"""

# pylint: disable=line-too-long
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from adjustText import adjust_text
from scipy.stats import gaussian_kde

from .find_indices import find_indices_closest_to_4_corners

def calculate_weights(adata, col, equal_weights=False):
    """
    Calculate weights based on the proportion of each category in adata.obs[col].

    Parameters
    -----------
    adata: AnnData object
    col: Column name in adata.obs (categorical)
    equal_weights: If True, return equal weights for all observations

    Returns
    --------
    weights: pandas Series with weights for each observation
    """

    if equal_weights:
        # Assign equal weight to each observation
        category_proportions = pd.Series(1, index=adata.obs[col].unique())
    else:
        # Calculate the proportion of each category
        category_proportions = adata.obs[col].value_counts(normalize=True)
        # Map the proportion to each observation
        # weights = adata.obs[col].map(category_proportions)
    return category_proportions


def plot_agreement(
    agreement_df, weights, show_legend=True, show_labels=True, normalize_kdes=True
):
    """
    Plots the agreement between ground truth and specified agreement, with points colored based on weight groups.
    KDEs are weighted by the provided weights and can be normalized so that all KDEs sum to 1.
    When normalize_kdes is False, KDEs are computed without weights (each point is weighted equally).

    Parameters
    ----------
    agreement_df: DataFrame containing 'ground_truth_agreement' and one other agreement column.
    weights: Series containing weights for each observation, indexed to match agreement_df.
    show_legend: Boolean to control whether legends are displayed.
    show_labels: Boolean to control whether point labels are displayed.
    normalize_kdes: Boolean to control whether to scale KDEs by their total weights and compute weighted KDEs.

    Returns
    -------
    fig: The matplotlib figure object.
    ax_scatter: The main axes object for the scatter plot.
    """
    # Identify the inter-rater agreement column (the one that's not 'ground_truth_agreement')
    inter_rater_agreement_col = [
        col for col in agreement_df.columns if col != "ground_truth_agreement"
    ][0]

    # Extract the agreement type from the column name
    agreement_type = inter_rater_agreement_col.replace("_agreement", "").capitalize()

    # Get indices of points closest to the four corners
    top_indices_dict = find_indices_closest_to_4_corners(agreement_df, n=3)

    # Find unique indices from the returned dictionary
    top_indices = list(set(sum(top_indices_dict.values(), [])))

    ground_truth_agreement = agreement_df["ground_truth_agreement"]
    specified_agreement = agreement_df[inter_rater_agreement_col]

    # Ensure that weights align with the data
    weights = weights.loc[agreement_df.index]

    # Split the data into three groups based on weights
    quantiles = weights.quantile([0, 1 / 3, 2 / 3, 1]).values
    group_labels = ["Low", "Medium", "High"]
    weights_group = pd.cut(
        weights, bins=quantiles, labels=group_labels, include_lowest=True
    )

    # Define colors for the groups
    group_colors = {"Low": "red", "Medium": "orange", "High": "green"}

    # Set up the figure and axes with gridspec
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(4, 4, hspace=0.05, wspace=0.05)
    ax_histx = fig.add_subplot(gs[0, :-1])
    ax_scatter = fig.add_subplot(gs[1:, :-1])
    ax_histy = fig.add_subplot(gs[1:, -1], sharey=ax_scatter)

    # Create the scatter plot, coloring points by weight groups
    scatter_plots = []
    for group in group_labels:
        idx = weights_group == group
        scatter = ax_scatter.scatter(
            ground_truth_agreement[idx],
            specified_agreement[idx],
            color=group_colors[group],
            alpha=0.6,
            label=f"{group} Weights" if show_legend else None,
            s=(
                (weights[idx] * 1000 + 10) if normalize_kdes else 10
            ),  # Adjust size based on weights
            edgecolors="black",
            linewidths=1,
        )
        scatter_plots.append(scatter)

    ax_scatter.set_xlabel("Agreement with Ground Truth")
    ax_scatter.set_ylabel(f"{agreement_type} Agreement")

    # Set the axis limits to go from (-0.1, -0.1) to (1.1, 1.1)
    ax_scatter.set_xlim(-0.1, 1.1)
    ax_scatter.set_ylim(-0.1, 1.1)

    # Add a legend to the scatter plot if requested
    if show_legend:
        ax_scatter.legend(title="Weight Groups", loc="upper left")

    # Label the points if requested
    if show_labels:
        texts = []
        for idx in top_indices:
            # Generate small random jitter for label positions
            jitter_x = random.uniform(-0.01, 0.01)
            jitter_y = random.uniform(-0.01, 0.01)

            # Apply jitter to the label positions
            label_x = ground_truth_agreement.loc[idx] + jitter_x
            label_y = specified_agreement.loc[idx] + jitter_y

            texts.append(
                ax_scatter.text(label_x, label_y, str(idx), color="black", fontsize=12)
            )

        # Adjust the text labels to minimize overlap
        adjust_text(
            texts,
            x=ground_truth_agreement.loc[top_indices],
            y=specified_agreement.loc[top_indices],
            arrowprops=dict(arrowstyle="->", color="gray"),
            ax=ax_scatter,
            lim=100,
        )

    # Draw quadrants
    ax_scatter.axhline(0.5, color="green", linestyle="--")  # Horizontal line at y=0.5
    ax_scatter.axvline(0.5, color="green", linestyle="--")  # Vertical line at x=0.5

    # Calculate total weights for normalization if needed
    if normalize_kdes:
        total_weight = weights.sum()
        group_total_weights = {}
        for group in group_labels:
            group_indices = weights_group == group
            group_total_weights[group] = weights[group_indices].sum()

    # Plot the marginal density plots for each group
    for group in group_labels:
        group_indices = weights_group == group
        color = group_colors[group]

        # Prepare data for KDE
        data_x = ground_truth_agreement[group_indices]
        data_y = specified_agreement[group_indices]

        # Set up KDEs
        if normalize_kdes:
            # Weighted KDEs
            kde_x = gaussian_kde(data_x, weights=weights[group_indices])
            kde_y = gaussian_kde(data_y, weights=weights[group_indices])
        else:
            # Unweighted KDEs
            kde_x = gaussian_kde(data_x)
            kde_y = gaussian_kde(data_y)

        kde_x.set_bandwidth(
            kde_x.factor * 0.5
        )  # Adjust bandwidth to match original bw_adjust=0.5
        kde_y.set_bandwidth(kde_y.factor * 0.5)

        x_range = np.linspace(-0.1, 1.1, 100)
        y_range = np.linspace(-0.1, 1.1, 100)

        density_x = kde_x(x_range)
        density_y = kde_y(y_range)

        # Normalize densities if requested
        if normalize_kdes:
            scaling_factor = group_total_weights[group] / total_weight
            density_x *= scaling_factor
            density_y *= scaling_factor

        # Plot the densities
        ax_histx.fill_between(
            x_range,
            density_x,
            color=color,
            alpha=0.6,
            label=f"{group} Weights" if show_legend else None,
        )
        ax_histy.fill_betweenx(y_range, density_y, color=color, alpha=0.6)

    # Hide the spines and tick labels for the marginal plots
    ax_histx.axis("off")
    ax_histy.axis("off")

    # Add legend to the marginal x-density plot if requested
    if show_legend:
        ax_histx.legend(loc="upper right")

    plt.tight_layout()
    plt.show()

    return fig, ax_scatter


def plot_agreement_simple(agreement_df):
    """Makes a scatter plot of agreement with ground truth vs. agreement with specified rater."""
    # Identify the inter-rater agreement column (the one that's not 'ground_truth_agreement')
    inter_rater_agreement_col = [
        col for col in agreement_df.columns if col != "ground_truth_agreement"
    ][0]

    # Extract the agreement type from the column name
    agreement_type = inter_rater_agreement_col.replace("_agreement", "").capitalize()

    # Get indices of points closest to the four corners
    top_indices = find_indices_closest_to_4_corners(agreement_df, n=3)

    # find_indices_closest_to_4_corners returns a dict, so get unique values from it
    top_indices = list(set(sum(top_indices.values(), [])))

    ground_truth_agreement = agreement_df["ground_truth_agreement"]
    specified_agreement = agreement_df[inter_rater_agreement_col]

    # Create a scatter plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(ground_truth_agreement, specified_agreement)
    ax.set_xlabel("Agreement with Ground Truth")
    ax.set_ylabel(f"{agreement_type} Agreement")
    ax.set_title(f"Agreement with Ground Truth vs {agreement_type} Agreement")

    # Set the axis limits to go from (0, 0) to (1, 1)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)

    texts = []
    # Label the points
    for idx in top_indices:
        # Generate small random jitter for label positions
        jitter_x = random.uniform(-0.01, 0.01)
        jitter_y = random.uniform(-0.01, 0.01)

        # Apply jitter to the label positions
        label_x = ground_truth_agreement.loc[idx] + jitter_x
        label_y = specified_agreement.loc[idx] + jitter_y

        texts.append(ax.text(label_x, label_y, str(idx), color="black", fontsize=12))

    # Adjust the text labels to minimize overlap
    adjust_text(
        texts,
        x=ground_truth_agreement.loc[top_indices],
        y=specified_agreement.loc[top_indices],
        arrowprops=dict(arrowstyle="->", color="gray"),
        # expand_text=(1.2, 1.2),
        # expand_points=(1.2, 1.2),
        # force_text=(0.5, 0.5),
        # force_points=(0.2, 0.2),
        lim=100,
    )

    # Draw quadrants
    ax.axhline(0.5, color="g", linestyle="-")  # Horizontal line at y=0.5
    ax.axvline(0.5, color="g", linestyle="-")  # Vertical line at x=0.5

    plt.legend()
    plt.tight_layout()
    plt.show()

    return fig, ax
