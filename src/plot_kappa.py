"""
Plots a clustermap of pairwise Kappa scores and a bar chart of average pairwise Kappa scores.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


def plot_pairwise_clustermap(
    kappa: dict,
    metric: str = "euclidean",
    method: str = "average"
) -> sns.matrix.ClusterGrid:
    """Plot a clustermap of pairwise Kappa scores.

    Available metrics for clustering:
    - ``'euclidean'``: Euclidean distance
    - ``'correlation'``: Correlation coefficient  
    - ``'cosine'``: Cosine similarity
    - ``'cityblock'``: Manhattan distance
    - ``'jaccard'``: Jaccard similarity coefficient

    Available methods for clustering:
    - ``'single'``: Nearest Point Algorithm
    - ``'complete'``: Farthest Point Algorithm  
    - ``'average'``: Unweighted Pair Group Method with Arithmetic Mean (UPGMA)
    - ``'weighted'``: Weighted Pair Group Method with Arithmetic Mean (WPGMA)
    - ``'centroid'``: Unweighted Pair Group Method with Centroid Averaging (UPGMC)
    - ``'median'``: Weighted Pair Group Method with Centroid Averaging (WPGMC)
    - ``'ward'``: Minimizes the variance of the clusters being merged

    Parameters
    ----------
    kappa
        A dictionary containing pairwise Kappa scores.

    metric
        The metric to use for clustering.

    method
        The method to use for clustering.

    Returns
    -------
    The plot as a seaborn ClusterGrid
    """

    # Extract unique model names from the keys
    models = set()
    for pair in kappa["pairwise"].keys():
        models.update(pair)
    models = sorted(list(models))

    # Create an empty DataFrame
    df = pd.DataFrame(index=models, columns=models, dtype=float)

    # Fill the DataFrame with pairwise scores
    for (model1, model2), score in kappa["pairwise"].items():
        df.loc[model1, model2] = float(score)
        df.loc[model2, model1] = float(score)  # Symmetry

    # Fill diagonal with 1 (for clustering purposes)
    np.fill_diagonal(df.values, 1)

    # Perform clustering
    row_linkage = linkage(pdist(df, metric=metric), method=method)
    col_linkage = linkage(pdist(df.T, metric=metric), method=method)

    # Create a mask for the diagonal
    mask = np.eye(df.shape[0], dtype=bool)

    # Create the clustermap
    plt.figure(figsize=(20, 18))
    g = sns.clustermap(
        df,
        cmap="viridis",
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        yticklabels=True,
        xticklabels=True,
        figsize=(18, 14),
        mask=mask,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cbar_pos=(0.02, 0.8, 0.05, 0.18),
        vmin=0,
        vmax=1,
    )

    # Color the diagonal grey
    for i in range(df.shape[0]):
        g.ax_heatmap.add_patch(
            plt.Rectangle(
                (i, i), 1, 1, fill=True, facecolor="lightgrey", edgecolor="none"
            )
        )

    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

    plt.title(
        f"Pairwise Kappa Scores Clustermap\nMetric: {metric}, Method: {method}", pad=50
    )
    plt.tight_layout()
    plt.show()

    return g


def plot_average_pairwise_barchart(
    kappa: dict
) -> tuple[plt.Figure, plt.Axes]:
    """
    Plot a bar chart of average pairwise Kappa scores.

    Parameters
    ----------
    kappa
        A dictionary containing pairwise Kappa scores.

    Returns
    -------
    The plot as a tuple of Figure and Axes.
    """
    # Convert the average_pairwise data to a DataFrame
    df = pd.DataFrame.from_dict(
        kappa["average_pairwise"], orient="index", columns=["Score"]
    )
    df = df.sort_values("Score", ascending=False)

    # Remove the None value if present
    df = df.dropna()

    # Create the bar chart
    fig, ax = plt.subplots(figsize=(7, 5))
    ax = sns.barplot(x=df.index, y="Score", data=df)
    ax.set_title("Average Pairwise Kappa Scores")
    plt.xticks(rotation=90)
    ax.set_ylabel("Kappa Score")

    # Add value labels on top of each bar
    for i, v in enumerate(df["Score"]):
        ax.text(i, v, f"{v:.2f}", ha="center", va="bottom")

    plt.tight_layout()

    # Return the figure and axes as a tuple
    return fig, ax
