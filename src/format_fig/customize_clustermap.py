"""
This module contains a function to customize a Seaborn clustermap plot.
"""
#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
#pylint: disable=too-many-branches
#pylint: disable=too-many-statements
#pylint: disable=line-too-long

import numpy as np
import matplotlib.pyplot as plt

from .replace_tick_labels import replace_tick_labels


def customize_clustermap(
    cluster_grid,
    remove_legend=False,
    x_tick_substrings=None,
    y_tick_substrings=None,
    new_ylabel=None,
    new_xlabel=None,
    fig_width=3.5,
    fig_height=3.5,
    remove_value_labels=False,
    new_tick_labels=None,
):

    """
    Customize a Seaborn clustermap plot.

    Parameters
    ----------
    cluster_grid : seaborn.matrix.ClusterGrid
        The Seaborn ClusterGrid object to customize.

    remove_legend : bool, optional
        Whether to remove the legend. Default is False.

    x_tick_substrings : list of str, optional
        A list of substrings to remove from x-axis tick labels. Default is None.

    y_tick_substrings : list of str, optional
        A list of substrings to remove from y-axis tick labels. Default is None.

    new_ylabel : str, optional
        The new y-axis label. Default is None.

    new_xlabel : str, optional
        The new x-axis label. Default is None.

    fig_width : float, optional
        The width of the figure. Default is 3.5.

    fig_height : float, optional
        The height of the figure. Default is 3.5.

    remove_value_labels : bool, optional
        Whether to remove value labels in each cell. Default is False.

    new_tick_labels : dict, optional
        A dictionary where the keys are existing tick labels (strings) and the values are the desired new labels.
        Default is None.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The customized figure and axes.
    """

    # Unpack the figure and axes from the ClusterGrid
    fig, ax = cluster_grid.fig, cluster_grid.ax_heatmap

    # Get the plot data
    plot_data = cluster_grid.data2d

    # Get the number of rows and columns
    _, num_cols = plot_data.shape

    # Set figure size
    fig.set_size_inches(fig_width, fig_height)

    # Set font properties
    font_size = 5
    font_weight = "bold"
    font_family = "sans-serif"

    # Remove plot title
    ax.set_title("")

    # Update x-axis and y-axis label font
    ax.xaxis.label.set_fontsize(font_size)
    ax.xaxis.label.set_fontweight(font_weight)
    ax.xaxis.label.set_fontfamily(font_family)
    ax.yaxis.label.set_fontsize(font_size)
    ax.yaxis.label.set_fontweight(font_weight)
    ax.yaxis.label.set_fontfamily(font_family)

    # Set ticks and labels for all columns
    ax.set_xticks(np.arange(num_cols) + 0.5)
    ax.set_xticklabels(plot_data.columns, rotation=90, ha="center")
    ax.tick_params(axis="x", which="major", labelsize=5, labelbottom=True)

    # Optionally update the y-axis label
    if new_ylabel is not None:
        ax.set_ylabel(
            new_ylabel, fontsize=font_size, fontweight=font_weight, fontname=font_family
        )

    # Optionally update the x-axis label
    if new_xlabel is not None:
        ax.set_xlabel(
            new_xlabel, fontsize=font_size, fontweight=font_weight, fontname=font_family
        )

    # Update tick label font sizes and optionally remove substrings from x-axis ticks
    x_ticklabels = ax.get_xticklabels()
    for tick in x_ticklabels:
        tick.set_fontsize(font_size)
        tick.set_fontweight(font_weight)
        tick.set_fontname(font_family)
        if x_tick_substrings is not None:
            tick_text = tick.get_text()
            for substring in x_tick_substrings:
                tick_text = tick_text.replace(substring, "")
            tick.set_text(tick_text)

    # Redraw the tick labels if modified
    if x_tick_substrings is not None:
        ax.set_xticklabels([tick.get_text() for tick in x_ticklabels])

    # Set x-axis tick labels to diagonal (45 degrees)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Update y-axis tick label font size, bold, and optionally remove substrings from y-axis ticks
    y_ticklabels = ax.get_yticklabels()
    for tick in y_ticklabels:
        tick.set_fontsize(font_size)
        tick.set_fontweight(font_weight)
        tick.set_fontname(font_family)
        if y_tick_substrings is not None:
            tick_text = tick.get_text()
            for substring in y_tick_substrings:
                tick_text = tick_text.replace(substring, "")
            tick.set_text(tick_text)

    # Redraw the y-axis tick labels if modified
    if y_tick_substrings is not None:
        ax.set_yticklabels([tick.get_text() for tick in y_ticklabels])

    # update tick labels with plotting labels
    if new_tick_labels:
        replace_tick_labels(ax, new_tick_labels, axis="x")
        replace_tick_labels(ax, new_tick_labels, axis="y")

    # Remove value labels in each cell if requested
    if remove_value_labels:
        for text in ax.texts:
            text.set_visible(False)
    else:
        # Update the font size of any text annotations (assumed to be bar value labels), but make them not bold
        for text in ax.texts:
            text.set_fontsize(font_size)
            text.set_fontweight("normal")  # Make bar value labels not bold
            text.set_fontname(font_family)

    # Remove x-axis label
    ax.set_xlabel("")

    # Set all lines (axes, gridlines, and plot lines) to 0.5pt
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    # Remove grid lines
    ax.grid(False)

    # Set tick line width to 0.5pt
    ax.tick_params(width=0.5)

    # Optionally modify or remove the legend
    legend = ax.get_legend()
    if legend is not None:
        if remove_legend:
            legend.set_visible(False)
        else:
            # Modify legend labels to remove substrings if present
            legend_labels = [text.get_text() for text in legend.get_texts()]
            if x_tick_substrings is not None:
                new_legend_labels = []
                for label in legend_labels:
                    for substring in x_tick_substrings:
                        label = label.replace(substring, "")
                    new_legend_labels.append(label)
                ax.legend(new_legend_labels)  # Update legend with modified labels

            # Set the legend font to bold and Arial (or fallback)
            for text in legend.get_texts():
                text.set_fontsize(font_size)
                text.set_fontweight(font_weight)
                text.set_fontname(font_family)

    # Automatically adjust the layout to make sure everything fits within the figure size
    fig.tight_layout(pad=0.5)

    return fig, ax
