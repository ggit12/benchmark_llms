"""
This module contains a function to customize a fig, ax tuple.
"""
#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
#pylint: disable=too-many-branches
#pylint: disable=too-many-statements
#pylint: disable=line-too-long

import numpy as np
import matplotlib.pyplot as plt

from .replace_tick_labels import replace_tick_labels

def customize_figure(
    fig_tuple,
    remove_legend=False,
    x_tick_substrings=None,
    new_ylabel=None,
    fig_width=3.5,
    fig_height=3.5,
    remove_bar_labels=False,
    new_tick_labels=None,
):
    """
    Customize a figure and axes object.

    Parameters
    ----------
    fig_tuple : tuple
        A tuple containing the figure and axes objects to customize.

    remove_legend : bool, optional
        Whether to remove the legend. Default is False.

    x_tick_substrings : list of str, optional
        A list of substrings to remove from x-axis tick labels. Default is None.

    new_ylabel : str, optional
        The new y-axis label. Default is None.

    fig_width : float, optional
        The width of the figure. Default is 3.5.

    fig_height : float, optional
        The height of the figure. Default is 3.5.

    remove_bar_labels : bool, optional
        Whether to remove bar value labels. Default is False.

    new_tick_labels : dict, optional
        A dictionary where the keys are existing tick labels (strings) and the values are the desired new labels.
        Default is None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The customized figure.
    """
    # Unpack the tuple to get the figure and axes
    fig, ax = fig_tuple

    # Set figure size to slightly less than 90mm wide and 90mm tall to include bounding box
    fig.set_size_inches(
        fig_width, fig_height
    )  # 90mm = 3.54331 inches, slightly smaller to account for labels

    # Set font size to 5pt, font weight to bold, and font family to Arial or sans-serif fallback
    font_size = 7
    font_weight = "bold"
    font_family = "sans-serif"

    # Set global font properties for the figure (with fallback)
    # plt.rcParams['font.family'] = font_family if font_family in plt.rcParams['font.family'] else 'sans-serif'

    # Remove plot title
    ax.set_title("")

    # Update x-axis and y-axis label font with size, bold, and Arial (or fallback)
    ax.xaxis.label.set_fontsize(font_size)
    ax.xaxis.label.set_fontweight(font_weight)
    ax.xaxis.label.set_fontfamily(font_family)
    ax.yaxis.label.set_fontsize(font_size)
    ax.yaxis.label.set_fontweight(font_weight)
    ax.yaxis.label.set_fontfamily(font_family)

    # Optionally update the y-axis label
    if new_ylabel is not None:
        ax.set_ylabel(
            new_ylabel, fontsize=font_size, fontweight=font_weight, fontname=font_family
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

    # Set x-axis tick labels to vertical (90 degrees)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # update tick labels with plotting labels
    if new_tick_labels:
        replace_tick_labels(ax, new_tick_labels, axis="x")

    # Update y-axis tick label font size, bold, and Arial (or fallback)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(font_size)
        tick.set_fontweight(font_weight)
        tick.set_fontname(font_family)

    # Update the font size of any text annotations (assumed to be bar value labels), but make them not bold
    if remove_bar_labels:
        for text in ax.texts:
            text.remove()
    else:
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

    # Set y-axis ticks at increments of 0.2
    # y_min, y_max = ax.get_ylim()
    ax.set_yticks(np.arange(0, 1.05, 0.2))
    ax.set_ylim(0, 1.05)

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

    # Add a horizontal dashed line at y=1
    ax.axhline(
        y=1, color="black", linestyle="--", linewidth=0.5
    )  # Set line width to 0.5pt

    # Automatically adjust the layout to make sure everything fits within the figure size
    fig.tight_layout(
        pad=0.5
    )  # Adjust padding so that the entire bounding box fits within the size

    return fig
