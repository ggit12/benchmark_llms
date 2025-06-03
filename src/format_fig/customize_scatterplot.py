"""
Customize a scatterplot figure with various formatting options.
"""

import matplotlib.pyplot as plt


def customize_scatterplot(
    fig_tuple,
    remove_legend=False,
    x_tick_substrings=None,
    new_ylabel=None,
    new_xlabel=None,
    fig_size=3.5,
    fig_width=None,
    fig_height=None,
    preserve_sizes=True,
):
    """
    Customize a scatterplot figure with various formatting options.

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

    new_xlabel : str, optional
        The new x-axis label. Default is None.

    fig_size : float, optional
        The size of the figure. Default is 3.5.

    preserve_sizes : bool, optional
        Whether to preserve the original dot sizes. Default is True.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The customized figure.
    ax : matplotlib.axes.Axes
        The customized axes.

    """
    # Unpack the tuple to get the figure and axes
    fig, ax = fig_tuple

    # Set figure size to slightly less than 90mm wide and 90mm tall to include bounding box
    if fig_width is None and fig_height is None:
        fig_width = fig_size
        fig_height = fig_size
    fig.set_size_inches(
        fig_width, fig_height
    )  # 90mm = 3.54331 inches, slightly smaller to account for labels

    # Set font size to 5pt, font weight to bold, and font family to Arial or sans-serif fallback
    font_size = 7
    font_weight = "bold"
    font_family = "sans-serif"

    # Set global font properties for the figure (with fallback)
    # plt.rcParams['font.family'] = font_family if font_family in plt.rcParams['font.family'] else 'sans-serif'

    # Modify the dot size setting
    if not preserve_sizes:
        for collection in ax.collections:
            collection.set_sizes([10])

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
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

    # Update y-axis tick label font size, bold, and Arial (or fallback)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(font_size)
        tick.set_fontweight(font_weight)
        tick.set_fontname(font_family)

    # Update the font size of any text annotations (assumed to be bar value labels), but make them not bold
    for text in ax.texts:
        text.set_fontsize(font_size)
        text.set_fontweight("normal")  # Make bar value labels not bold
        text.set_fontname(font_family)

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

    # Add a horizontal dashed line at y=1
    ax.axhline(y=1, color="black", linestyle="--", linewidth=0.5)
    ax.axvline(x=1, color="black", linestyle="--", linewidth=0.5)

    # Automatically adjust the layout to make sure everything fits within the figure size
    fig.tight_layout(
        pad=0.5
    )  # Adjust padding so that the entire bounding box fits within the size

    return fig, ax
