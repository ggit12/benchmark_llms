"""
This function customizes the legend labels of a matplotlib Axes object.
"""

def customize_legend_labels(
    fig_tuple,
    remove_legend: bool = False,
    label_map: dict = None,
):
    """
    Customize a figure and axes object.

    Parameters
    ----------
    fig_tuple : tuple
        A tuple containing the figure and axes objects to customize.

    remove_legend : bool, optional
        Whether to remove the legend. Default is False.

    label_map : dict, optional
        A dictionary where the keys are existing legend labels (strings) and the values are the desired new labels.
        Default is None, which means no changes to the legend labels.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The customized figure.
    """
    # Unpack the tuple to get the figure and axes
    fig, ax = fig_tuple

    # Optionally modify or remove the legend
    legend = ax.get_legend()
    if legend is not None:
        if remove_legend:
            legend.set_visible(False)
        else:
            # Modify legend labels to remove substrings if present
            legend_labels = [text.get_text() for text in legend.get_texts()]
            if label_map is not None:
                # Map existing labels to new labels using the provided label_map
                new_legend_labels = [label_map[label] for label in legend_labels]
                ax.legend(new_legend_labels)  # Update legend with modified labels

    return fig
