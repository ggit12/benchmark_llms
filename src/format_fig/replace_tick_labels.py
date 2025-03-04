"""
This module contains a function to replace string tick labels on a specific axis with custom labels.
"""

def replace_tick_labels(
    ax,
    tick_dict,
    axis='x'
) -> None:
    """
    Replace string tick labels on a specific axis with custom labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis on which to replace tick labels.

    tick_dict : dict
        A dictionary where the keys are existing tick labels (strings) and the values are the desired new labels.

    axis : str, optional
        Which axis to modify ('x' or 'y'). Default is 'x'.
    """
    if axis == 'x':
        tick_labels = [tick.get_text() for tick in ax.get_xticklabels()]
        new_tick_labels = [tick_dict.get(label, label) for label in tick_labels]
        ax.set_xticklabels(new_tick_labels)
    elif axis == 'y':
        tick_labels = [tick.get_text() for tick in ax.get_yticklabels()]
        new_tick_labels = [tick_dict.get(label, label) for label in tick_labels]
        ax.set_yticklabels(new_tick_labels)
