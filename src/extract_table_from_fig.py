"""
Functions for extracting data from plots.
"""

import pandas as pd

def extract_table_from_fig(
    fig_tuple: tuple,
    value_col_name: str | None = None
) -> pd.DataFrame:
    """
    Extracts data from a bar plot and returns it as a pandas :class:`DataFrame`.

    Parameters
    ----------
    fig_tuple
        Tuple containing the figure and axes objects.

    value_col_name
        Name of the column in the returned :class:`DataFrame` that contains the bar values.

    Returns
    -------
    DataFrame containing the extracted data with model names as index
    and values from the bars as columns.

    Notes
    -----
    If multiple bar groups exist, column names will be formatted as
    ``"{value_col_name}_{legend_label}"``
    """

    if not value_col_name:
        value_col_name = "value"

    _, ax = fig_tuple

    # Grab the tick labels from the x-axis
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]

    # Get legend handles and labels (assume they match up, in order, with the bar containers)
    _, labels = ax.get_legend_handles_labels()

    # Handle empty labels (i.e. when labels is `[]`)
    if not labels:
        labels = [f"container_{i+1}" for i in range(len(ax.containers))]

    # Prepare a dict starting with the 'Model' column
    data_dict = {"Model": xtick_labels}

    # Zip the containers and legend labels, then iterate over them
    for label, container in zip(labels, ax.containers):

        # Extract bar heights
        heights = [bar_i.get_height() for bar_i in container]

        # Store the column using the legend label
        column_name = f"{value_col_name}_{label}" if len(labels) != 1 else value_col_name
        data_dict[column_name] = heights

    # Create the DataFrame
    df = pd.DataFrame(data_dict).set_index("Model")
    return df
