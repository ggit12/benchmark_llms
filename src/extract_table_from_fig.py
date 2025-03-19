"""
Functions for extracting data from plots.
"""

import pandas as pd

def extract_table_from_fig(
    fig_tuple: tuple,
    value_col_name: str | None = None,
    make_percent: bool = True,
    round_to: int | None = 3
) -> pd.DataFrame:
    """
    Extracts data from a bar plot and returns it as a pandas :class:`DataFrame`.

    Parameters
    ----------
    fig_tuple
        Tuple containing the figure and axes objects.

    value_col_name
        Name of the column in the returned :class:`DataFrame` that contains the bar values.

    make_percent
        Whether to convert the values to percentages.
    
    round_to
        Number of decimal places to round the values to.
        If not ``None``, values will be returned as strings.
        If ``None``, values will not be rounded and will be returned as floats.

    Returns
    -------
    DataFrame containing the extracted data with model names as index
    and values from the bars as columns.

    Notes
    -----
    - If multiple bar groups exist, column names will be formatted as
      ``"{value_col_name}_{legend_label}"``
    - If ``make_percent`` is ``True``, values will be multiplied by 100 and 
      formatted with "%" suffix
    - If ``round_to`` is not ``None``, values will be rounded to specified 
      decimal places
    - ``round_to`` controls the return type of the value column(s).
    """

    if not value_col_name:
        value_col_name = "value"

    _, ax = fig_tuple

    # Grab the tick labels from the x-axis
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]

    # Get legend handles and labels (assume they match up, in order, with the bar containers)
    _, labels = ax.get_legend_handles_labels()

    # Prepare a dict starting with the 'Model' column
    data_dict = {"Model": xtick_labels}

    # Zip the containers and legend labels, then iterate over them
    for label, container in zip(labels, ax.containers):
        # Use a fallback if the label was never set (or is _nolegend_)
        if not label or label == "_nolegend_":
            label = value_col_name

        # Extract bar heights
        heights = [bar.get_height() for bar in container]

        # Format values
        if make_percent:
            formatted_values = [f"{h * 100:.2f}%" for h in heights]
        else:
            if round_to is not None:
                format_spec = f"{{:.{round_to}f}}"
                formatted_values = [format_spec.format(h) for h in heights]
            else:
                formatted_values = heights

        # Store the column using the legend label
        column_name = f"{value_col_name}_{label}" if len(labels) != 1 else value_col_name
        data_dict[column_name] = formatted_values

    # Create the DataFrame
    df = pd.DataFrame(data_dict).set_index("Model")
    return df
