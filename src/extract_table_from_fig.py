"""
Functions for extracting data from plots.
"""

import pandas as pd

def extract_table_from_fig(fig_tuple, value_col_name=None, make_percent=True):
    """
    Extracts the data from a bar plot and returns it as a pandas DataFrame.
    """

    if not value_col_name:
        value_col_name = 'value'
    # Unpack the tuple to get the figure and axes
    _, ax = fig_tuple

    # Extract x-tick labels (the names)
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]

    # Extract the heights of the bars (bar values)
    bar_heights = [bar.get_height() for bar in ax.patches]

    if make_percent:
        # Multiply by 100, round to 2 decimal places, and add '%' sign
        formatted_values = [f"{height * 100:.2f}%" for height in bar_heights]
    else:
        # Round to 3 decimal places
        formatted_values = [f"{height:.3f}" for height in bar_heights]

    # Create the table as a pandas DataFrame
    df = pd.DataFrame({'Model': xtick_labels, value_col_name: formatted_values}).set_index('Model')

    # Convert DataFrame to an HTML table
    # html_table = df.to_html(index=False)  # Do not include the index in the HTML

    return df
