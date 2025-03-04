"""
This module finds indices that meet criteria in a dataframe.
"""

import numpy as np

def find_indices_closest_to_4_corners(agreement_df, n):
    """
    This function finds the indices of the n rows in the agreement_df that are closest to the four corners of the plot.
    """
    # Define the four corners with their names
    corners = {
        "top_left": (0, 1),
        "bottom_left": (0, 0),
        "top_right": (1, 1),
        "bottom_right": (1, 0),
    }

    # Identify the inter-rater agreement column (the one that's not 'ground_truth_agreement')
    inter_rater_agreement_col = [
        col for col in agreement_df.columns if col != "ground_truth_agreement"
    ][0]

    # Initialize the result dictionary
    result = {}

    # Calculate the Euclidean distance from each corner and get indices
    for name, corner in corners.items():
        distances = np.sqrt(
            (agreement_df["ground_truth_agreement"] - corner[0]) ** 2
            + (agreement_df[inter_rater_agreement_col] - corner[1]) ** 2
        )
        # Get the indices of the n smallest distances
        indices = distances.nsmallest(n).index.tolist()
        # Assign to the result dictionary
        result[name] = indices

    # Return the dictionary of indices
    return result
