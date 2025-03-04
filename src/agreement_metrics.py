"""
Module to calculate agreement metrics between raters and a ground truth column in an AnnData object.
"""
# pylint: disable=line-too-long

from typing import Literal
import inspect

import numpy as np
import pandas as pd


def cell_type_by_plurality(adata, rater_cols, new_col_name="plurality_of_raters"):
    """This is a function to calculate a label for a cell based on the plurality vote of raters"""

    # Extract the relevant rater columns from obs
    raters = adata.obs[rater_cols]

    # Stack the DataFrame to transform it into a Series with a MultiIndex (row_index, column_name)
    stacked = raters.stack()

    # Compute the counts of each value per row
    counts = stacked.groupby(level=0).value_counts()

    # For each row, find the value with the highest count
    plurality_vote = counts.groupby(level=0).idxmax().apply(lambda x: x[1])

    # Add the plurality vote as a new column in adata.obs
    adata.obs[new_col_name] = plurality_vote


def row_wise_agreement_with_ground_truth(adata, rater_cols, manual_col):
    """This function calculates the average agreement of raters with a ground truth col"""
    def agreement_with_ground_truth(row):
        manual_label = row[manual_col]
        rater_labels = row[rater_cols]
        agreements = rater_labels == manual_label
        return agreements.mean()

    # Apply the agreement calculation row-wise
    agreement_results = adata.obs.apply(agreement_with_ground_truth, axis=1)

    return agreement_results


# The following functions are different ways to assess agreement amongst the raters
# for each type of agreement, there are two functions:
# 1) a function that operates on a single row of anndata
# 2) a function that caluclates the average within each category of a column in anndata

# There are three types of agreement (each calculated on a single row)
# 1) Plurality: The proportion of raters that agreed with the plurality of the raters
# 2) Gini: gini coefficient of all raters


def row_wise_plurality_agreement(adata, rater_cols):
    """This function calculates the agreement with the plurality label for each row"""
    def agreement_with_plurality(row, plurality_label):
        return (row == plurality_label).mean()

    plurality_label_col = "temp_plurality_label"
    cell_type_by_plurality(adata, rater_cols, new_col_name=plurality_label_col)

    # Extract the rater columns and plurality label from adata.obs
    ratings = adata.obs[rater_cols]
    plurality_labels = adata.obs[plurality_label_col]

    # Calculate agreement for each row
    plurality_agreements = ratings.apply(
        lambda row: agreement_with_plurality(row, plurality_labels.loc[row.name]),
        axis=1,
    )

    # Remove temp plurality label col
    adata.obs.drop(columns=[plurality_label_col], inplace=True)

    return plurality_agreements


def row_wise_gini(adata, rater_cols):
    """This function calculates the Gini coefficient for each row of ratings"""
    def gini(x):
        # Convert categorical labels to numeric
        _, numeric_x = np.unique(x, return_inverse=True)
        return 1 - np.sum(np.square(np.bincount(numeric_x) / len(x)))

    # Extract the rater columns from adata.obs
    ratings = adata.obs[rater_cols]

    # Calculate Gini coefficient for each row
    gini_coefficients = ratings.apply(gini, axis=1)

    return gini_coefficients


def apply_agreement_and_average(
    adata, rater_cols, manual_col, agreement_func, agreement_col
):
    """This function applies an agreement func to adata"""
    # Get the parameter names of the agreement_func
    func_params = list(inspect.signature(agreement_func).parameters.keys())

    # Prepare the arguments based on the function's parameters
    func_args = {}
    if "adata" in func_params:
        func_args["adata"] = adata
    if "rater_cols" in func_params:
        func_args["rater_cols"] = rater_cols
    if "manual_col" in func_params:
        func_args["manual_col"] = manual_col

    # Apply the agreement function and assign results to adata.obs
    try:
        adata.obs[agreement_col] = agreement_func(**func_args)
    except TypeError as e:
        # pylint: disable=raise-missing-from
        raise ValueError(
            f"Error calling agreement function: {str(e)}. "
            f"Expected parameters: {func_params}"
        )

    # Group by manual_col and calculate the mean of agreement_col
    grouped_averages = adata.obs.groupby(manual_col)[agreement_col].mean()

    return grouped_averages


# This is the final function that will calculate the agreements
def compute_agreement_df(
    adata,
    rater_cols,
    manual_col,
    agreement_type: Literal["plurality", "gini"],
    normalize_values=False,
):
    """
    This function calculates the agreement between raters and a ground truth column in an AnnData object using a the specified agreement metric.
    """

    # Define a dictionary to map agreement_type to the corresponding function
    agreement_functions = {
        "plurality": row_wise_plurality_agreement,
        "gini": row_wise_gini,
    }

    # 1. Calculate agreement with ground truth
    ground_truth_agreement = apply_agreement_and_average(
        adata,
        rater_cols,
        manual_col,
        agreement_func=row_wise_agreement_with_ground_truth,
        agreement_col="avg_ground_truth_agreement",
    )

    # 2. Calculate agreement as specified by agreement_type
    specified_agreement = apply_agreement_and_average(
        adata,
        rater_cols,
        manual_col,
        agreement_func=agreement_functions[agreement_type],
        agreement_col=f"{agreement_type}_agreement",
    )

    # Normalize if true
    if normalize_values:
        ground_truth_agreement = ground_truth_agreement / max(ground_truth_agreement)
        specified_agreement = specified_agreement / max(specified_agreement)

    # Create a DataFrame of agreements
    agreement_df = pd.DataFrame(
        {
            "ground_truth_agreement": ground_truth_agreement,
            "inter_rater_agreement": specified_agreement,
        }
    )

    return agreement_df
