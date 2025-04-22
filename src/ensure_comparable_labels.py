"""
This function creates new copies of the specified columns in the DataFrame, 
ensuring that all columns share the same set of categories.
"""
import pandas as pd
from anndata import AnnData

def ensure_comparable_labels_adata(
    adata: AnnData,
    columns: list[str],
    new_col_prefix: str = "shared_categories"
) -> AnnData:
    """
    Ensure all specified columns in the :class:`AnnData` object have the same set of categories.

    Parameters
    ----------
    adata
        The :class:`AnnData` object containing the data.

    columns
        The list of column names to be processed.

    new_col_prefix : str
        The prefix for the new categorical columns created.

    Returns
    -------
    The :class:`AnnData` object with new categorical columns added.
    """

    # Call the main function on the adata.obs DataFrame
    ensure_comparable_labels_main(adata.obs, columns, new_col_prefix)
    return adata

def ensure_comparable_labels_main(
    df: pd.DataFrame,
    columns: list[str],
    new_col_prefix: str
) -> pd.DataFrame:
    """
    Ensure all specified columns in the :class:`DataFrame` have the same set of categories.

    Parameters
    ----------
    df
        The DataFrame containing the columns to be processed.

    columns
        The list of column names to be processed.

    new_col_prefix : str
        The prefix for the new categorical columns created.

    Returns
    -------
    The :class:`DataFrame` with new categorical columns added.
    """
    # Gather the union of all unique categories from the specified columns
    all_categories = set()
    for col in columns:
        all_categories.update(df[col].unique())
    all_categories = sorted(all_categories)

    # For each column, create a new categorical column with the union of categories
    for col in columns:
        new_col_name = f"{new_col_prefix}_{col}"
        df[new_col_name] = pd.Categorical(df[col], categories=all_categories)

    return df
