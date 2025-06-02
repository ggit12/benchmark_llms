"""
This function compares cell type labels by direct string comparison.
"""

from anndata import AnnData

def direct_compare_cell_type_labels_pairwise(
    adata: AnnData,
    cols1: list[str],
    cols2: list[str],
    new_col_prefix: str = 'agreement'
) -> None:
    """
    Compare cell type labels by finding unique combinations between 
    labels in ``cols1`` and ``cols2``, applying the comparison, and 
    mapping the results back to ``adata.obs``.

    Parameters
    -----------
    adata
        an :class:`AnnData`.

    cols1:
        :class:`List` of columns to compare against cols2.

    cols2:
        :class:`List` of columns to compare with cols1.

    new_col_prefix:
        The base name for the new comparison result columns.

    Returns
    --------
    ``None``

    Notes
    -----
    This function modifies ``adata`` in-place by creating new columns in ``adata.obs`` for each unique
    combination of columns in ``cols1`` and ``cols2``. The new columns
    are named using the format ``"{new_col_prefix}_{col1}_{col2}"``.
    The values in these columns are boolean, indicating whether the
    labels in the two columns are the same.
    """
    # Loop over all combinations of columns
    # and create an agreement column for each pair
    for col1 in cols1:
        for col2 in cols2:
            if col1 == col2:
                continue

            new_col_name = f"{new_col_prefix}_{col1}_{col2}"
            agreement = adata.obs[col1].astype(str) == adata.obs[col2].astype(str)
            adata.obs[new_col_name] = agreement
