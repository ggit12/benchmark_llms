"""
This Module contains functions that to plot model agreement comapred to each other and to ground truth.
"""

# pylint: disable=line-too-long
# pylint: disable=unused-variable
# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from anndata import AnnData


def plot_model_agreement_unweighted(
    adata, group_by, sub_group_by, model_cols, granularity=2
):
    """
    Plots the average values of specified model columns across varying levels of granularity,
    ensuring that each group (e.g., cell type) contributes equally to the average regardless of its size.

    Parameters:
    - adata: AnnData object containing the data.
    - group_by: str, key in adata.obs for the main grouping (e.g., 'cell_type').
    - sub_group_by: str, key in adata.obs for the sub-grouping (e.g., 'tissue').
    - model_cols: list of str, column names for the models (e.g., ['agreement_model_1', 'agreement_model_2']).
    - granularity: int, level of detail in the plot
                   (0 = models only,
                    1 = models within cell types,
                    2 = models within cell types and tissues).
    """
    # Validate input columns
    missing_cols = [col for col in model_cols if col not in adata.obs]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in adata.obs.")
    if group_by not in adata.obs:
        raise ValueError(f"Group key '{group_by}' not found in adata.obs.")
    if sub_group_by and sub_group_by not in adata.obs:
        raise ValueError(f"Sub-group key '{sub_group_by}' not found in adata.obs.")

    # Melt the DataFrame to long format
    melted = adata.obs.melt(
        id_vars=[group_by, sub_group_by] if sub_group_by else [group_by],
        value_vars=model_cols,
        var_name="model_name",
        value_name="agreement",
    )

    # Function to compute unweighted mean per model
    def compute_unweighted_mean(df, group_keys):
        """
        Computes the unweighted mean by first averaging within group_keys and then averaging those means.

        Parameters:
        - df: DataFrame to group.
        - group_keys: list of columns to group by first.

        Returns:
        - Series with the unweighted mean for each model.
        """
        # First, compute the mean agreement per model and group_keys
        per_group_mean = (
            df.groupby(["model_name"] + group_keys)["agreement"].mean().reset_index()
        )
        # Then, compute the mean across the group_keys for each model
        unweighted_mean = per_group_mean.groupby("model_name")["agreement"].mean()
        return unweighted_mean

    # Function to compute unweighted mean for heatmap (granularity=2)
    def compute_unweighted_heatmap_mean(df, group_keys):
        """
        Computes the unweighted mean for heatmap by averaging within group_keys and then reshaping.

        Parameters:
        - df: DataFrame to group.
        - group_keys: list of columns to group by first.

        Returns:
        - Pivoted DataFrame suitable for heatmap.
        """
        # First, compute the mean agreement per model and group_keys
        per_group_mean = (
            df.groupby(["model_name"] + group_keys)["agreement"].mean().reset_index()
        )
        # Then, pivot to have models as rows and group_keys as columns
        pivot_cols = group_keys[-1] if len(group_keys) == 2 else group_keys
        pivot = per_group_mean.pivot_table(
            index="model_name",
            columns=group_keys[1] if len(group_keys) == 2 else group_keys[0],
            values="agreement",
        )
        # Finally, compute the mean across the pivoted groups for each model
        unweighted_pivot = pivot.mean(axis=1)
        return unweighted_pivot

    if granularity == 0:
        # Unweighted average across all cell types
        # Step 1: Average per model per cell type
        grouped_means = compute_unweighted_mean(melted, [group_by])
        # Step 2: Overall average per model
        grouped_means = grouped_means.sort_values(ascending=False)

        # Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        grouped_means.plot(kind="bar", ax=ax, colormap="Paired")

        # Add value labels on top of each bar
        for i, v in enumerate(grouped_means):
            ax.text(i, v, f"{v * 100:.0f}%", ha="center", va="bottom")

        # Set labels and title
        ax.set_xlabel("Model")
        ax.set_ylabel("Unweighted Average Agreement")
        ax.set_title("Unweighted Average Model Agreement")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center")

        plt.tight_layout()
        return fig, ax

    elif granularity == 1:
        # Unweighted average within each cell type
        # Step 1: Average per model per cell type
        per_group_mean = (
            melted.groupby(["model_name", group_by])["agreement"].mean().reset_index()
        )
        # Step 2: Pivot to have models as rows and cell types as columns
        pivot = per_group_mean.pivot(
            index="model_name", columns=group_by, values="agreement"
        )
        # Step 3: Plotting the unweighted averages
        pivot.plot(kind="bar", figsize=(14, 8), colormap="Paired")

        # Add value labels on top of each bar for each cell type
        for idx, model in enumerate(pivot.index):
            for cell_type in pivot.columns:
                value = pivot.loc[model, cell_type]
                ax = plt.gca()
                x = idx
                y = pivot.loc[model, cell_type]
                ax.text(x, y, f"{y * 100:.0f}%", ha="center", va="bottom", fontsize=8)

        # Set labels and title
        plt.xlabel("Model")
        plt.ylabel("Unweighted Average Agreement by Cell Type")
        plt.title(f"Unweighted Average Model Agreement by {group_by.capitalize()}")
        plt.xticks(rotation=90, ha="center")
        plt.legend(
            title=group_by.capitalize(), bbox_to_anchor=(1.05, 1), loc="upper left"
        )

        plt.tight_layout()
        fig = plt.gcf()
        ax = plt.gca()
        return fig, ax

    elif granularity == 2:
        if not sub_group_by:
            raise ValueError("sub_group_by must be provided for granularity=2.")

        # Step 1: Average per model per cell type per tissue
        per_group_mean = (
            melted.groupby(["model_name", group_by, sub_group_by])["agreement"]
            .mean()
            .reset_index()
        )

        # Step 2: Compute unweighted averages across cell types and tissues
        # This involves averaging across all unique (cell_type, tissue) combinations
        # For visualization, we'll create a pivot table with models as rows and (cell_type, tissue) as columns
        pivot = per_group_mean.pivot_table(
            index="model_name", columns=[group_by, sub_group_by], values="agreement"
        )

        # Optional: If you want to compute an overall unweighted average per model across all groups
        # unweighted_pivot = pivot.mean(axis=1)

        # Create a mask for NaN values
        mask = pivot.isnull()

        # Define a color palette for tissues
        unique_tissues = pivot.columns.get_level_values(1).unique()
        num_tissues = len(unique_tissues)
        if num_tissues > 30:
            raise ValueError(
                f"Number of tissues ({num_tissues}) exceeds the supported number of colors (30)."
            )

        # Generate distinct colors for tissues
        tissue_palette = sns.color_palette("hsv", num_tissues)
        tissue_colors = {
            tissue: color for tissue, color in zip(unique_tissues, tissue_palette)
        }

        # Create a list of colors corresponding to each (cell_type, tissue) column based on tissue
        col_colors = [tissue_colors[tissue] for _, tissue in pivot.columns]

        # Plotting the heatmap
        cmap = sns.cm.rocket_r
        cmap.set_bad(color="lightgrey")  # Color for NaN values

        g = sns.clustermap(
            pivot,
            cmap=cmap,
            mask=mask,
            linewidths=0.5,
            figsize=(16, 10),
            col_colors=col_colors,
            dendrogram_ratio=(0.1, 0.2),
            cbar_kws={"label": "Unweighted Average Agreement"},
            xticklabels=True,
            yticklabels=True,
        )

        # Add vertical lines to separate tissues
        current_pos = 0
        for i, tissue in enumerate(unique_tissues):
            num_cols = len([col for col in pivot.columns if col[1] == tissue])
            current_pos += num_cols
            g.ax_col_dendrogram.axvline(current_pos, color="black", linewidth=0.5)

        # Create a legend for tissues
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                label=tissue,
                markersize=10,
                markerfacecolor=color,
            )
            for tissue, color in tissue_colors.items()
        ]
        g.ax_col_dendrogram.legend(
            handles=handles,
            title=sub_group_by.capitalize(),
            bbox_to_anchor=(1, 1),
            loc="upper left",
        )

        # Set titles
        plt.title("Unweighted Average Model Agreement by Cell Type and Tissue", pad=100)

        plt.tight_layout()
        return g

    else:
        raise ValueError("Granularity must be 0, 1, or 2.")


def plot_model_agreement_categorical_unweighted(
    adata: AnnData,
    group_by: str,
    sub_group_by: str | None,
    model_cols: list[str],
    granularity: int = 2,
) -> tuple:
    """
    Plots the relative proportions of categories within specified model columns across varying levels of granularity,
    ensuring that each group (e.g., cell type and tissue) contributes equally to the proportions regardless of group size.

    Parameters
    -----------
    adata
        An :class:`AnnData`.

    group_by
        key in adata.obs for the main grouping (e.g., 'cell_type').

    sub_group_by
        key in adata.obs for the sub-grouping (e.g., 'tissue'). Set to None if not used.

    model_cols
        column names for the models (e.g., ['model_1', 'model_2']). These should be categorical.

    granularity
        level of detail in the plot
                   (0 = models only,
                    1 = models within cell types,
                    2 = models within cell types and tissues).
    """
    # Validate input columns
    missing_cols = [col for col in model_cols if col not in adata.obs]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in adata.obs.")
    if group_by not in adata.obs:
        raise ValueError(f"Group key '{group_by}' not found in adata.obs.")
    if granularity >= 1 and sub_group_by not in adata.obs:
        raise ValueError(f"Sub-group key '{sub_group_by}' not found in adata.obs.")

    # Ensure that model_cols are categorical or convert numeric types to categories
    for col in model_cols:
        if not pd.api.types.is_categorical_dtype(adata.obs[col]):
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                adata.obs[col] = adata.obs[col].astype("category")
            else:
                raise ValueError(
                    f"Column '{col}' must be categorical or convertible to categorical."
                )

    # Melt the dataframe to long format
    id_vars = [group_by]
    if granularity == 2:
        id_vars.append(sub_group_by)
    elif granularity == 1:
        pass  # Only group_by is needed
    melted = adata.obs.melt(
        id_vars=id_vars,
        value_vars=model_cols,
        var_name="model_name",
        value_name="agreement",
    )

    # Ensure 'agreement' is categorical and reverse the order of categories
    if not pd.api.types.is_categorical_dtype(melted["agreement"]):
        melted["agreement"] = melted["agreement"].astype("category")

    # Reverse the order of 'agreement' categories for consistent plotting
    original_categories = melted["agreement"].cat.categories.tolist()
    reversed_categories = original_categories[::-1]
    melted["agreement"] = melted["agreement"].cat.reorder_categories(
        reversed_categories, ordered=True
    )

    if granularity == 0:
        # Granularity 0: Models only
        # Treat each group_by category equally by computing proportions per group and then averaging

        # Step 1: Compute counts per model, group_by, and agreement
        counts = (
            melted.groupby(["model_name", group_by, "agreement"])
            .size()
            .reset_index(name="count")
        )

        # Step 2: Compute proportions within each model and group_by
        counts["proportion"] = counts.groupby(["model_name", group_by])[
            "count"
        ].transform("sum")
        counts["proportion"] = counts["count"] / counts["proportion"]

        # Step 3: Compute unweighted average proportions across all group_by categories
        avg_counts = (
            counts.groupby(["model_name", "agreement"])["proportion"]
            .mean()
            .reset_index()
        )

        # Step 4: Sort models based on the average proportion of the highest agreement category
        sort_category = reversed_categories[0]
        sort_order = avg_counts[avg_counts["agreement"] == sort_category].set_index(
            "model_name"
        )["proportion"]
        sorted_models = sort_order.sort_values(ascending=False).index.tolist()
        avg_counts["model_name"] = pd.Categorical(
            avg_counts["model_name"], categories=sorted_models, ordered=True
        )

        # Step 5: Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        sns.barplot(
            data=avg_counts,
            x="model_name",
            y="proportion",
            hue="agreement",
            hue_order=reversed_categories,
            order=sorted_models,
            ax=ax,
        )

        # Add proportion labels on top of each bar
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.text(
                    p.get_x() + p.get_width() / 2.0,
                    height + 0.01,
                    f"{height * 100:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        # Customize plot
        ax.set_xlabel("Model")
        ax.set_ylabel("Unweighted Average Proportion")
        ax.set_title("Unweighted Proportion of Agreement Categories by Model")
        ax.set_ylim(0, 1.05)
        ax.legend(
            title="Agreement Categories", bbox_to_anchor=(1.05, 1), loc="upper left"
        )
        plt.xticks(rotation=90)
        plt.tight_layout()
        return fig, ax

    elif granularity == 1:
        # Granularity 1: Models within cell types
        # Ensure that each cell type contributes equally by averaging proportions across cell types

        # Step 1: Compute counts per model, group_by, and agreement
        counts = (
            melted.groupby(["model_name", group_by, "agreement"])
            .size()
            .reset_index(name="count")
        )

        # Step 2: Compute proportions within each model and group_by
        counts["proportion"] = counts.groupby(["model_name", group_by])[
            "count"
        ].transform("sum")
        counts["proportion"] = counts["count"] / counts["proportion"]

        # Step 3: Compute unweighted average proportions across all group_by categories
        avg_counts = (
            counts.groupby(["model_name", "agreement"])["proportion"]
            .mean()
            .reset_index()
        )

        # Step 4: Sort models based on the average proportion of the highest agreement category
        sort_category = reversed_categories[0]
        sort_order = avg_counts[avg_counts["agreement"] == sort_category].set_index(
            "model_name"
        )["proportion"]
        sorted_models = sort_order.sort_values(ascending=False).index.tolist()
        avg_counts["model_name"] = pd.Categorical(
            avg_counts["model_name"], categories=sorted_models, ordered=True
        )

        # Step 5: Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        sns.barplot(
            data=avg_counts,
            x="model_name",
            y="proportion",
            hue="agreement",
            hue_order=reversed_categories,
            order=sorted_models,
            ax=ax,
        )

        # Add proportion labels on top of each bar
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.text(
                    p.get_x() + p.get_width() / 2.0,
                    height + 0.01,
                    f"{height * 100:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        # Customize plot
        ax.set_xlabel("Model")
        ax.set_ylabel("Unweighted Average Proportion")
        ax.set_title(
            f"Unweighted Average Proportion of Agreement Categories by Model and {group_by.capitalize()}"
        )
        ax.set_ylim(0, 1.05)
        ax.legend(
            title="Agreement Categories", bbox_to_anchor=(1.05, 1), loc="upper left"
        )
        plt.xticks(rotation=90)
        plt.tight_layout()
        return fig, ax

    elif granularity == 2:
        # Granularity 2: Models within cell types and tissues
        if not sub_group_by:
            raise ValueError("sub_group_by must be provided for granularity=2.")

        # Step 1: Compute counts per model, group_by, sub_group_by, and agreement
        counts = (
            melted.groupby(["model_name", group_by, sub_group_by, "agreement"])
            .size()
            .reset_index(name="count")
        )

        # Step 2: Compute proportions within each model, group_by, and sub_group_by
        counts["proportion"] = counts.groupby(["model_name", group_by, sub_group_by])[
            "count"
        ].transform("sum")
        counts["proportion"] = counts["count"] / counts["proportion"]

        # Step 3: Compute unweighted average proportions across all (group_by, sub_group_by) combinations
        avg_counts = (
            counts.groupby(["model_name", "agreement"])["proportion"]
            .mean()
            .reset_index()
        )

        # Step 4: Sort models based on the average proportion of the highest agreement category
        sort_category = reversed_categories[0]
        sort_order = avg_counts[avg_counts["agreement"] == sort_category].set_index(
            "model_name"
        )["proportion"]
        sorted_models = sort_order.sort_values(ascending=False).index.tolist()
        avg_counts["model_name"] = pd.Categorical(
            avg_counts["model_name"], categories=sorted_models, ordered=True
        )

        # Step 5: Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        sns.barplot(
            data=avg_counts,
            x="model_name",
            y="proportion",
            hue="agreement",
            hue_order=reversed_categories,
            order=sorted_models,
            ax=ax,
        )

        # Add proportion labels on top of each bar
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.text(
                    p.get_x() + p.get_width() / 2.0,
                    height + 0.01,
                    f"{height * 100:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        # Customize plot
        ax.set_xlabel("Model")
        ax.set_ylabel("Unweighted Average Proportion")
        ax.set_title(
            f"Unweighted Average Proportion of Agreement Categories by Model, {group_by.capitalize()}, and {sub_group_by.capitalize()}"
        )
        ax.set_ylim(0, 1.05)
        ax.legend(
            title="Agreement Categories", bbox_to_anchor=(1.05, 1), loc="upper left"
        )
        plt.xticks(rotation=90)
        plt.tight_layout()
        return fig, ax

    else:
        raise ValueError("Granularity must be 0, 1, or 2.")
