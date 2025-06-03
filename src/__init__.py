"""
init file for src of benchmark_llms
"""

from . import providers
from . import run_multiple_providers
from . import extract_table_from_fig
from . import format_fig

from .providers import (
    PROVIDERS,
    ENDPOINTS,
    MODEL_TICK_LABELS,
    REMOVE_TICK_LABELS,
)

from .run_multiple_providers import (
    run_multiple_providers_models,
)

# pylint: disable=reimported
from .extract_table_from_fig import (
    extract_table_from_fig,
)

from .find_indices import (
    find_indices_closest_to_4_corners,
)

from .format_fig import (
    customize_figure,
    customize_clustermap,
    customize_scatterplot,
    replace_tick_labels,
    customize_legend_labels,
)

from .direct_compare_cell_type_labels_pairwise import (
    direct_compare_cell_type_labels_pairwise,
)

from .agreement_plots import (
    plot_model_agreement_unweighted,
    plot_model_agreement_categorical_unweighted
)

from .agreement_metrics import (
    cell_type_by_plurality,
    compute_agreement_df,
)

from .agreement_scatterplot import (
    calculate_weights,
    plot_agreement,
    plot_agreement_simple,
)

from .plot_genes_umap import (
    plot_genes_umap,
)

from .plot_kappa import (
    plot_pairwise_clustermap,
    plot_average_pairwise_barchart,
)

from .ensure_label_consistency_legacy import (
    ensure_label_consistency_adata,
    ensure_label_consistency_main,
)


__all__ = [
    # providers
    "PROVIDERS",
    "ENDPOINTS",
    "MODEL_TICK_LABELS",
    "REMOVE_TICK_LABELS",

    # run_multiple_providers
    "run_multiple_providers_models",

    # extract_table_from_fig
    "extract_table_from_fig",

    # find_indices
    "find_indices_closest_to_4_corners",

    # format_fig
    "customize_figure",
    "customize_clustermap",
    "customize_scatterplot",
    "replace_tick_labels",
    "customize_legend_labels",

    # direct_compare_cell_type_labels_pairwise
    "direct_compare_cell_type_labels_pairwise",

    # agreement_plots
    "plot_model_agreement_unweighted",
    "plot_model_agreement_categorical_unweighted",

    # agreement_metrics
    "cell_type_by_plurality",
    "compute_agreement_df",

    # agreement_scatterplot
    "calculate_weights",
    "plot_agreement",
    "plot_agreement_simple",

    # plot_genes_umap
    "plot_genes_umap",

    # plot_kappa
    "plot_pairwise_clustermap",
    "plot_average_pairwise_barchart",

    # ensure_label_consistency_legacy
    "ensure_label_consistency_adata",
    "ensure_label_consistency_main",
]
