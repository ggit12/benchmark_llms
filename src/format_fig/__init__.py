"""
init for format fig functions
"""

from .customize_figure import (
    customize_figure,
)

from .customize_clustermap import (
    customize_clustermap,
)

from .customize_scatterplot import (
    customize_scatterplot,
)

from .replace_tick_labels import (
    replace_tick_labels,
)

from .customize_legend_labels import (
    customize_legend_labels,
)

__all__ = [
    # customize_figure
    "customize_figure",

    # customize_clustermap
    "customize_clustermap",

    # customize_scatterplot
    "customize_scatterplot",

    # replace_tick_labels
    "replace_tick_labels",

    # customize_legend_labels
    "customize_legend_labels",

]
