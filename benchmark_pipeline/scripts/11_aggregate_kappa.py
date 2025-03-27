"""
This rule aggregates kappa across runs and plots.
"""

import sys
import os
import argparse
import pickle
import numpy as np
import matplotlib


from dotenv import load_dotenv

load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    plot_pairwise_clustermap,
    customize_clustermap,
)

from src import MODEL_TICK_LABELS as model_tick_labels

# Use a non-GUI backend
matplotlib.use('Agg')


def load_pickle(file_path):
    """Load a pickle file."""
    with open(file_path, "rb") as f:
        return pickle.load(f)


def save_pickle(obj, file_path):
    """Save an object to a pickle file."""
    with open(file_path, "wb") as f:
        pickle.dump(obj, f)


def aggregate_kappa_objects(kappa_files):
    """
    Aggregate kappa objects from multiple runs into mean and standard deviation.

    Parameters:
    -----------
    kappa_files : list
        List of file paths to kappa pickle files.

    Returns:
    --------
    tuple
        (kappa_mean, kappa_sd) dictionaries with same structure as input kappa objects.
    """
    print(f"Aggregating {len(kappa_files)} kappa files...")

    # Load all kappa objects
    kappas = [load_pickle(f) for f in kappa_files]

    # Initialize result dictionaries
    kappa_mean = {}
    kappa_sd = {}

    # Process each key in the kappa objects
    for key in kappas[0].keys():
        if key == "pairwise":
            # For pairwise data, which is a nested dictionary
            kappa_mean[key] = {}
            kappa_sd[key] = {}

            # Get all possible pairs from the first kappa object
            all_pairs = kappas[0][key].keys()

            for pair in all_pairs:
                # Extract values for this pair from all kappa objects
                values = [k[key][pair] for k in kappas if pair in k[key]]
                if values:
                    kappa_mean[key][pair] = np.mean(values)
                    kappa_sd[key][pair] = np.std(values)
        elif key == "average_pairwise":
            # For average_pairwise data, also a dictionary
            kappa_mean[key] = {}
            kappa_sd[key] = {}

            # Get all possible raters from the first kappa object
            all_raters = kappas[0][key].keys()

            for rater in all_raters:
                # Extract values for this rater from all kappa objects
                values = [k[key][rater] for k in kappas if rater in k[key]]
                if values:
                    kappa_mean[key][rater] = np.mean(values)
                    kappa_sd[key][rater] = np.std(values)
        else:
            # For scalar values like 'fleiss'
            values = [k[key] for k in kappas]
            kappa_mean[key] = np.mean(values)
            kappa_sd[key] = np.std(values)

    return kappa_mean, kappa_sd


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate kappa objects from multiple runs."
    )
    parser.add_argument(
        "--input-kappas", nargs="+", required=True, help="Input kappa pickle files"
    )
    args = parser.parse_args()

    # Process kappa files
    kappa_mean, kappa_sd = aggregate_kappa_objects(args.input_kappas)

    # Save aggregated kappa objects
    save_pickle(kappa_mean, "res/11_aggregated_kappa/kappa_mean.pkl")
    save_pickle(kappa_sd, "res/11_aggregated_kappa/kappa_sd.pkl")

    # Plot kappa_mean
    kappa_clustermap_mean = plot_pairwise_clustermap(
        kappa_mean, metric="cosine", method="centroid"
    )
    kappa_clustermap_mean = customize_clustermap(
        kappa_clustermap_mean,
        remove_legend=True,
        x_tick_substrings=["consistent_including_manual_", "_simplified_ai_cell_type"],
        y_tick_substrings=["consistent_including_manual_", "_simplified_ai_cell_type"],
        new_ylabel="",
        fig_width=3,
        fig_height=3.7,
        remove_value_labels=True,
        new_tick_labels=model_tick_labels,
    )
    kappa_clustermap_mean[0].savefig(
        "res/11_aggregated_kappa/kappa_mean.svg", format="svg"
    )

    # Plot kappa_sd
    kappa_clustermap_sd = plot_pairwise_clustermap(
        kappa_sd, metric="cosine", method="centroid"
    )
    kappa_clustermap_sd = customize_clustermap(
        kappa_clustermap_sd,
        remove_legend=True,
        x_tick_substrings=["consistent_including_manual_", "_simplified_ai_cell_type"],
        y_tick_substrings=["consistent_including_manual_", "_simplified_ai_cell_type"],
        new_ylabel="",
        fig_width=3,
        fig_height=3.7,
        remove_value_labels=True,
        new_tick_labels=model_tick_labels,
    )
    kappa_clustermap_sd[0].savefig("res/11_aggregated_kappa/kappa_sd.svg", format="svg")

if __name__ == "__main__":
    main()
