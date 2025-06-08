"""
This script can be used to calculate the correlations between the 
performance metrics of separate aggregate runs.
Used to produce Supplementary Tables 4 and 6.
"""
# pylint: disable=line-too-long, invalid-name

import pickle
import pandas as pd
from scipy import stats

#define paths
base_path = "/Path/to/data"

agg_metrics_1_path = f"{base_path}/first_run/res/10_aggregated_tables/aggregated_performance_table_mean.pkl"
agg_metrics_2_path = f"{base_path}/second_run/res/10_aggregated_tables/aggregated_performance_table_mean.pkl"


# Load the means from each run
mean1 = pickle.load(open(agg_metrics_1_path, "rb"))
mean2 = pickle.load(open(agg_metrics_2_path, "rb"))

#Calculate Pearson correlation coefficients between the mean performance metrics of the two aggregated results
corr_results = pd.DataFrame(index=mean1.columns, columns=["r", "p"])

for col in corr_results.index:
    r, p = stats.pearsonr(mean1[col], mean2[col])
    corr_results.loc[col] = [r, p]

# Save
corr_results.to_pickle("correlations_of_performance_metrics.pkl")

# Convert the p values to scientific notation
corr_results["p"] = corr_results["p"].apply(lambda x: f"{x:.1e}")
# Round the r values to 2 decimal places
corr_results["r"] = corr_results["r"].apply(lambda x: f"{x:.2f}")

# rename the columns to indicate that these are postprocessed by claude
corr_results.columns = ["Pearson's R", "P"]

corr_results.to_html("correlations_of_performance_metrics.html")
