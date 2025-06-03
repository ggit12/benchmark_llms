"""
This script aggregates performance tables across multiple runs.
"""

import argparse
import pickle

import pandas as pd
import sigfig


def format_with_sd(mean, sd, make_percent=False):
    """Format mean and standard deviation as 'mean (sd)' with appropriate rounding."""
    if make_percent:
        mean = mean * 100
        sd = sd * 100

    # Round mean to same decimal place as sd
    mean_rounded = sigfig.round(mean, uncertainty=sd, format='PDG')

    return mean_rounded


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Aggregate performance tables across multiple runs."
    )
    parser.add_argument(
        "--input-tables",
        nargs="+",
        required=True,
        help="List of input table pickle files",
    )
    args = parser.parse_args()

    # Load all tables
    tables = []
    for table_path in args.input_tables:
        with open(table_path, "rb") as f:
            table = pickle.load(f)
            tables.append(table)

    # Calculate mean and standard deviation across all tables
    df_concat = pd.concat(tables)
    mean_table = df_concat.groupby(df_concat.index).mean()
    std_table = df_concat.groupby(df_concat.index).std()

    # Write the mean and standard deviation tables to pickle files
    with open(
        "./res/10_aggregated_tables/aggregated_performance_table_mean.pkl", "wb"
    ) as file:
        pickle.dump(mean_table, file)
    with open(
        "./res/10_aggregated_tables/aggregated_performance_table_std.pkl", "wb"
    ) as file:
        pickle.dump(std_table, file)

    # Create the aggregated table
    aggregated_table = mean_table.copy()

    # Define columns to format
    format_as_percent = [
        "Overall Binary (% of Cells)",
        "Overall Binary (% of Cell Types)",
        "Perfect Match (% of Cells)",
        "Perfect Match (% of Cell Types)",
        "Exact String Match (% of Cells)",
        "Exact String Match (% of Cell Types)"
    ]

    format_as_rounded = [
        "Kappa with Manual Annotations",
        "Average Kappa with Other LLMs",
    ]

    # Format each cell as "mean (sd)" with appropriate rounding
    for col in aggregated_table.columns:
        if col in format_as_percent:
            aggregated_table[col] = [
                format_with_sd(
                    mean_table.loc[idx, col], std_table.loc[idx, col], make_percent=True
                )
                for idx in mean_table.index
            ]
        elif col in format_as_rounded:
            aggregated_table[col] = [
                format_with_sd(
                    mean_table.loc[idx, col], std_table.loc[idx, col], make_percent=False
                )
                for idx in mean_table.index
            ]

    # Drop extra columns
    for col in aggregated_table.columns:
        if col not in format_as_percent + format_as_rounded:
            aggregated_table.drop(col, axis=1, inplace=True)

    # Reset the index to make 'Model' a column
    aggregated_table.reset_index(inplace=True)

    # Convert the DataFrame to HTML without the index and header
    html_table = aggregated_table.to_html(
        index=False,
        header=False,
        classes="sortable performance-table",
        border=0,
        escape=False,
    )

    # Strip the unnecessary <table> and </table> tags from the output, leaving just the rows
    html_table = html_table.split("<tbody>")[1].split("</tbody>")[0]

    # Manually add the custom headers and wrap it in a proper table structure
    html_content = f"""
<table class="sortable performance-table" id="leaderboard">
  <thead>
    <tr>
      <th data-column="0" data-numeric="false"></th>
      <th data-column="1" data-numeric="true">Overall Binary (% of Cells)</th>
      <th data-column="2" data-numeric="true">Overall Binary (% of Cell Types)</th>
      <th data-column="3" data-numeric="true">Perfect Match (% of Cells)</th>
      <th data-column="4" data-numeric="true">Perfect Match (% of Cell Types)</th>
      <th data-column="5" data-numeric="true">Exact String Match (% of Cells)</th>
      <th data-column="6" data-numeric="true">Exact String Match (% of Cell Types)</th>
      <th data-column="7" data-numeric="true">Kappa with Manual Annotations</th>
      <th data-column="8" data-numeric="true">Average Kappa with Other LLMs</th>
    </tr>
  </thead>
  <tbody>
    {html_table}
  </tbody>
</table>
"""

    # Write the modified HTML content to a file
    with open(
        "./res/10_aggregated_tables/aggregated_performance_table.html",
        "w",
        encoding="utf-8",
    ) as file:
        file.write(html_content)

    print("Aggregated table saved to HTML file.")


if __name__ == "__main__":
    main()
