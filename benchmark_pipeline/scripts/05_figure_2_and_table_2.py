"""
This script generates results presented in table 2. Also generates the information as figures (figures not used).

Files saved:
1. './res/05_figure_2_and_table_2/agreement_plot_overall_binary.svg'
2. './res/05_figure_2_and_table_2/agreement_plot_overall_binary_unweighted.svg'
3. './res/05_figure_2_and_table_2/agreement_plot_overall_categorical_perfect_only.svg'
4. './res/05_figure_2_and_table_2/agreement_plot_overall_categorical.svg'
5. './res/05_figure_2_and_table_2/agreement_plot_overall_categorical_perfect_only_unweighted.svg'
6. './res/05_figure_2_and_table_2/agreement_plot_overall_categorical_unweighted.svg'
7. './res/05_figure_2_and_table_2/kappa_clustermap_overall.svg'
8. './res/05_figure_2_and_table_2/average_pairwise_kappa.svg'
9. './res/05_figure_2_and_table_2/performance_table.html'

Note: "unweighted" means that the cell types have not been weighted by their size (number of cells), 
while "weighted" means that the metric reflects the proportion of cells in total (not cell types). 
The "weighted" suffix is dropped in the code below, 
only "unweighted" is used (so the absence of "unweighted" means "weighted").
"""
# pylint: disable=line-too-long
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

import pickle
import os
import sys

import anndict as adt
import scanpy as sc
import pandas as pd
import matplotlib

from dotenv import load_dotenv
load_dotenv()

source_dir = os.environ["SOURCE_DIR"]  # retrieve the path to src from env var
sys.path.append(os.path.join(source_dir))  # add to Python path

# pylint: disable=wrong-import-position
from src import (
    customize_figure,
    customize_clustermap,
    plot_pairwise_clustermap,
    plot_average_pairwise_barchart,
    extract_table_from_fig,
    plot_model_agreement_unweighted,
    plot_model_agreement_categorical_unweighted,
)

from src import MODEL_TICK_LABELS as model_tick_labels

# Use a non-GUI backend
matplotlib.use('Agg')


# Read the results
adata = sc.read_h5ad('./res/04_postprocess_results/adt_de_novo_llm_annotated.h5ad')

# And the various column names
manual_cell_type_col = pickle.load(open('../../dat/manual_cell_type_col.pkl', 'rb'))
llm_celltype_cols = pickle.load(open('./res/04_postprocess_results/llm_celltype_cols.pkl', 'rb'))
binary_agreement_cols = pickle.load(open('./res/04_postprocess_results/binary_agreement_cols.pkl', 'rb'))
categorical_agreement_cols = pickle.load(open('./res/04_postprocess_results/categorical_agreement_cols.pkl', 'rb'))
perfect_only_categorical_agreement_cols = pickle.load(open('./res/04_postprocess_results/perfect_only_categorical_agreement_cols.pkl', 'rb'))
direct_string_agreement_cols = pickle.load(open('./res/04_postprocess_results/direct_string_agreement_cols.pkl', 'rb'))

#plot model binary agreement overall
agreement_plot_overall_binary = adt.plot_model_agreement(adata, group_by=manual_cell_type_col, sub_group_by='tissue', agreement_cols=binary_agreement_cols, granularity=0)

# Get the top n models based on the overall binary agreement, not including the plurality vote column
# Extract from the plot object
top_n = 4
exclude_model = 'binary_agreement_consistent_including_manual_' + manual_cell_type_col + '_cell_type_by_plurality'
xtick_labels = [tick.get_text() for tick in agreement_plot_overall_binary[1].get_xticklabels()]
bar_heights = [bar.get_height() for bar in agreement_plot_overall_binary[1].patches]
df = pd.DataFrame({'Model': xtick_labels, 'agreement': bar_heights}).set_index('Model').drop(exclude_model, errors='ignore')
binary_agreement_cols_top_models = df.sort_values(by='agreement', ascending=False).head(top_n).index.tolist()

# Also get the corresponding col names for the top models from categorical_agreement_cols, perfect_only_categorical_agreement_cols, direct_string_agreement_cols, and llm_celltype_cols_top_models
categorical_agreement_cols_top_models = ['categorical_agreement_' + col.replace('binary_agreement_', '') for col in binary_agreement_cols_top_models]
# This can be done by stripping the 'binary_agreement_' prefix from the model names and replacing it with 'perfect_only_categorical_agreement_'
perfect_only_categorical_agreement_cols_top_models = ['perfect_only_categorical_agreement_' + col.replace('binary_agreement_', '') for col in binary_agreement_cols_top_models]
# And by stripping the 'binary_agreement_' prefix from the model names and replacing it with 'direct_string_agreement_'
direct_string_agreement_cols_top_models = ['direct_string_agreement_' + col.replace('binary_agreement_', '') for col in binary_agreement_cols_top_models]
# And by stripping the 'binary_agreement_consistent_including_manual_' + manual_cell_type_col prefix from the model names
llm_celltype_cols_top_models = [col.replace('binary_agreement_consistent_including_manual_' + manual_cell_type_col + '_', '') for col in binary_agreement_cols_top_models]

# Write these model cols for later use
pickle.dump(binary_agreement_cols_top_models, open('./res/05_figure_2_and_table_2/binary_agreement_cols_top_models.pkl', 'wb'))
pickle.dump(categorical_agreement_cols_top_models, open('./res/05_figure_2_and_table_2/categorical_agreement_cols_top_models.pkl', 'wb'))
pickle.dump(perfect_only_categorical_agreement_cols_top_models, open('./res/05_figure_2_and_table_2/perfect_only_categorical_agreement_cols_top_models.pkl', 'wb'))
pickle.dump(direct_string_agreement_cols_top_models, open('./res/05_figure_2_and_table_2/direct_string_agreement_cols_top_models.pkl', 'wb'))
pickle.dump(llm_celltype_cols_top_models, open('./res/05_figure_2_and_table_2/llm_celltype_cols_top_models.pkl', 'wb'))


# Get the tick labels and define custom ones
# [tick.get_text() for tick in agreement_plot_overall_binary[1].get_xticklabels()]

# Customize the figure
customize_figure(agreement_plot_overall_binary, remove_legend = True, x_tick_substrings=['binary_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (yes/no)',
                 new_tick_labels=model_tick_labels)

 #Save the plot as an SVG file
agreement_plot_overall_binary[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_binary.svg', format='svg')

agreement_table_overall_binary = extract_table_from_fig(agreement_plot_overall_binary, value_col_name='Overall Binary (% of Cells)')


# Same plot as above but unweighted
agreement_plot_overall_binary_unweighted = plot_model_agreement_unweighted(adata, group_by=manual_cell_type_col, sub_group_by='tissue', model_cols=binary_agreement_cols, granularity=0)
customize_figure(agreement_plot_overall_binary_unweighted, remove_legend = True, x_tick_substrings=['binary_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (yes/no)', new_tick_labels=model_tick_labels)
agreement_plot_overall_binary_unweighted[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_binary_unweighted.svg', format='svg')

# Get info as a table
agreement_table_overall_binary_unweighted = extract_table_from_fig(agreement_plot_overall_binary_unweighted, value_col_name="Overall Binary (% of Cell Types)")


#plot the model categorical agreement overall (perfect only)
agreement_plot_overall_categorical_perfect = adt.plot_model_agreement(adata, group_by=manual_cell_type_col, sub_group_by='tissue', agreement_cols=perfect_only_categorical_agreement_cols, granularity=0)
customize_figure(agreement_plot_overall_categorical_perfect, remove_legend = True, x_tick_substrings=['perfect_only_categorical_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (% perfect matches)', new_tick_labels=model_tick_labels)
agreement_plot_overall_categorical_perfect[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_categorical_perfect_only.svg', format='svg')

agreement_table_overall_categorical_perfect = extract_table_from_fig(agreement_plot_overall_categorical_perfect, value_col_name="Perfect Match (% of Cells)")

#and separated by level of agreement
#This is Figure 2A
agreement_plot_overall_categorical = adt.plot_model_agreement_categorical(adata, group_by=manual_cell_type_col, sub_group_by='tissue', agreement_cols=categorical_agreement_cols, granularity=0)
customize_figure(agreement_plot_overall_categorical, remove_legend = True,
                 x_tick_substrings=['categorical_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'],
                 new_ylabel='Agreement with Manual Annotation (by level of agreement)', remove_bar_labels=True,
                 fig_width=2.4,
                 fig_height=3,
                 new_tick_labels=model_tick_labels)
agreement_plot_overall_categorical[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_categorical.svg', format='svg')

agreement_table_overall_categorical = extract_table_from_fig(agreement_plot_overall_categorical, value_col_name="Categorical Agreement (% of Cells)")

#unweighted versions of the plots
#(perfect only)
agreement_plot_overall_categorical_perfect_unweighted = plot_model_agreement_unweighted(adata, group_by=manual_cell_type_col, sub_group_by='tissue', model_cols=perfect_only_categorical_agreement_cols, granularity=0)

#and separated by level of agreement
#This is Figure 2B
agreement_plot_overall_categorical_unweighted = plot_model_agreement_categorical_unweighted(adata, group_by=manual_cell_type_col, sub_group_by='tissue', model_cols=categorical_agreement_cols, granularity=0)

customize_figure(agreement_plot_overall_categorical_perfect_unweighted, remove_legend = True, x_tick_substrings=['perfect_only_categorical_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (% perfect matches)', new_tick_labels=model_tick_labels)
customize_figure(agreement_plot_overall_categorical_unweighted, remove_legend = True, x_tick_substrings=['categorical_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (by level of agreement)', remove_bar_labels=True, fig_width=2.4, fig_height=3, new_tick_labels=model_tick_labels)

agreement_plot_overall_categorical_perfect_unweighted[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_categorical_perfect_only_unweighted.svg', format='svg')
agreement_plot_overall_categorical_unweighted[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_categorical_unweighted.svg', format='svg')

# Extract Table
agreement_table_categorical_perfect_unweighted = extract_table_from_fig(agreement_plot_overall_categorical_perfect_unweighted, value_col_name="Perfect Match (% of Cell Types)")
agreement_table_categorical_unweighted = extract_table_from_fig(agreement_plot_overall_categorical_unweighted, value_col_name="Categorical Agreement (% of Cell Types)")

#Same plots as above but for direct string comparison
agreement_plot_overall_direct = adt.plot_model_agreement(adata, group_by=manual_cell_type_col, sub_group_by='tissue', agreement_cols=direct_string_agreement_cols, granularity=0)
agreement_plot_overall_direct_unweighted = plot_model_agreement_unweighted(adata, group_by=manual_cell_type_col, sub_group_by='tissue', model_cols=direct_string_agreement_cols, granularity=0)

customize_figure(agreement_plot_overall_direct, remove_legend = True, x_tick_substrings=['direct_string_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (% exact string match)', new_tick_labels=model_tick_labels)
customize_figure(agreement_plot_overall_direct_unweighted, remove_legend = True, x_tick_substrings=['direct_string_agreement_consistent_including_manual_' + manual_cell_type_col + '_', 'consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Agreement with Manual Annotation (% exact string match)', new_tick_labels=model_tick_labels)

agreement_plot_overall_direct[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_direct_string.svg', format='svg')
agreement_plot_overall_direct_unweighted[0].savefig('./res/05_figure_2_and_table_2/agreement_plot_overall_direct_string_unweighted.svg', format='svg')

# Extract Table
agreement_table_overall_direct = extract_table_from_fig(agreement_plot_overall_direct, value_col_name="Exact String Match (% of Cells)")
agreement_table_overall_direct_unweighted = extract_table_from_fig(agreement_plot_overall_direct_unweighted, value_col_name="Exact String Match (% of Cell Types)")

#calculate kappa
kappa = adt.kappa_adata(adata, llm_celltype_cols)

# save for later aggregation
pickle.dump(kappa, open('./res/05_figure_2_and_table_2/kappa.pkl', 'wb'))

# calculate kappa including manual column
consistent_manual_cell_type_col = 'consistent_including_manual_' + manual_cell_type_col
kappa_with_manual = adt.kappa_adata(adata, llm_celltype_cols + [consistent_manual_cell_type_col])

# calculate alpha
# alpha = adt.krippendorff_alpha_adata(adata, llm_celltype_cols)

# For kappa values
def extract_kappa_pairs(data_dict, column_name, substring_to_remove, replace_dict):
    """Extracts the kappa values from the pairwise dictionary and returns them as a DataFrame. """
    pairwise_data = data_dict['pairwise']
    extracted_data = []

    # Iterate through the pairwise dictionary
    for (model1, model2), kappa in pairwise_data.items():
        # Check if the column_name is in either model1 or model2
        if column_name in model1:
            other_model = model2
        elif column_name in model2:
            other_model = model1
        else:
            continue

        # Append the other model and kappa value
        extracted_data.append([other_model, kappa])

    # Convert to DataFrame with the specified column names
    df = pd.DataFrame(extracted_data, columns=['Model', 'Kappa with Manual Annotations']).drop_duplicates(subset='Model', keep='first')
    df['Model'] = df['Model'].replace(substring_to_remove, regex=True)

    # Replace values in 'Model' column based on the replace_dict argument
    df['Model'] = df['Model'].replace(replace_dict)

    # Round non-'Model' column (Kappa) to 3 decimal places
    df['Kappa with Manual Annotations'] = df['Kappa with Manual Annotations'].round(3)

    # Set 'Model' as the index
    df.set_index('Model', inplace=True)

    return df

kappa_with_manual_df = extract_kappa_pairs(kappa_with_manual, consistent_manual_cell_type_col, substring_to_remove={'consistent_including_manual_':'', '_simplified_ai_cell_type':''},
                                           replace_dict=model_tick_labels)


# Figure 2C:
# plot kappa
kappa_clustermap = plot_pairwise_clustermap(kappa, metric='cosine', method='centroid')

kappa_clustermap = customize_clustermap(kappa_clustermap, remove_legend=True, x_tick_substrings=['consistent_including_manual_', '_simplified_ai_cell_type'],
                                        y_tick_substrings=['consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='', fig_width=3, fig_height=3.7, remove_value_labels=True,
                                        new_tick_labels=model_tick_labels)

kappa_clustermap[0].savefig('./res/05_figure_2_and_table_2/kappa_clustermap_overall.svg', format='svg')

# same info as table, but as a bar plot
average_pairwise_kappa = plot_average_pairwise_barchart(kappa)
customize_figure(average_pairwise_kappa, remove_legend = True, x_tick_substrings=['consistent_including_manual_', '_simplified_ai_cell_type'], new_ylabel='Average Kappa', new_tick_labels=model_tick_labels)

#reformat plot and save
average_pairwise_kappa[0].savefig('./res/05_figure_2_and_table_2/average_pairwise_kappa.svg', format='svg')

# Extract Table
kappa_table = extract_table_from_fig(average_pairwise_kappa, value_col_name="Average Kappa with Other LLMs")

#write this for debugging
kappa_table.to_pickle('./res/05_figure_2_and_table_2/kappa_table.pkl')

#merge all the value tables
performance_table = pd.concat([
                                agreement_table_overall_binary, agreement_table_overall_binary_unweighted,
                                agreement_table_overall_categorical_perfect, agreement_table_categorical_perfect_unweighted,
                                agreement_table_overall_direct, agreement_table_overall_direct_unweighted,
                                agreement_table_overall_categorical, agreement_table_categorical_unweighted,
                                kappa_with_manual_df, kappa_table
                              ], axis=1)

# Write the performance table as .pkl for later aggregation
performance_table.to_pickle('./res/05_figure_2_and_table_2/performance_table.pkl')

# Create a copy of performance_table
performance_table_for_html = performance_table.copy()

# Format the values in the table as percentages and/or round
format_as_percent = ['Overall Binary (% of Cells)', 'Overall Binary (% of Cell Types)',
                    'Perfect Match (% of Cells)', 'Perfect Match (% of Cell Types)',
                    'Exact String Match (% of Cells)', 'Exact String Match (% of Cell Types)'
                    ]

format_as_rounded = ['Kappa with Manual Annotations', 'Average Kappa with Other LLMs']

for col in format_as_percent:
    performance_table_for_html[col] = performance_table_for_html[col].apply(lambda x: f"{x:.2%}")

for col in format_as_rounded:
    performance_table_for_html[col] = performance_table_for_html[col].apply(lambda x: f"{x:.3f}")

# Drop extra columns
for col in performance_table_for_html.columns:
    if col not in format_as_percent + format_as_rounded:
        performance_table_for_html.drop(col, axis=1, inplace=True)

# Table 2 is performance_table_for_html (its value at this point in the code). Below, we write it to an HTML file.

# Reset the index on the copy to make 'Model' a column
performance_table_for_html.reset_index(inplace=True)

# Convert the copied DataFrame to HTML without the index and header
# This gives us only the table body (rows) part of the table
html_table = performance_table_for_html.to_html(index=False, header=False, classes='sortable performance-table', border=0, escape=False)

# Strip the unnecessary <table> and </table> tags from the output, leaving just the rows
html_table = html_table.split('<tbody>')[1].split('</tbody>')[0]

# Manually add the custom headers and wrap it in a proper table structure
html_content = f'''
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
'''

# Write the modified HTML content to a file
with open("./res/05_figure_2_and_table_2/performance_table.html", "w", encoding="utf-8") as file:
    file.write(html_content)
