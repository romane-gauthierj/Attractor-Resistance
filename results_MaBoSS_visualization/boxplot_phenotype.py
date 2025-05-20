# Create Boxplots of phenotype distribution


import pandas as pd
import ast
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



# Step 1: Convert the dataframes to long format and min-max normalization
def melt_df(df, group_label):
    records = []
    phenotypes = df.index
   
    conditions = df.columns
    for phenotype in phenotypes:
        all_values = []
        for condition in conditions:
            values = ast.literal_eval(df.loc[phenotype, condition])
            all_values.extend(values)
        # Min-max scaling
        min_val, max_val = min(all_values), max(all_values)
        normalized_values = [(v - min_val) / (max_val - min_val + 1e-9) for v in all_values]  # Avoid division by zero

        i = 0
        for condition in conditions:
            values = ast.literal_eval(df.loc[phenotype, condition])
            for _ in values:
                records.append({
                    'Phenotype': phenotype,
                    'Condition': condition,
                    'Group': group_label,
                    'Expression': normalized_values[i],
                    'Condition_Group': f'{condition} ({group_label})'
                })
                i += 1
    return pd.DataFrame(records)



def create_boxplot(folder, patient_res_values, patient_sens_values, data_two_sides):

    patient_res_values.set_index('Unnamed: 0', inplace=True)
    patient_sens_values.set_index('Unnamed: 0', inplace=True)
    df_res_long = melt_df(patient_res_values, 'Resistant')
    df_sens_long = melt_df(patient_sens_values, 'Sensitive')

    df_combined = pd.concat([df_res_long, df_sens_long], ignore_index=True)
    df_combined = df_combined.loc[:, ~df_combined.columns.duplicated()]
    df_combined['Condition_Group'] = df_combined['Condition'] + ' (' + df_combined['Group'] + ')'


    significant_pairs = {
    (row['Phenotype'], row['Condition']): row['Star_significant']
    for _, row in data_two_sides.iterrows()
    if row['Star_significant'] != ''
}
    
# if only want to keep the phenotype-condition significant 

#     df_combined = df_combined[
#     df_combined.apply(lambda row: (row['Phenotype'], row['Condition']) in significant_pairs, axis=1)
# ]
    print('the df_combined is:')
    print(df_combined)
    print(df_combined['Phenotype'].unique())

    # Step 2: Plot — one figure per phenotype
    g = sns.catplot(
        data=df_combined,
        kind='box',
        x='Condition',
        y='Expression',
        hue='Group',            # Separates Resistant vs Sensitive
        col='Phenotype',        # One subplot (facet) per phenotype
        col_wrap=4,             # Adjust this value depending on the number of phenotypes
        sharey=False,           # Each subplot gets its own y-axis range if needed
        palette={'Resistant': '#FF7F0E', 'Sensitive': '#008000'},
        height=6,               # Height of each subplot
        aspect=0.8              # Aspect ratio of each subplot (width = height * aspect)
    )

# Step 3: Apply significance annotation to each subplot
    for ax in g.axes.flat:
             
        # Get the current phenotype from the subplot title
        phenotype = ax.get_title().split(' = ')[-1] if ' = ' in ax.get_title() else ax.get_title()


# comment here
        for (sig_phenotype, condition), stars in significant_pairs.items():
            if sig_phenotype == phenotype:
                # Get the data for just this condition and phenotype
                subset = df_combined[(df_combined['Phenotype'] == phenotype) & (df_combined['Condition'] == condition)]
                
                if not subset.empty:
                    max_expression = subset['Expression'].max()
                    
                    # Find the tick positions and labels
                    tick_labels = [t.get_text() for t in ax.get_xticklabels()]
                    try:
                        xpos = tick_labels.index(condition)
                        ax.text(
                        xpos,
                        max_expression + 0.01,  # much closer to the box
                        stars,
                        ha='center',
                        va='bottom',
                        fontsize=16,
                        fontweight='bold',
                        color='red'
                        )
                    except ValueError:
                        # Condition not found among x-tick labels for this subplot
                        continue

    # Improve formatting
    g.set_xticklabels(rotation=45, ha='right', fontsize=10)
    g.set_axis_labels("", "Expression") 

    g.set_titles("{col_name}", size=14, fontweight='bold')
    g.fig.subplots_adjust(top=0.9)

    # Move the legend to the top-right corner

    g.legend.set_bbox_to_anchor((0.6, 0.9))  # Position the legend outside of the plot

    # Set the font size for the legend
    for legend_text in g.legend.get_texts():
        legend_text.set_fontsize(10)

    # Add the main title for the entire figure
    #g.fig.suptitle("Boxplot of Expression Levels by Condition and Group\nfor Each Phenotype", fontsize=16)

    # Tighten the layout to avoid overlap
    plt.tight_layout()

    output_path = f'{folder}/sensitive_resistant_results/boxplot_expression_per_phenotype.png'  # or use .pdf/.svg
    g.legend.set_bbox_to_anchor((0.6, 0.9))
    g.legend.get_frame().set_edgecolor('black')
    g.legend.get_frame().set_linewidth(1)
    g.legend.get_frame().set_boxstyle('round')
    g.legend.set_title("Group", prop={'weight': 'bold', 'size': 12})

    for legend_text in g.legend.get_texts():
        legend_text.set_fontsize(10)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    # Show the plot (optional if you're running this headlessly)
    plt.show()




















