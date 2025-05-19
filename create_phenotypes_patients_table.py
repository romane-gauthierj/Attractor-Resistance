import matplotlib.pyplot as plt
import pandas as pd

# Assuming you already have:
# df_resistant_mean and df_sensitive_mean as described
# We'll select only the 3 phenotypes: <nil>, Apoptosis, Proliferation



def vizualise_table_phenotype_condition(folder, patient_resistant_mean, patient_sensitive_mean):
    patient_resistant_mean = patient_resistant_mean.round(2)
    patient_sensitive_mean = patient_sensitive_mean.round(2)
    patient_resistant_mean = patient_resistant_mean.set_index(patient_resistant_mean.columns[0])
    patient_sensitive_mean = patient_sensitive_mean.set_index(patient_sensitive_mean.columns[0])
    

    outcomes = ['<nil>', 'Apoptosis', 'Proliferation', 'Metastasis']
    conditions = list(patient_resistant_mean.index)

    # Build full header with two levels
    top_header = []
    sub_header = []

    for outcome in outcomes:
        top_header.extend([outcome, ''])  # Span 2 columns per outcome
        sub_header.extend(['Resistant', 'Sensitive'])

    # Build the table content (rows = conditions)
    cell_text = []
    for condition in conditions:
        row = []
        for outcome in outcomes:
            val_res = patient_resistant_mean.at[condition, outcome]
            val_sens = patient_sensitive_mean.at[condition, outcome]
            row.extend([val_res, val_sens])
        cell_text.append(row)

    # Add the sub-header row
    cell_text.insert(0, sub_header)

    df = pd.DataFrame(cell_text[1:], columns=top_header)
    #df.to_excel('phenotype_comparison.xlsx', index=True)

    row_labels = [''] + conditions  # Empty label for the sub-header row

    # Plotting
    fig, ax = plt.subplots(figsize=(18, 6))
    ax.axis('off')

    # Add table
    table = ax.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=top_header,
        colWidths=[0.06] * len(top_header),
        cellLoc='center',
        loc='center'
    )

    # Styling
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 2.2)

    # Optional: bold top row
    for col in range(len(top_header)):
        cell = table[0, col]
        cell.set_fontsize(9)
        cell.set_text_props(weight='bold')

    # Optional: shade the top two rows
    for col in range(len(top_header)):
        table[0, col].set_facecolor('#add8e6')
        table[1, col].set_facecolor('#d3e0ea')


    plt.title("Cell Fate Outcomes by Input Condition\nGrouped by Phenotype with Resistant/Sensitive Comparison", pad=14)
    plt.tight_layout()
    
    output_path = f'{folder}/sensitive_resistant_results/table_expression_per_phenotype.png' 
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    df.to_csv(f'{folder}/sensitive_resistant_results/table_expression_per_phenotype.csv', index=False)

    plt.show()

