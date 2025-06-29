import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os


def vizualise_table_phenotype_condition(
    folder_results, patient_resistant_mean, patient_sensitive_mean
):
    """Visualizes a table comparing resistant and sensitive phenotypes across different conditions.
    Args:
        folder (str): The folder where the results will be saved.
        patient_resistant_mean (pd.DataFrame): DataFrame containing mean values for resistant phenotypes
        patient_sensitive_mean (pd.DataFrame): DataFrame containing mean values for sensitive phenotypes
    """

    conditions = list(patient_resistant_mean.index)
    phenotypes = list(patient_resistant_mean.columns)
    # Build full header with two levels
    top_header = []
    sub_header = []

    for phenotype in phenotypes:
        top_header.extend([phenotype, ""])  # Span 2 columns per phenotype
        sub_header.extend(["Resistant", "Sensitive"])

    # Build the table content (rows = conditions)
    cell_text = []
    for condition in conditions:
        row = []
        for phenotype in phenotypes:
            val_res = round(patient_resistant_mean.at[condition, phenotype], 2)
            val_sens = round(patient_sensitive_mean.at[condition, phenotype], 2)
            row.extend([val_res, val_sens])
        cell_text.append(row)

    # round the values in the DataFrame
    patient_resistant_mean = patient_resistant_mean.round(2)
    patient_sensitive_mean = patient_sensitive_mean.round(2)

    # Add the sub-header row
    cell_text.insert(0, sub_header)

    df = pd.DataFrame(cell_text[1:], columns=top_header)
    # df.to_excel('phenotype_comparison.xlsx', index=True)

    row_labels = [""] + conditions  # Empty label for the sub-header row

    # Plotting
    fig, ax = plt.subplots(figsize=(18, 6))
    ax.axis("off")

    # Add table
    table = ax.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=top_header,
        colWidths=[0.06] * len(top_header),
        cellLoc="center",
        loc="center",
    )

    # Styling
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 2.2)

    # Optional: bold top row
    for col in range(len(top_header)):
        cell = table[0, col]
        cell.set_fontsize(9)
        cell.set_text_props(weight="bold")

    # Optional: shade the top two rows
    for col in range(len(top_header)):
        table[0, col].set_facecolor("#add8e6")
        table[1, col].set_facecolor("#d3e0ea")

    plt.title(
        "Cell Fate Outcomes by Input Condition\nGrouped by Phenotype with Resistant/Sensitive Comparison",
        pad=14,
    )
    plt.tight_layout()

    output_path = f"{folder_results}/output"
    os.makedirs(output_path, exist_ok=True)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    df.to_csv(
        f"{output_path}/table_expression_per_phenotype.csv",
        index=False,
    )

    plt.show()


def plot_side_by_side_heatmaps(resistant_mean, sensitive_mean, folder_results):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sns.heatmap(
        resistant_mean.astype(float),
        annot=False,
        fmt=".2f",
        cmap="RdYlGn_r",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Mean Value"},
        ax=axes[0],
    )
    axes[0].set_title("Resistant Mean Phenotypes", fontsize=14)
    axes[0].set_ylabel("Condition", fontsize=12)
    axes[0].set_xlabel("Phenotype", fontsize=12)
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].tick_params(axis="y", labelsize=10)

    sns.heatmap(
        sensitive_mean.astype(float),
        annot=False,
        fmt=".2f",
        cmap="RdYlGn_r",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Mean Value"},
        ax=axes[1],
    )
    axes[1].set_title("Sensitive Mean Phenotypes", fontsize=14)
    axes[1].set_ylabel("Condition", fontsize=12)
    axes[1].set_xlabel("Phenotype", fontsize=12)
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].tick_params(axis="y", labelsize=10)

    plt.tight_layout()

    output_path = f"{folder_results}/output"
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(
        f"{output_path}/heatmap_resistant_vs_sensitive.png",
        dpi=300,
    )
    plt.show()
