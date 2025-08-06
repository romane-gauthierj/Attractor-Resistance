# Create Boxplots of phenotype distribution


import pandas as pd
import ast
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.patches as mpatches
import os
import logging

logger = logging.getLogger(__name__)

# Step 1: Convert the dataframes to long format and min-max normalization
# def melt_df(df, group_label):
#     records = []
#     phenotypes = df.index

#     conditions = df.columns
#     for phenotype in phenotypes:
#         all_values = []
#         for condition in conditions:
#             values = ast.literal_eval(df.loc[phenotype, condition])
#             all_values.extend(values)
#         # Min-max scaling
#         min_val, max_val = min(all_values), max(all_values)
#         normalized_values = [
#             (v - min_val) / (max_val - min_val + 1e-9) for v in all_values
#         ]  # Avoid division by zero

#         i = 0
#         for condition in conditions:
#             values = ast.literal_eval(df.loc[phenotype, condition])
#             for _ in values:
#                 records.append(
#                     {
#                         "Phenotype": phenotype,
#                         "Condition": condition,
#                         "Group": group_label,
#                         "Expression": normalized_values[i],
#                         "Condition_Group": f"{condition} ({group_label})",
#                     }
#                 )
#                 i += 1
#     return pd.DataFrame(records)


def create_boxplot(folder_result, res_data, sens_data, significant_df):
    # if only significatif results printed:
    # if significant_df is None or significant_df.empty:
    #     logger.debug("No significant results to plot. Skipping boxplot creation.")
    #     return

    # phenotypes = significant_df["Phenotype"].unique()
    # conditions = significant_df.index.values
    # # only unique values
    # conditions = list(set(conditions))

    phenotypes = list(res_data.columns)
    conditions = list(res_data.index)

    n_phenotypes = len(phenotypes)
    fig, axes = plt.subplots(
        1, n_phenotypes, figsize=(10 * n_phenotypes, 8), sharey=True
    )

    for i, phenotype in enumerate(phenotypes):
        data = []
        group_positions = []
        group_labels = []

        # Store individual box positions (for plotting)
        resistant_positions = []
        sensitive_positions = []
        resistant_data = []
        sensitive_data = []

        tick = 1  # x position for groups

        for condition in conditions:
            res_value = res_data.loc[condition, phenotype]
            res_value_list = ast.literal_eval(res_value)
            # ensure each value is a list
            if not isinstance(res_value_list, list):
                logger.debug(
                    f"Warning: Resistant value for {condition}, {phenotype} is not a list: {res_value_list}"
                )

            resistant_data.append(res_value_list)
            resistant_positions.append(tick - 0.2)

            sens_value = sens_data.loc[condition, phenotype]
            sens_value_list = ast.literal_eval(sens_value)
            if not isinstance(sens_value_list, list):
                logger.debug(
                    f"Warning: Sensitive value for {condition}, {phenotype} is not a list: {sens_value_list}"
                )

            sensitive_data.append(sens_value_list)
            sensitive_positions.append(tick + 0.2)

            group_positions.append(tick)
            group_labels.append(condition)

            tick += 1

        ax = axes[i] if n_phenotypes > 1 else axes

        # Plot resistant and sensitive as separate sets with offset
        box1 = ax.boxplot(
            resistant_data, positions=resistant_positions, patch_artist=True
        )
        box2 = ax.boxplot(
            sensitive_data, positions=sensitive_positions, patch_artist=True
        )

        # Color the boxes
        for patch in box1["boxes"]:
            patch.set_facecolor("#FF7F0E")  # Resistant = orange
        for patch in box2["boxes"]:
            patch.set_facecolor("#008000")  # Sensitive = green
        
          # Set median lines to black
        for median in box1["medians"]:
            median.set_color("black")
            median.set_linewidth(2)  # Optional: make median line thicker
        for median in box2["medians"]:
            median.set_color("black")
            median.set_linewidth(2)  #

        # Set x-ticks at the center of each group
        ax.set_xticks(group_positions)
        ax.set_xticklabels(group_labels, rotation=45, ha="right")
        ax.set_xlabel("Condition")
        ax.set_title(phenotype)
        if i == 0:
            ax.set_ylabel("Expression Values")

        # Add significance stars (above the higher of the two boxplots)
        for j, condition in enumerate(conditions):
            row = significant_df[
                (significant_df.index.values == condition)
                & (significant_df["Phenotype"] == phenotype)
            ]
            if not row.empty:
                star = row["Star_significant"].values[0]
                if star:
                    max_val = max(max(resistant_data[j]), max(sensitive_data[j]))
                    ax.text(
                        group_positions[j],
                        max_val + 0.01,
                        star,
                        ha="center",
                        va="bottom",
                        fontsize=16,
                        color="red",
                    )

    # Add legend
    resistant_patch = mpatches.Patch(color="#FF7F0E", label="Resistant")
    sensitive_patch = mpatches.Patch(color="#008000", label="Sensitive")
    fig.legend(
        handles=[resistant_patch, sensitive_patch],
        loc="upper right",
        fontsize=12,
        title="Group",
    )

    output_path = f"{folder_result}/output"
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(
        f"{output_path}/boxplot_expression_per_phenotype.png",
        dpi=400,
        bbox_inches="tight",
    )

    plt.tight_layout()
    plt.show()


def create_boxplot_three_groups(
    folder_result, res_data, sens_data, healthy_data, data_kruskal_significant
):
    if data_kruskal_significant is None or data_kruskal_significant.empty:
        return
    # Ensure all input dataframes are indeed pandas DataFrames
    # This is a crucial check to prevent the AttributeError: 'str' object has no attribute 'loc'
    if (
        not isinstance(res_data, pd.DataFrame)
        or not isinstance(sens_data, pd.DataFrame)
        or not isinstance(healthy_data, pd.DataFrame)
    ):
        return

    if data_kruskal_significant.empty:
        logger.debug(" No significant phenotypes found for plotting.")
        return

    phenotypes = data_kruskal_significant["Phenotype"].unique().tolist()

    n_phenotypes = len(phenotypes)

    fig_width = max(10, 5 * n_phenotypes)
    fig, axes = plt.subplots(1, n_phenotypes, figsize=(fig_width, 8), sharey=True)

    if n_phenotypes == 1:
        axes = [axes]  # Ensure axes is always an iterable

    for i, phenotype in enumerate(phenotypes):
        current_phenotype_significant_data = data_kruskal_significant[
            data_kruskal_significant["Phenotype"] == phenotype
        ]

        if current_phenotype_significant_data.empty:
            logger.debug(f" No significant data for phenotype: {phenotype}. Skipping.")
            continue

        significant_conditions = current_phenotype_significant_data.index.tolist()

        group_positions = []
        group_labels = []

        resistant_data_for_plot = []
        sensitive_data_for_plot = []
        healthy_data_list_for_plot = []

        resistant_positions = []
        sensitive_positions = []
        healthy_positions = []

        tick = 1  # x-axis position tracker

        for condition in significant_conditions:
            try:
                # --- Resistant ---
                # Check if condition and phenotype exist in res_data's index/columns
                if condition in res_data.index and phenotype in res_data.columns:
                    res_value = res_data.loc[condition, phenotype]
                    res_value_list = ast.literal_eval(
                        str(res_value)
                    )  # Ensure it's a string before literal_eval
                    resistant_data_for_plot.append(res_value_list)
                    resistant_positions.append(tick - 0.3)
                else:
    
                    continue

                # --- Sensitive ---
                if condition in sens_data.index and phenotype in sens_data.columns:
                    sens_value = sens_data.loc[condition, phenotype]
                    sens_value_list = ast.literal_eval(str(sens_value))
                    sensitive_data_for_plot.append(sens_value_list)
                    sensitive_positions.append(tick)
                else:
                    logger.debug(
                        f" Missing sensitive data for condition '{condition}', phenotype '{phenotype}'. Skipping for this condition."
                    )
                    continue

                # --- Healthy ---
                if (
                    condition in healthy_data.index
                    and phenotype in healthy_data.columns
                ):
                    healthy_value = healthy_data.loc[condition, phenotype]
                    healthy_value_list = ast.literal_eval(str(healthy_value))
                    healthy_data_list_for_plot.append(healthy_value_list)
                    healthy_positions.append(tick + 0.3)
                else:
                    logger.debug(
                        f" Missing healthy data for condition '{condition}', phenotype '{phenotype}'. Skipping for this condition."
                    )
                    continue

                group_positions.append(tick)
                group_labels.append(condition)
                tick += 1

            except (ValueError, SyntaxError) as e:
                logger.debug(
                    f" Error processing data for condition '{condition}', phenotype '{phenotype}': {e}"
                )
                logger.debug(f"Problematic res_value: '{res_data.loc[condition, phenotype]}'")
                # Optionally, you might choose to skip this specific condition or handle it differently
                continue  # Skip to the next condition if there's a parsing error

        if (
            not resistant_data_for_plot
            and not sensitive_data_for_plot
            and not healthy_data_list_for_plot
        ):
            logger.debug(
                f" No valid data to plot for phenotype: {phenotype}. Skipping plot for this phenotype."
            )
            continue

        ax = axes[i]

        # Plot all three categories
        # Check if lists are not empty before plotting to avoid error with empty data
        if resistant_data_for_plot:
            box1 = ax.boxplot(
                resistant_data_for_plot,
                positions=resistant_positions,
                patch_artist=True,
                widths=0.2,
            )
            for patch in box1["boxes"]:
                patch.set_facecolor("#FF7F0E")  # Orange = Resistant

            # Set median lines to black
            for median in box1["medians"]:
                median.set_color("black")
                median.set_linewidth(2)

        else:
            box1 = {"boxes": []}  # Empty dict to prevent error in legend loop

        if sensitive_data_for_plot:
            box2 = ax.boxplot(
                sensitive_data_for_plot,
                positions=sensitive_positions,
                patch_artist=True,
                widths=0.2,
            )
            for patch in box2["boxes"]:
                patch.set_facecolor("#008000")  # Green = Sensitive

             # Set median lines to black
            for median in box2["medians"]:
                median.set_color("black")
                median.set_linewidth(2)
        else:
            box2 = {"boxes": []}

        if healthy_data_list_for_plot:
            box3 = ax.boxplot(
                healthy_data_list_for_plot,
                positions=healthy_positions,
                patch_artist=True,
                widths=0.2,
            )
            for patch in box3["boxes"]:
                patch.set_facecolor("#1f77b4")  # Blue = Healthy
            for median in box3["medians"]:
                median.set_color("black")
                median.set_linewidth(2)
        else:
            box3 = {"boxes": []}

        # Set x-ticks and labels
        ax.set_xticks(group_positions)
        ax.set_xticklabels(group_labels, rotation=45, ha="right")
        ax.set_xlabel("Condition")
        ax.set_title(phenotype)
        if i == 0:
            ax.set_ylabel("Expression Values")

        # ---------- ADD SIGNIFICANCE STARS ----------
        # Re-get the row as current_phenotype_significant_data might have been filtered
        # It's better to use the already filtered dataframe here
        for j, condition in enumerate(significant_conditions):
            # Find the specific row for the current condition and phenotype
            row = current_phenotype_significant_data[
                current_phenotype_significant_data.index == condition
            ]
            if not row.empty:
                p_val = row["Mannwhitneyu_P_value_Resistant"].values[0]
                star = ""
                if p_val < 0.001:
                    star = "***"
                elif p_val < 0.01:
                    star = "**"
                elif p_val < 0.05:
                    star = "*"

                if star:
                    # Calculate max_y dynamically based on the actual plotted data for this group
                    # Ensure indices are valid for the lists holding data for plot
                    current_group_max_y_values = []
                    if j < len(resistant_data_for_plot) and resistant_data_for_plot[j]:
                        current_group_max_y_values.extend(resistant_data_for_plot[j])
                    if j < len(sensitive_data_for_plot) and sensitive_data_for_plot[j]:
                        current_group_max_y_values.extend(sensitive_data_for_plot[j])
                    if (
                        j < len(healthy_data_list_for_plot)
                        and healthy_data_list_for_plot[j]
                    ):
                        current_group_max_y_values.extend(healthy_data_list_for_plot[j])

                    if current_group_max_y_values:
                        max_y = max(current_group_max_y_values)
                        ax.text(
                            group_positions[j],
                            max_y * 1.05,  # Place star slightly above the highest point
                            star,
                            ha="center",
                            va="bottom",
                            fontsize=16,
                            color="red",
                        )

    # Add legend to the figure, not an individual axis
    resistant_patch = mpatches.Patch(color="#FF7F0E", label="Resistant")
    sensitive_patch = mpatches.Patch(color="#008000", label="Sensitive")
    healthy_patch = mpatches.Patch(color="#1f77b4", label="Healthy")
    fig.legend(
        handles=[resistant_patch, sensitive_patch, healthy_patch],
        loc="upper right",
        bbox_to_anchor=(1.05, 1),  # Position legend outside the plot area
        fontsize=12,
        title="Group",
    )

    # Adjust layout before saving to prevent labels from overlapping
    plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust rect to make space for legend

    output_path = f"{folder_result}/output"
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(
        f"{output_path}/boxplot_expression_per_phenotype_three_groups.png",
        dpi=300,
        bbox_inches="tight",  # Ensure everything is saved
    )

    plt.show()
