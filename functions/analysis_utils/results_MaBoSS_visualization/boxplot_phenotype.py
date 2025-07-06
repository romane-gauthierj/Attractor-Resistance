# Create Boxplots of phenotype distribution


import pandas as pd
import ast
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.patches as mpatches
import os

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
    # res_data.rename(columns={"Unnamed: 0": "Condition"}, inplace=True)
    # sens_data.rename(columns={"Unnamed: 0": "Condition"}, inplace=True)
    if significant_df is None or significant_df.empty:
        print("No significant results to plot. Skipping boxplot creation.")
        return

    phenotypes = significant_df["Phenotype"].unique()
    conditions = pd.unique(significant_df.index.values)

    # phenotypes = list(res_data.columns[1:])
    # conditions = list(res_data.index)

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
                print(
                    f"Warning: Resistant value for {condition}, {phenotype} is not a list: {res_value_list}"
                )

            resistant_data.append(res_value_list)
            resistant_positions.append(tick - 0.2)

            sens_value = sens_data.loc[condition, phenotype]
            sens_value_list = ast.literal_eval(sens_value)
            if not isinstance(sens_value_list, list):
                print(
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
        dpi=300,
        bbox_inches="tight",
    )

    plt.tight_layout()
    plt.show()


def create_boxplot_three_groups(
    folder_result, res_data, sens_data, healthy_data, data_kruskal_significant
):
    # Ensure all input dataframes have the same structure

    if data_kruskal_significant.empty:
        print("⚠️ No significant phenotypes found for plotting.")
        return
    phenotypes = list(set(data_kruskal_significant["Phenotype"].unique()))
    conditions = list(set(data_kruskal_significant["Condition"].unique()))

    n_phenotypes = len(phenotypes)
    fig, axes = plt.subplots(
        1, n_phenotypes, figsize=(10 * n_phenotypes, 8), sharey=True
    )

    if n_phenotypes == 1:
        axes = [axes]

    for i, phenotype in enumerate(phenotypes):
        if data_kruskal_significant is not None and not data_kruskal_significant.empty:
            significant_conditions = data_kruskal_significant[
                data_kruskal_significant["Phenotype"] == phenotype
            ]["Condition"].unique()
        else:
            continue  # Skip phenotype if no significant data

        if len(significant_conditions) == 0:
            continue

        group_positions = []
        group_labels = []

        resistant_data = []
        sensitive_data = []
        healthy_data_list = []

        resistant_positions = []
        sensitive_positions = []
        healthy_positions = []

        tick = 1  # x-axis position tracker

        for condition in significant_conditions:
            # --- Resistant ---
            res_value = res_data.loc[condition, phenotype]
            res_value_list = ast.literal_eval(res_value)
            resistant_data.append(res_value_list)
            resistant_positions.append(tick - 0.3)

            # --- Sensitive ---
            sens_value = sens_data.loc[condition, phenotype]
            sens_value_list = ast.literal_eval(sens_value)
            sensitive_data.append(sens_value_list)
            sensitive_positions.append(tick)

            # --- Healthy ---
            healthy_value = healthy_data.loc[condition, phenotype]
            healthy_value_list = ast.literal_eval(healthy_value)
            healthy_data_list.append(healthy_value_list)
            healthy_positions.append(tick + 0.3)

            group_positions.append(tick)
            group_labels.append(condition)
            tick += 1

        ax = axes[i]

        # Plot all three categories
        box1 = ax.boxplot(
            resistant_data, positions=resistant_positions, patch_artist=True
        )
        box2 = ax.boxplot(
            sensitive_data, positions=sensitive_positions, patch_artist=True
        )
        box3 = ax.boxplot(
            healthy_data_list, positions=healthy_positions, patch_artist=True
        )

        # Color the boxes
        for patch in box1["boxes"]:
            patch.set_facecolor("#FF7F0E")  # Orange = Resistant
        for patch in box2["boxes"]:
            patch.set_facecolor("#008000")  # Green = Sensitive
        for patch in box3["boxes"]:
            patch.set_facecolor("#1f77b4")  # Blue = Healthy

        # Set x-ticks and labels
        ax.set_xticks(group_positions)
        ax.set_xticklabels(group_labels, rotation=45, ha="right")
        ax.set_xlabel("Condition")
        ax.set_title(phenotype)
        if i == 0:
            ax.set_ylabel("Expression Values")

        # ---------- ADD SIGNIFICANCE STARS ----------
        if data_kruskal_significant is not None and not data_kruskal_significant.empty:
            for j, condition in enumerate(significant_conditions):
                row = data_kruskal_significant[
                    (data_kruskal_significant["Condition"] == condition)
                    & (data_kruskal_significant["Phenotype"] == phenotype)
                ]
                if not row.empty:
                    p_val = row["Adjusted_P_Value"].values[0]
                    if p_val < 0.001:
                        star = "***"
                    elif p_val < 0.01:
                        star = "**"
                    elif p_val < 0.05:
                        star = "*"
                    else:
                        star = ""

                    if star:
                        max_y = max(
                            max(resistant_data[j]),
                            max(sensitive_data[j]),
                            max(healthy_data_list[j]),
                        )
                        ax.text(
                            group_positions[j],
                            max_y + 0.05 * max_y,
                            star,
                            ha="center",
                            va="bottom",
                            fontsize=16,
                            color="red",
                        )

    # Add legend
    resistant_patch = mpatches.Patch(color="#FF7F0E", label="Resistant")
    sensitive_patch = mpatches.Patch(color="#008000", label="Sensitive")
    healthy_patch = mpatches.Patch(color="#1f77b4", label="Healthy")
    fig.legend(
        handles=[resistant_patch, sensitive_patch, healthy_patch],
        loc="upper right",
        fontsize=12,
        title="Group",
    )

    output_path = f"{folder_result}/output"
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(
        f"{output_path}/boxplot_expression_per_phenotype_three_groups.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.tight_layout()
    plt.show()


# def create_boxplot_three_groups(folder_result, res_data, sens_data, healthy_data, data_kruskal_significant):
#     # Ensure all input dataframes have the same structure
#     phenotypes = res_data.columns
#     conditions = res_data.index

#     n_phenotypes = len(phenotypes)
#     fig, axes = plt.subplots(
#         1, n_phenotypes, figsize=(10 * n_phenotypes, 8), sharey=True
#     )

#     for i, phenotype in enumerate(phenotypes):
#         group_positions = []
#         group_labels = []

#         resistant_data = []
#         sensitive_data = []
#         healthy_data_list = []

#         resistant_positions = []
#         sensitive_positions = []
#         healthy_positions = []

#         tick = 1  # x-axis position tracker

#         for condition in conditions:
#             # --- Resistant ---
#             res_value = res_data.loc[condition, phenotype]
#             res_value_list = ast.literal_eval(res_value)
#             resistant_data.append(res_value_list)
#             resistant_positions.append(tick - 0.3)

#             # --- Sensitive ---
#             sens_value = sens_data.loc[condition, phenotype]
#             sens_value_list = ast.literal_eval(sens_value)
#             sensitive_data.append(sens_value_list)
#             sensitive_positions.append(tick)

#             # --- Healthy ---
#             healthy_value = healthy_data.loc[condition, phenotype]
#             healthy_value_list = ast.literal_eval(healthy_value)
#             healthy_data_list.append(healthy_value_list)
#             healthy_positions.append(tick + 0.3)

#             group_positions.append(tick)
#             group_labels.append(condition)
#             tick += 1

#         ax = axes[i] if n_phenotypes > 1 else axes

#         # Plot all three categories
#         box1 = ax.boxplot(resistant_data, positions=resistant_positions, patch_artist=True)
#         box2 = ax.boxplot(sensitive_data, positions=sensitive_positions, patch_artist=True)
#         box3 = ax.boxplot(healthy_data_list, positions=healthy_positions, patch_artist=True)

#         # Color the boxes
#         for patch in box1["boxes"]:
#             patch.set_facecolor("#FF7F0E")  # Orange = Resistant
#         for patch in box2["boxes"]:
#             patch.set_facecolor("#008000")  # Green = Sensitive
#         for patch in box3["boxes"]:
#             patch.set_facecolor("#1f77b4")  # Blue = Healthy

#         # Set x-ticks and labels
#         ax.set_xticks(group_positions)
#         ax.set_xticklabels(group_labels, rotation=45, ha="right")
#         ax.set_xlabel("Condition")
#         ax.set_title(phenotype)
#         if i == 0:
#             ax.set_ylabel("Expression Values")

#     # Add legend
#     resistant_patch = mpatches.Patch(color="#FF7F0E", label="Resistant")
#     sensitive_patch = mpatches.Patch(color="#008000", label="Sensitive")
#     healthy_patch = mpatches.Patch(color="#1f77b4", label="Healthy")
#     fig.legend(
#         handles=[resistant_patch, sensitive_patch, healthy_patch],
#         loc="upper right",
#         fontsize=12,
#         title="Group",
#     )

#     output_path = f"{folder_result}/output"
#     os.makedirs(output_path, exist_ok=True)
#     plt.savefig(
#         f"{output_path}/boxplot_expression_per_phenotype.png",
#         dpi=300,
#         bbox_inches="tight",
#     )

#     plt.tight_layout()
#     plt.show()


# patient_res_values = pd.read_csv(
#     "../results/Refametinib_PAN_CANCER/resistant_results/only_gene_expression/single_input_on/phenotype_distribution_patients/combined_results.csv"
# )
# patient_sens_values = pd.read_csv(
#     "../results/Refametinib_PAN_CANCER/sensitive_results/only_gene_expression/single_input_on/phenotype_distribution_patients/combined_results.csv"
# )
# data_greater_side = pd.read_csv(
#     "../results/Refametinib_PAN_CANCER/sensitive_resistant_results/p_values_df_mannwhitneyu_greater_sign_Refametinib.csv"
# )
# data_less_side = pd.read_csv(
#     "../results/Refametinib_PAN_CANCER/sensitive_resistant_results/p_values_df_mannwhitneyu_less_sign_Refametinib.csv"
# )
# data_both_side = pd.read_csv(
#     "../results/Refametinib_PAN_CANCER/sensitive_resistant_results/p_values_df_mannwhitneyu_two_sides_Refametinib.csv"
# )

# signif_df = pd.concat([data_greater_side, data_less_side, data_both_side], axis=0)


# test_create_boxplot(patient_res_values, patient_sens_values, signif_df)

# def create_boxplot(folder, patient_res_values, patient_sens_values, data_two_sides):
#     patient_res_values.set_index("Unnamed: 0", inplace=True)
#     patient_sens_values.set_index("Unnamed: 0", inplace=True)
#     df_res_long = melt_df(patient_res_values, "Resistant")
#     df_sens_long = melt_df(patient_sens_values, "Sensitive")

#     df_combined = pd.concat([df_res_long, df_sens_long], ignore_index=True)
#     df_combined = df_combined.loc[:, ~df_combined.columns.duplicated()]
#     df_combined["Condition_Group"] = (
#         df_combined["Condition"] + " (" + df_combined["Group"] + ")"
#     )

#     significant_pairs = {
#         (row["Phenotype"], row["Condition"]): row["Star_significant"]
#         for _, row in data_two_sides.iterrows()
#         if row["Star_significant"] != ""
#     }

#     # if only want to keep the phenotype-condition significant

#     #     df_combined = df_combined[
#     #     df_combined.apply(lambda row: (row['Phenotype'], row['Condition']) in significant_pairs, axis=1)
#     # ]

#     # Step 2: Plot — one figure per phenotype
#     g = sns.catplot(
#         data=df_combined,
#         kind="box",
#         x="Condition",
#         y="Expression",
#         hue="Group",  # Separates Resistant vs Sensitive
#         col="Phenotype",  # One subplot (facet) per phenotype
#         col_wrap=4,  # Adjust this value depending on the number of phenotypes
#         sharey=False,  # Each subplot gets its own y-axis range if needed
#         palette={"Resistant": "#FF7F0E", "Sensitive": "#008000"},
#         height=6,  # Height of each subplot
#         aspect=0.8,  # Aspect ratio of each subplot (width = height * aspect)
#     )

#     # Step 3: Apply significance annotation to each subplot
#     for ax in g.axes.flat:
#         # Get the current phenotype from the subplot title
#         phenotype = (
#             ax.get_title().split(" = ")[-1]
#             if " = " in ax.get_title()
#             else ax.get_title()
#         )
#         for (sig_phenotype, condition), stars in significant_pairs.items():
#             if sig_phenotype == phenotype:
#                 # Get the data for just this condition and phenotype
#                 subset = df_combined[
#                     (df_combined["Phenotype"] == phenotype)
#                     & (df_combined["Condition"] == condition)
#                 ]
#                 if not subset.empty:
#                     max_expression = subset["Expression"].max()

#                     # Find the tick positions and labels
#                     tick_labels = [t.get_text() for t in ax.get_xticklabels()]
#                     try:
#                         xpos = tick_labels.index(condition)
#                         ax.text(
#                             xpos,
#                             max_expression + 0.01,  # much closer to the box
#                             stars,
#                             ha="center",
#                             va="bottom",
#                             fontsize=16,
#                             fontweight="bold",
#                             color="red",
#                         )
#                     except ValueError:
#                         # Condition not found among x-tick labels for this subplot
#                         continue

#     # Improve formatting
#     g.set_xticklabels(rotation=45, ha="right", fontsize=10)
#     g.set_axis_labels("", "Expression")

#     g.set_titles("{col_name}", size=14, fontweight="bold")
#     g.fig.subplots_adjust(top=0.9)

#     # Move the legend to the top-right corner

#     g.legend.set_bbox_to_anchor((0.6, 0.9))  # Position the legend outside of the plot

#     # Set the font size for the legend
#     for legend_text in g.legend.get_texts():
#         legend_text.set_fontsize(10)

#     # Add the main title for the entire figure
#     # g.fig.suptitle("Boxplot of Expression Levels by Condition and Group\nfor Each Phenotype", fontsize=16)

#     # Tighten the layout to avoid overlap
#     plt.tight_layout()

#     output_path = f"{folder}/sensitive_resistant_results/boxplot_expression_per_phenotype.png"  # or use .pdf/.svg
#     g.legend.set_bbox_to_anchor((0.6, 0.9))
#     g.legend.get_frame().set_edgecolor("black")
#     g.legend.get_frame().set_linewidth(1)
#     g.legend.get_frame().set_boxstyle("round")
#     g.legend.set_title("Group", prop={"weight": "bold", "size": 12})

#     for legend_text in g.legend.get_texts():
#         legend_text.set_fontsize(10)
#     plt.savefig(output_path, dpi=300, bbox_inches="tight")

#     # Show the plot (optional if you're running this headlessly)
#     plt.show()


# # drug_interest = "Refametinib"  # Pictilisib, 'Avagacestat' AZD8931
# # tissue_interest = "PAN_CANCER"
# # tissue_remove = "Haematopoietic and Lymphoid"
# # folder_result = f"../results/{drug_interest}_{tissue_interest}"

# # patient_res_values = pd.read_csv(
# #     "../results/Refametinib_PAN_CANCER/resistant_results/only_gene_expression/single_input_on/phenotype_distribution_patients/combined_trans_results.csv",
# #     header=1,
# # )
# # patient_sens_values = pd.read_csv(
# #     "../results/Refametinib_PAN_CANCER/sensitive_results/only_gene_expression/single_input_on/phenotype_distribution_patients/combined_trans_results.csv",
# #     header=1,
# # )
# # data_greater_side = pd.read_csv(
# #     "../results/Refametinib_PAN_CANCER/sensitive_resistant_results/p_values_df_mannwhitneyu_greater_sign_Refametinib.csv"
# # )
# # print(patient_sens_values)


# patient_res_values_transposed = patient_res_values.T
# patient_sens_values_transposed = patient_sens_values.T

# create_boxplot(
#     folder_result,
#     patient_res_values,
#     patient_sens_values,
#     data_greater_side,
# )
