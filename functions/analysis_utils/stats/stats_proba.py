import pandas as pd
import json
import numpy as np
from scipy.stats import shapiro, ttest_ind, mannwhitneyu, kruskal
import ast
import os
import seaborn as sns
import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower


def compute_kruskal_test_means(group_files, folder_results, group_categories):
    # Load all groups into dict of DataFrames
    group_dfs = {}
    for group, path in group_files.items():
        df = pd.read_csv(path, index_col=0)
        df = df.applymap(ast.literal_eval)
        group_dfs[group] = df

    # Get all inputs and phenotypes from one dataframe (assuming all share the same shape)
    inputs = group_dfs[group_categories[0]].index
    phenotypes = group_dfs[group_categories[0]].columns
    # Prepare result storage
    kruskal_results = pd.DataFrame(index=inputs, columns=phenotypes)

    # Run Kruskal-Wallis test for each (input, phenotype)
    for input_name in inputs:
        for phenotype in phenotypes:
            data_low = group_dfs[group_categories[0]].at[input_name, phenotype]
            data_medium = group_dfs[group_categories[1]].at[input_name, phenotype]
            data_high = group_dfs[group_categories[2]].at[input_name, phenotype]

            # Run the Kruskal-Wallis test only if all groups have data
            if data_low and data_medium and data_high:
                stat, pvalue = kruskal(data_low, data_medium, data_high)
                kruskal_results.at[input_name, phenotype] = pvalue
            else:
                kruskal_results.at[input_name, phenotype] = None

    pvals = kruskal_results.values.flatten()
    pvals = [p for p in pvals if p is not None]

    # Adjust using BH method
    _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")

    adjusted_df = kruskal_results.copy()

    # Fill adjusted p-values sequentially where there was a non-None p-value
    idx = 0
    for i in adjusted_df.index:
        for j in adjusted_df.columns:
            if adjusted_df.at[i, j] is not None:
                adjusted_df.at[i, j] = pvals_adj[idx]
                idx += 1
    # Save to CSV
    # adjusted_df.to_csv(
    #     os.path.join(folder_results_adj_df, "kruskal_pvalues_adjusted.csv")
    # )

    # Create a long-format DataFrame with only significant p-values
    significant_df = (
        adjusted_df.stack()
        .reset_index()
        .rename(
            columns={
                "level_0": "Condition",
                "level_1": "Phenotype",
                0: "Adjusted_P_Value",
            }
        )
    )

    significant_df = significant_df[significant_df["Adjusted_P_Value"] <= 0.05]
    significant_df = significant_df.dropna()

    return adjusted_df, significant_df


# compute normality (Shapiro-Wilk test p-value ≤ 0.05 → Reject H₀ → Data is not normally distributed.)
def compute_mannwhitneyu_test_means(
    folder,
    group_1_stats_values,
    group_2_stats_values,
    drug_interest,
    gene=None,
):
    group_1_stats_values.index.name = "Condition"
    group_2_stats_values.index.name = "Condition"

    phenotypes_list_res = group_1_stats_values.columns
    conditions_list_res = group_1_stats_values.index
    phenotypes_list_sens = group_2_stats_values.columns
    conditions_list_sens = group_2_stats_values.index

    p_values_records_mannwhitneyu_two_sides = []
    p_values_records_mannwhitneyu_greater = []

    for phenotype in phenotypes_list_res:
        if phenotype in phenotypes_list_sens:
            for condition in conditions_list_res:
                if condition in conditions_list_sens:
                    stats_data_res = group_1_stats_values.loc[condition, phenotype]
                    stats_data_sens = group_2_stats_values.loc[condition, phenotype]

                    stats_data_res = ast.literal_eval(stats_data_res)
                    stats_data_sens = ast.literal_eval(stats_data_sens)

                    statistic, p_value_res_sens = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="two-sided"
                    )

                    p_values_records_mannwhitneyu_two_sides.append(
                        {
                            "Condition": condition,
                            "Phenotype": phenotype,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens,
                        }
                    )

    # Create DataFrames from the records
    p_values_df_mannwhitneyu_two_sides = pd.DataFrame(
        p_values_records_mannwhitneyu_two_sides
    )
    p_values_df_mannwhitneyu_two_sides.set_index(
        ["Condition", "Phenotype"], inplace=True
    )

    # check which one has a pvalue < 0.05

    p_values_df_mannwhitneyu_two_sides = p_values_df_mannwhitneyu_two_sides[
        p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] <= 0.05
    ]

    choices = ["***", "**", "*"]
    conditions = [
        p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] < 0.001,
        (p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] >= 0.001)
        & (p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] < 0.01),
        (p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] >= 0.01)
        & (
            p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] <= 0.05
        ),
    ]

    p_values_df_mannwhitneyu_two_sides.loc[:, "Star_significant"] = np.select(
        conditions, choices, default=""
    )
    # p_values_df_mannwhitneyu_two_sides.to_csv(
    #     f"{folder}/p_values_df_mannwhitneyu_two_sides_{drug_interest}.csv",
    #     index=True,
    # )
    return p_values_df_mannwhitneyu_two_sides


def compute_power_calculation(df_res_combined, df_sens_combined, folder_results_temp):
    """
    Compute the number of patients required for a power of 0.8 and alpha of 0.05
    """

    conditions = df_res_combined.index
    phenotypes = df_res_combined.columns

    nb_patients_required = pd.DataFrame(index=conditions, columns=phenotypes)

    for condition, values in df_res_combined.iterrows():
        for phenotype in df_res_combined.columns:
            data_res = ast.literal_eval(df_res_combined.loc[condition, phenotype])
            data_sens = ast.literal_eval(df_sens_combined.loc[condition, phenotype])

            res_patients_mean = np.mean(data_res)
            sens_patients_mean = np.mean(data_sens)

            res_patients_std = np.std(data_res)
            sens_patients_std = np.std(data_sens)

            n_res_group = len(data_res)
            n_sens_group = len(data_sens)
            sd_pooled = np.sqrt(
                (
                    (n_res_group - 1) * res_patients_std**2
                    + (n_sens_group - 1) * sens_patients_std**2
                )
                / (n_res_group + n_sens_group - 2)
            )

            d = (res_patients_mean - sens_patients_mean) / sd_pooled

            analysis = TTestIndPower()
            sample_size = analysis.solve_power(
                effect_size=d, power=0.8, alpha=0.05, alternative="two-sided"
            )
            nb_patients_required.loc[condition, phenotype] = sample_size


    output_path = os.path.join(
        folder_results_temp, "output", "power_calculation_nb_patients_table.png"
    )
    # visually respresentation
    nb_patients_required = nb_patients_required.astype(float).round(1)

    # Create figure
    fig, ax = plt.subplots(
        figsize=(
            nb_patients_required.shape[1] * 1.5,
            nb_patients_required.shape[0] * 0.6,
        )
    )

    # Hide axes
    ax.axis("off")

    # Create the table
    table = ax.table(
        cellText=nb_patients_required.values,
        rowLabels=nb_patients_required.index,
        colLabels=nb_patients_required.columns,
        cellLoc="center",
        rowLoc="center",
        loc="center",
    )

    # Set general table styling
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)

    # Set header background color
    header_color = "#f0f0f0"  # Light gray
    row_header_color = "#f7f7f7"

    # Apply styling to headers only
    for (row, col), cell in table.get_celld().items():
        if row == 0 and col == -1:
            # Top-left corner (blank) cell
            cell.set_facecolor("white")
            cell.set_text_props(weight="bold")
        elif row == 0:
            # Column headers
            cell.set_facecolor(header_color)
            cell.set_text_props(weight="bold")
        elif col == -1:
            # Row headers
            cell.set_facecolor(row_header_color)
            cell.set_text_props(weight="bold")
        else:
            # All other cells (keep clean)
            cell.set_facecolor("white")

    # Set global font (e.g., Arial) if available
    plt.rcParams["font.family"] = "Arial"

    # Add title
    plt.title(
        "Required Patients per Input-Phenotype (Power=0.8, α=0.05)", fontsize=12, pad=20
    )

    # Save the image
   
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close()


def compute_power_calculation_genes(
    gene, rna_seq_data_filtered, top_resistant_ids, top_sensitive_ids
):
    # power calculation for FOXA1 -> do it also for E2F1

    group_gene_res = list(
        rna_seq_data_filtered[
            (rna_seq_data_filtered["model_id"].isin(top_resistant_ids))
            & (rna_seq_data_filtered["gene_symbol"] == gene)
        ]["rsem_tpm"]
    )

    group_gene_sens = list(
        rna_seq_data_filtered[
            (rna_seq_data_filtered["model_id"].isin(top_sensitive_ids))
            & (rna_seq_data_filtered["gene_symbol"] == gene)
        ]["rsem_tpm"]
    )

    res_patients_mean = np.mean(group_gene_res)
    sens_patients_mean = np.mean(group_gene_sens)

    res_patients_std = np.std(group_gene_res)
    sens_patients_std = np.std(group_gene_sens)

    n_res_group = len(top_resistant_ids)
    n_sens_group = len(top_sensitive_ids)

    sd_pooled = np.sqrt(
        (
            (n_res_group - 1) * res_patients_std**2
            + (n_sens_group - 1) * sens_patients_std**2
        )
        / (n_res_group + n_sens_group - 2)
    )

    d = res_patients_mean - sens_patients_mean / sd_pooled

    # Set parameters
    effect_size = d  # Cohen's d
    alpha = 0.05  # Significance level
    power = 0.80  # Desired power
    alternative = "two-sided"  # Equivalent to "two.sample" in R

    # Initialize power analysis object
    analysis = TTestIndPower()

    # Compute required sample size per group
    sample_size = analysis.solve_power(
        effect_size=effect_size, power=power, alpha=alpha, alternative=alternative
    )
    return sample_size






# Mannwhiteney stats test before and after KO 
def compute_stats_test_after_ko(res_values, sens_values, res_values_knockout_ctnnb1, sens_values_knockout_ctnnb1, output_folder , input_interest, phenotype_interest):



    res_values_knockout = res_values_knockout_ctnnb1.loc[input_interest, phenotype_interest]
    res_values = res_values.loc[input_interest, phenotype_interest]
    res_values_knockout = ast.literal_eval(res_values_knockout)
    res_values = ast.literal_eval(res_values)



    sens_values_knockout = sens_values_knockout_ctnnb1.loc[input_interest, phenotype_interest]
    sens_values = sens_values.loc[input_interest, phenotype_interest]
    sens_values_knockout = ast.literal_eval(sens_values_knockout)
    sens_values = ast.literal_eval(sens_values)

    stat_sens, p_value_sens = mannwhitneyu(sens_values_knockout, sens_values, alternative='two-sided')
    stat_res, p_value_res = mannwhitneyu(res_values_knockout, res_values, alternative='two-sided')


    results_pval_ko = pd.DataFrame(index = ['resistant', 'sensitive'], columns = ['stat', 'mannwhitneyu_pval', 'mean_change_KO_vs_baseline'])
    results_pval_ko.loc['resistant','stat'] = stat_res
    results_pval_ko.loc['sensitive','stat'] = stat_sens


    results_pval_ko.loc['resistant','mannwhitneyu_pval'] = p_value_res
    results_pval_ko.loc['sensitive','mannwhitneyu_pval'] = p_value_sens


    results_pval_ko.loc['resistant','mean_change_KO_vs_baseline'] = np.mean(sens_values_knockout) - np.mean(sens_values)
    results_pval_ko.loc['sensitive','mean_change_KO_vs_baseline'] = np.mean(res_values_knockout) - np.mean(res_values)

    results_pval_ko.to_csv(
        f'{output_folder}/pval_test_after_ko_{phenotype_interest}.csv'
    )
















