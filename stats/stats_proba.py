import pandas as pd
import json
import numpy as np
from scipy.stats import shapiro, ttest_ind, mannwhitneyu, kruskal
import ast
import os


def compute_kruskal_test_means_validation(data_groups):
    results = []

    for (input_name, phenotype), group_vals in data_groups.items():
        groupA = group_vals["GroupA"]
        groupB = group_vals["GroupB"]

        if len(groupA) > 0 and len(groupB) > 0:
            stat, p = kruskal(groupA, groupB)
            results.append(
                {
                    "Input": input_name,
                    "Phenotype": phenotype,
                    "Kruskal_Stat": stat,
                    "p_value": p,
                    "Significant": p <= 0.05,
                }
            )

        df_results = pd.DataFrame(results)
        significant_results = df_results[df_results["Significant"] == True]
    return significant_results


# compute normality (Shapiro-Wilk test p-value ≤ 0.05 → Reject H₀ → Data is not normally distributed.)
def compute_mannwhitneyu_test_means(
    folder, patient_res_stats_values, patient_sens_stats_values, drug_interest
):
    # compute the statistics test and assign to the related star
    patient_res_stats_values.set_index(
        patient_res_stats_values.columns[0], inplace=True
    )
    patient_sens_stats_values.set_index(
        patient_sens_stats_values.columns[0], inplace=True
    )

    patient_res_stats_values.index.name = "Phenotype"
    patient_sens_stats_values.index.name = "Phenotype"

    phenotypes_list_res = patient_res_stats_values.index
    conditions_list_res = patient_res_stats_values.columns
    phenotypes_list_sens = patient_sens_stats_values.index
    conditions_list_sens = patient_sens_stats_values.columns

    # p_values_records_shapiro = []
    p_values_records_mannwhitneyu_two_sides = []
    p_values_records_mannwhitneyu_less = []
    p_values_records_mannwhitneyu_greater = []

    for phenotype in phenotypes_list_res:
        if phenotype in phenotypes_list_sens:
            for condition in conditions_list_res:
                if condition in conditions_list_sens:
                    stats_data_res = patient_res_stats_values.loc[phenotype, condition]

                    stats_data_sens = patient_sens_stats_values.loc[
                        phenotype, condition
                    ]

                    stats_data_res = ast.literal_eval(stats_data_res)
                    stats_data_sens = ast.literal_eval(stats_data_sens)

                    # statistic, p_value_res = shapiro(stats_data_res)
                    # statistic, p_value_sens = shapiro(stats_data_sens)

                    statistic, p_value_res_sens = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="two-sided"
                    )
                    statistic, p_value_res_sens_less = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="less"
                    )
                    statistic, p_value_res_sens_greater = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="greater"
                    )

                    # Shapiro: Append the row as a dictionary
                    # p_values_records_shapiro.append({
                    #     'Phenotype': phenotype,
                    #     'Condition': condition,
                    #     'Shapiro_P_value_Resistant': p_value_res,
                    #     'Shapiro_P_value_Sensitive': p_value_sens
                    # })

                    # Mannwhitneyu: Append the row as a dictionary
                    p_values_records_mannwhitneyu_two_sides.append(
                        {
                            "Phenotype": phenotype,
                            "Condition": condition,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens,
                        }
                    )

                    p_values_records_mannwhitneyu_less.append(
                        {
                            "Phenotype": phenotype,
                            "Condition": condition,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens_less,
                        }
                    )

                    p_values_records_mannwhitneyu_greater.append(
                        {
                            "Phenotype": phenotype,
                            "Condition": condition,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens_greater,
                        }
                    )

                    # print(f"Shapiro-Wilk Test for Resistant: p-value = {p_value_res}")
                    # print(f"Shapiro-Wilk Test for Sensitive: p-value = {p_value_sens}")

    # Create the DataFrames
    # p_values_df_shapiro = pd.DataFrame(p_values_records_shapiro)
    # p_values_df_shapiro.set_index(['Phenotype', 'Condition'], inplace=True)

    # print(p_values_df_shapiro)

    p_values_df_mannwhitneyu_two_sides = pd.DataFrame(
        p_values_records_mannwhitneyu_two_sides
    )
    p_values_df_mannwhitneyu_two_sides.set_index(
        ["Phenotype", "Condition"], inplace=True
    )

    p_values_df_mannwhitneyu_less = pd.DataFrame(p_values_records_mannwhitneyu_less)
    p_values_df_mannwhitneyu_less.set_index(["Phenotype", "Condition"], inplace=True)

    p_values_df_mannwhitneyu_greater = pd.DataFrame(
        p_values_records_mannwhitneyu_greater
    )
    p_values_df_mannwhitneyu_greater.set_index(["Phenotype", "Condition"], inplace=True)

    # check which one has a pvalue < 0.05

    p_values_df_mannwhitneyu_two_sides = p_values_df_mannwhitneyu_two_sides[
        p_values_df_mannwhitneyu_two_sides["Mannwhitneyu_P_value_Resistant"] <= 0.05
    ]
    p_values_df_mannwhitneyu_less_sign = p_values_df_mannwhitneyu_less[
        p_values_df_mannwhitneyu_less["Mannwhitneyu_P_value_Resistant"] <= 0.05
    ].copy()
    p_values_df_mannwhitneyu_greater_sign = p_values_df_mannwhitneyu_greater[
        p_values_df_mannwhitneyu_greater["Mannwhitneyu_P_value_Resistant"] <= 0.05
    ].copy()

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

    conditions = [
        p_values_df_mannwhitneyu_less_sign["Mannwhitneyu_P_value_Resistant"] < 0.001,
        (p_values_df_mannwhitneyu_less_sign["Mannwhitneyu_P_value_Resistant"] >= 0.001)
        & (p_values_df_mannwhitneyu_less_sign["Mannwhitneyu_P_value_Resistant"] < 0.01),
        (p_values_df_mannwhitneyu_less_sign["Mannwhitneyu_P_value_Resistant"] >= 0.01)
        & (
            p_values_df_mannwhitneyu_less_sign["Mannwhitneyu_P_value_Resistant"] <= 0.05
        ),
    ]

    p_values_df_mannwhitneyu_less_sign.loc[:, "Star_significant"] = np.select(
        conditions, choices, default=""
    )

    conditions = [
        p_values_df_mannwhitneyu_greater_sign["Mannwhitneyu_P_value_Resistant"] < 0.001,
        (
            p_values_df_mannwhitneyu_greater_sign["Mannwhitneyu_P_value_Resistant"]
            >= 0.001
        )
        & (
            p_values_df_mannwhitneyu_greater_sign["Mannwhitneyu_P_value_Resistant"]
            < 0.01
        ),
        (
            p_values_df_mannwhitneyu_greater_sign["Mannwhitneyu_P_value_Resistant"]
            >= 0.01
        )
        & (
            p_values_df_mannwhitneyu_greater_sign["Mannwhitneyu_P_value_Resistant"]
            <= 0.05
        ),
    ]

    p_values_df_mannwhitneyu_greater_sign.loc[:, "Star_significant"] = np.select(
        conditions, choices, default=""
    )

    p_values_df_mannwhitneyu_two_sides.to_csv(
        f"{folder}/sensitive_resistant_results/p_values_df_mannwhitneyu_two_sides_{drug_interest}.csv",
        index=True,
    )
    p_values_df_mannwhitneyu_less_sign.to_csv(
        f"{folder}/sensitive_resistant_results/p_values_df_mannwhitneyu_less_sign_{drug_interest}.csv",
        index=True,
    )
    p_values_df_mannwhitneyu_greater_sign.to_csv(
        f"{folder}/sensitive_resistant_results/p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv",
        index=True,
    )
