import pandas as pd
import json
import numpy as np
from scipy.stats import shapiro, ttest_ind, mannwhitneyu, kruskal
import ast
import os
from statsmodels.stats.power import TTestIndPower


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


def compute_kruskal(inputs_list, phenotype_interest, df_res_combined, df_sens_combined):
    kruskal_results = pd.DataFrame(index=inputs_list, columns=phenotype_interest)

    for input_name in inputs_list:
        for phenotype in phenotype_interest:
            group1 = df_res_combined.at[input_name, phenotype]
            group2 = df_sens_combined.at[input_name, phenotype]
            group1 = ast.literal_eval(group1)
            group2 = ast.literal_eval(group2)

            if group1 and group2:
                stat, pvalue = kruskal(group1, group2)
                kruskal_results.at[input_name, phenotype] = pvalue
            else:
                kruskal_results.at[input_name, phenotype] = None

    significant_kruskal_results = kruskal_results.apply(
        lambda x: x if (x is not None and float(x) < 0.05) else np.nan
    )
    return significant_kruskal_results


# compute normality (Shapiro-Wilk test p-value ≤ 0.05 → Reject H₀ → Data is not normally distributed.)
def compute_mannwhitneyu_test_means(
    folder,
    patient_res_stats_values,
    patient_sens_stats_values,
    drug_interest,
    gene=None,
):
    patient_res_stats_values.index.name = "Condition"
    patient_sens_stats_values.index.name = "Condition"

    phenotypes_list_res = patient_res_stats_values.columns
    conditions_list_res = patient_res_stats_values.index
    phenotypes_list_sens = patient_sens_stats_values.columns
    conditions_list_sens = patient_sens_stats_values.index

    p_values_records_mannwhitneyu_two_sides = []
    p_values_records_mannwhitneyu_less = []
    p_values_records_mannwhitneyu_greater = []

    for phenotype in phenotypes_list_res:
        if phenotype in phenotypes_list_sens:
            for condition in conditions_list_res:
                if condition in conditions_list_sens:
                    stats_data_res = patient_res_stats_values.loc[condition, phenotype]
                    stats_data_sens = patient_sens_stats_values.loc[
                        condition, phenotype
                    ]

                    stats_data_res = ast.literal_eval(stats_data_res)
                    stats_data_sens = ast.literal_eval(stats_data_sens)

                    statistic, p_value_res_sens = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="two-sided"
                    )
                    statistic, p_value_res_sens_less = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="less"
                    )
                    statistic, p_value_res_sens_greater = mannwhitneyu(
                        stats_data_res, stats_data_sens, alternative="greater"
                    )

                    p_values_records_mannwhitneyu_two_sides.append(
                        {
                            "Condition": condition,
                            "Phenotype": phenotype,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens,
                        }
                    )

                    p_values_records_mannwhitneyu_less.append(
                        {
                            "Condition": condition,
                            "Phenotype": phenotype,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens_less,
                        }
                    )

                    p_values_records_mannwhitneyu_greater.append(
                        {
                            "Condition": condition,
                            "Phenotype": phenotype,
                            "Mannwhitneyu_P_value_Resistant": p_value_res_sens_greater,
                        }
                    )
    # Create DataFrames from the records
    p_values_df_mannwhitneyu_two_sides = pd.DataFrame(
        p_values_records_mannwhitneyu_two_sides
    )
    p_values_df_mannwhitneyu_two_sides.set_index(
        ["Condition", "Phenotype"], inplace=True
    )

    p_values_df_mannwhitneyu_less = pd.DataFrame(p_values_records_mannwhitneyu_less)
    p_values_df_mannwhitneyu_less.set_index(["Condition", "Phenotype"], inplace=True)

    p_values_df_mannwhitneyu_greater = pd.DataFrame(
        p_values_records_mannwhitneyu_greater
    )
    p_values_df_mannwhitneyu_greater.set_index(["Condition", "Phenotype"], inplace=True)

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
        f"{folder}/p_values_df_mannwhitneyu_two_sides_{drug_interest}.csv",
        index=True,
    )
    p_values_df_mannwhitneyu_less_sign.to_csv(
        f"{folder}/p_values_df_mannwhitneyu_less_sign_{drug_interest}.csv",
        index=True,
    )

    filename = f"{folder}/p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv"
    if gene is not None:
        filename = (
            f"{folder}/{gene}_p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv"
        )

    p_values_df_mannwhitneyu_greater_sign.to_csv(filename, index=True)


def compute_power_calculation(df_res_combined, df_sens_combined):
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

    return nb_patients_required


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
    # print(f"Required sample size per group: {float(sample_size):.2f}")
