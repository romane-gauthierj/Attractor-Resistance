import os
import pandas as pd
import ast
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

from functions.analysis_utils.MaBoSS_simulation.maboss_phenotype_patient import (
    compute_phenotype_table,
    compute_phenotype_mean_group_validation,
)

from functions.analysis_utils.MaBoSS_simulation.maboss_phenotype_patient import (
    collect_group_data,
)

from functions.analysis_utils.results_MaBoSS_visualization.create_phenotypes_patients_table import (
    plot_three_side_by_side_heatmaps,
)


def validation_analysis(
    folder_results,
    folder_validation,
    patients_id,
    patients_groups,
    inputs_list,
    phenotype_interest,
    groups,
    context_label,
    type_models,
):
    for patient in patients_id:
        group = patients_groups[patients_groups["sampleID"] == patient][
            "Gleason_group"
        ].iloc[0]
        folder_results_temp = os.path.join(folder_results, group)
        os.makedirs(folder_results_temp, exist_ok=True)

        folder_validation_temp = os.path.join(folder_validation, group)

        # compute the phenotype distribution for each condition for all the group
        compute_phenotype_table(
            folder_results_temp,
            folder_validation_temp,
            patient,
            inputs_list,
            phenotype_interest,
            context_label,
        )

    # compute the mean of each attractor for each condition for all th group
    compute_phenotype_mean_group_validation(groups, folder_results)

    # heatmap of the mean of each group
    mean1 = pd.read_csv(
        f"analysis/validation/{type_models}/results/high_aggressive/high_aggressive_mean_phenotype_values.csv",
        index_col=0,
    )
    mean2 = pd.read_csv(
        f"analysis/validation/{type_models}/results/middle_aggressive/middle_aggressive_mean_phenotype_values.csv",
        index_col=0,
    )
    mean3 = pd.read_csv(
        f"analysis/validation/{type_models}/results/low_aggressive/low_aggressive_mean_phenotype_values.csv",
        index_col=0,
    )

    plot_three_side_by_side_heatmaps(
        mean1, mean2, mean3, folder_results, labels=["high", "middle", "low"]
    )

    # combine all the data for each group
    for group in groups:
        folder_group_temp = os.path.join(folder_results, group)
        collect_group_data(folder_group_temp, patients_id)

    # compute stats test (Kruskal, p adjusted values)
    group_files = {
        "low": os.path.join(folder_results, "low_aggressive", "combined_results.csv"),
        "medium": os.path.join(
            folder_results, "middle_aggressive", "combined_results.csv"
        ),
        "high": os.path.join(folder_results, "high_aggressive", "combined_results.csv"),
    }

    # Load all groups into dict of DataFrames
    group_dfs = {}
    for group, path in group_files.items():
        df = pd.read_csv(path, index_col=0)
        df = df.applymap(ast.literal_eval)
        group_dfs[group] = df

    # Get all inputs and phenotypes from one dataframe (assuming all share the same shape)
    inputs = group_dfs["low"].index
    phenotypes = group_dfs["low"].columns
    # Prepare result storage
    kruskal_results = pd.DataFrame(index=inputs, columns=phenotypes)

    # Run Kruskal-Wallis test for each (input, phenotype)
    for input_name in inputs:
        for phenotype in phenotypes:
            data_low = group_dfs["low"].at[input_name, phenotype]
            data_medium = group_dfs["medium"].at[input_name, phenotype]
            data_high = group_dfs["high"].at[input_name, phenotype]

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
    folder_results_adj_df = os.path.join(folder_results, "output")
    os.makedirs(folder_results_adj_df, exist_ok=True)

    adjusted_df.to_csv(
        os.path.join(folder_results_adj_df, "kruskal_pvalues_adjusted.csv")
    )
