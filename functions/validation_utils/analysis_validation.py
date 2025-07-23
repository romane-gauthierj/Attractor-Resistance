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
    plot_two_stacked_heatmaps,
)
from functions.analysis_utils.stats.stats_proba import compute_kruskal_test_means


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
    group_categories,
):
    for patient in patients_id:
        group = patients_groups[patients_groups["sampleID"] == patient][
            "Gleason_group"
        ].iloc[0]
        folder_results_temp = os.path.join(folder_results, group)
        os.makedirs(folder_results_temp, exist_ok=True)

        folder_validation_temp = os.path.join(folder_validation, group)

        # compute the phenotype distribution for each condition for all the group
        # compute_phenotype_table(
        #     folder_results_temp,
        #     folder_validation_temp,
        #     patient,
        #     inputs_list,
        #     phenotype_interest,
        #     context_label,
        # )

    # compute the mean of each attractor for each condition for all th group
    # compute_phenotype_mean_group_validation(groups, folder_results)

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

    # plot_three_side_by_side_heatmaps(
    #     mean1, mean2, mean3, folder_results, labels=["high", "middle", "low"]
    # )

    # can also just look at low and high gleason score
    plot_two_stacked_heatmaps(mean1, mean3, folder_results, labels=["high", "low"])

    # combine all the data for each group
    for group in groups:
        folder_group_temp = os.path.join(folder_results, group)
        collect_group_data(folder_group_temp, patients_id)

    # compute stats test (Kruskal, p adjusted values)
    group_files = {
        group_categories[0]: os.path.join(
            folder_results, group_categories[0], "combined_results.csv"
        ),
        group_categories[1]: os.path.join(
            folder_results, group_categories[1], "combined_results.csv"
        ),
        group_categories[2]: os.path.join(
            folder_results, group_categories[2], "combined_results.csv"
        ),
    }

    adjusted_df = compute_kruskal_test_means(
        group_files, folder_results, group_categories
    )
    return adjusted_df
