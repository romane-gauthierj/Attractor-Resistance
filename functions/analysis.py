import pandas as pd
import os
import tempfile
import shutil
import ast
from scipy.stats import mannwhitneyu
import numpy as np


import logging

logger = logging.getLogger(__name__)

from functions.analysis_utils.MaBoSS_simulation.maboss_phenotype_patient import (
    compute_phenotype_table,
    collect_group_data,
    compute_mean_phenotype_values,
    compute_phenotype_table_custom_inputs,
)
from functions.analysis_utils.stats.stats_proba import (
    compute_mannwhitneyu_test_means,
    compute_power_calculation,
)

from functions.analysis_utils.results_MaBoSS_visualization.boxplot_phenotype import (
    create_boxplot,
    create_boxplot_three_groups,
)

from functions.analysis_utils.results_MaBoSS_visualization.create_phenotypes_patients_table import (
    plot_side_by_side_heatmaps,
    plot_three_side_by_side_heatmaps,
)

from functions.analysis_utils.genes_intervention.pers_interventions import (
    tailor_bnd_genes_intervention,
)

from functions.analysis_utils.results_MaBoSS_visualization.patients_ids_phenotype_table import (
    create_combined_table_patients,
)
from functions.analysis_utils.gene_enrichment.genes_signature import (
    create_results_gene_enrichment,
)
from functions.analysis_utils.identify_cell_lines import get_cell_lines

from functions.analysis_utils.stats.stats_proba import compute_kruskal_test_means, compute_stats_test_after_ko


def downstream_analysis(
    folder_name, 
    subdir,
    folder_result,
    folder_models,
    drug_interest,
    top_resistant_ids,
    top_sensitive_ids,
    top_healthy_ids,
    patients_categ,
    inputs_list,
    phenotype_interest,
    annotations_models,
    intervention_gene=None,
    list_active_inputs=None,
):
    # patients_categ = ["resistant", "sensitive", "healthy"]

    for patient_categ in patients_categ:
        folder_models_temp = f"{folder_models}/{patient_categ}"
        folder_models_temp_single_input = f"{folder_models}/{patient_categ}/pers_models"
        if list_active_inputs is None:
            folder_results_temp = f"{folder_result}/results/{patient_categ}"

        else:
            folder_results_temp = f"{folder_result}/results_more_inputs/{patient_categ}"

        os.makedirs(folder_results_temp, exist_ok=True)

        if patient_categ == "resistant":
            top_patients_ids = top_resistant_ids
        elif patient_categ == "sensitive":
            top_patients_ids = top_sensitive_ids
        elif patient_categ == "healthy":
            top_patients_ids = top_healthy_ids

        # get cell line distribution for each group of patients
        get_cell_lines(
            top_patients_ids, patient_categ, annotations_models, folder_result
        )

        # # compute phenotype table for each patient (attractor distribution)
        for patient in top_patients_ids:
            
            if list_active_inputs is None:
                # run for only one input ON at a time
                compute_phenotype_table(
                    folder_results_temp,
                    folder_models_temp_single_input,
                    patient,
                    inputs_list,
                    phenotype_interest,
                    drug_interest,
                )
            else:
                compute_phenotype_table_custom_inputs(
                    folder_results_temp,
                    folder_models_temp,
                    patient,
                    inputs_list,
                    phenotype_interest,
                    drug_interest,
                    list_active_inputs,
                )

        # # Group data together for each group of patient
        # top_patients_ids only to have the prefix
        collect_group_data(folder_results_temp, top_patients_ids)

    if list_active_inputs is None:
        folder_results_temp = f"{folder_result}/results"
    else:
        folder_results_temp = f"{folder_result}/results_more_inputs"

    folder_result_resistant = f"{folder_results_temp}/resistant"
    folder_result_sensitive = f"{folder_results_temp}/sensitive"
    folder_result_healthy = f"{folder_results_temp}/healthy"

    df_res_combined = pd.read_csv(
        f"{folder_result_resistant}/combined_results.csv", index_col=0
    )  # index: inputs, columns: phenotypes
    df_sens_combined = pd.read_csv(
        f"{folder_result_sensitive}/combined_results.csv", index_col=0
    )
    df_healthy_combined = pd.read_csv(
        f"{folder_result_healthy}/combined_results.csv", index_col=0
    )

    #  STATISTICS TESTS
    folder_result_temp_stats = f"{folder_results_temp}/statistics"
    if not os.path.exists(folder_result_temp_stats):
        os.makedirs(folder_result_temp_stats)

    # # two sides
    p_values_df_res_sens = compute_mannwhitneyu_test_means(
        folder_result_temp_stats,
        df_res_combined,
        df_sens_combined,
        drug_interest,
        gene=None,
    )

    p_values_df_res_healthy = compute_mannwhitneyu_test_means(
        folder_result_temp_stats,
        df_res_combined,
        df_healthy_combined,
        drug_interest,
        gene=None,
    )

    # # two sides p test
    if p_values_df_res_sens is not None:
        p_values_df_res_sens.to_csv(
            f"{folder_result_temp_stats}/p_values_df_mannwhitneyu_res_sens_{drug_interest}.csv",
            index=True,
        )
        data_res_sens_stats_sign = pd.read_csv(
            f"{folder_result_temp_stats}/p_values_df_mannwhitneyu_res_sens_{drug_interest}.csv",
            index_col=0,
        )
        create_boxplot(
            folder_results_temp,
            df_res_combined,
            df_sens_combined,
            data_res_sens_stats_sign,
        )

    else:
        logger.debug("Skipping boxplot for res_sensitive due to no p-values.")

    if p_values_df_res_healthy is not None:
        p_values_df_res_healthy.to_csv(
            f"{folder_result_temp_stats}/p_values_df_mannwhitneyu_res_healthy_{drug_interest}.csv",
            index=True,
        )
        data_res_healthy = pd.read_csv(
            f"{folder_result_temp_stats}/p_values_df_mannwhitneyu_res_healthy_{drug_interest}.csv",
            index_col=0,
        )
        create_boxplot_three_groups(
            folder_results_temp,
            df_res_combined,
            df_sens_combined,
            df_healthy_combined,
            data_res_healthy,
        )
    else:
        logger.debug("Skipping boxplot for res_healthy due to no p-values.")

    # Create a heatmap of phenotype distribution under diff growth condition between resistant and sensitive

    folder_result_resistant = f"{folder_results_temp}/resistant"
    folder_result_sensitive = f"{folder_results_temp}/sensitive"
    folder_result_healthy = f"{folder_results_temp}/healthy"

    df_res_combined = pd.read_csv(
        f"{folder_result_resistant}/combined_results.csv", index_col=0
    )  # index: inputs, columns: phenotypes
    df_sens_combined = pd.read_csv(
        f"{folder_result_sensitive}/combined_results.csv", index_col=0
    )
    df_healthy_combined = pd.read_csv(
        f"{folder_result_healthy}/combined_results.csv", index_col=0
    )

    patient_resistant_mean = compute_mean_phenotype_values(df_res_combined)
    patient_sensitive_mean = compute_mean_phenotype_values(df_sens_combined)
    patient_healthy_mean = compute_mean_phenotype_values(df_healthy_combined)




    # resistant and sensitive patients
    plot_side_by_side_heatmaps(
        patient_resistant_mean, patient_sensitive_mean, folder_results_temp
    )

    # plot_side_by_side_heatmaps(
    #     patient_resistant_median, patient_sensitive_median, folder_results_temp
    # )

    # resistant, sensitive and healthy patients
    plot_three_side_by_side_heatmaps(
        patient_resistant_mean,
        patient_sensitive_mean,
        patient_healthy_mean,
        folder_results_temp,
        labels=patients_categ,
    )



    rna_seq_data_filtered = pd.read_csv(
        f"analysis/{drug_interest}/{folder_name}/{subdir}/data_filtered/rna_seq_data_filtered.csv"
    )


    res_tables_path = f"{folder_result_resistant}"
    sens_tables_path = f"{folder_result_sensitive}"

    greater_sign_results = pd.read_csv(
        f"{folder_results_temp}/statistics/p_values_df_mannwhitneyu_res_sens_{drug_interest}.csv"
    )

    conditions_phenotypes_df = greater_sign_results[["Condition", "Phenotype"]]

    # combine all the patients together (resistant and sensitive)
    patients_phenot_table = create_combined_table_patients(
        res_tables_path, sens_tables_path, folder_results_temp
    )

    # do power calculation (assess the required size of each condition-phenotype pair)
    compute_power_calculation(df_res_combined, df_sens_combined, folder_results_temp)

    # gene differential expression analysis between resistant and sensitive
    create_results_gene_enrichment(
        rna_seq_data_filtered,
        patients_phenot_table,
        top_resistant_ids,
        top_sensitive_ids,
        conditions_phenotypes_df,
        folder_results_temp,
        annotations_models,
    )


    # if intervention genes (mannwhiteney test stats)
    # downstream analysis - CTNNB1 knockout 

    # input_inter = 'EGF'
    # phenotype_inter = 'Invasion'
    # if intervention_gene is not None:
    #     compute_stats_test_after_ko(drug_interest, discrete_variable,continuous_variable, drug_target, subdir, norm_technique ,  'EGF', 'Invasion')

   


    # compute kruskal test for healthy, resistant and sensitive patients
    # group_files = {
    #     patients_categ[0]: os.path.join(
    #         folder_result_resistant, "combined_results.csv"
    #     ),
    #     patients_categ[1]: os.path.join(
    #         folder_result_sensitive, "combined_results.csv"
    #     ),
    #     patients_categ[2]: os.path.join(folder_result_healthy, "combined_results.csv"),
    # }

    # kruskal_adjusted_df, significant_df = compute_kruskal_test_means(
    #     group_files, folder_result_temp_stats, patients_categ
    # )
    # kruskal_adjusted_df.to_csv(
    #     f"{folder_result_temp_stats}/kruskal_stats_healthy_res_sens_{drug_interest}.csv",
    #     index=True,
    # )
    # significant_df.to_csv(
    #     f"{folder_result_temp_stats}/kruskal_significant_{drug_interest}.csv",
    #     index=True,
    # )
