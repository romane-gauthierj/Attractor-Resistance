import pandas as pd
import os
import tempfile
import shutil

from functions.analysis_utils.MaBoSS_simulation.maboss_phenotype_patient import (
    compute_phenotype_table,
    collect_group_data,
    compute_mean_phenotype_values,
)
from functions.analysis_utils.stats.stats_proba import compute_mannwhitneyu_test_means

from functions.analysis_utils.results_MaBoSS_visualization.boxplot_phenotype import (
    create_boxplot,
)

from functions.analysis_utils.results_MaBoSS_visualization.create_phenotypes_patients_table import (
    plot_side_by_side_heatmaps,
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


def downstream_analysis(
    folder_result,
    drug_interest,
    top_resistant_ids,
    top_sensitive_ids,
    inputs_list,
    phenotype_interest,
):
    patients_categ = ["resistant", "sensitive"]

    for patient_categ in patients_categ:
        folder_results_temp = f"{folder_result}/{patient_categ}"
        folder_models_temp = f"analysis/{drug_interest}/models/{patient_categ}"

        top_patients_ids = (
            top_resistant_ids if patient_categ == "resistant" else top_sensitive_ids
        )

    # compute phenotype table for each patient (attractor distribution)
    for patient in top_patients_ids:
        compute_phenotype_table(
            folder_results_temp,
            folder_models_temp,
            patient,
            inputs_list,
            phenotype_interest,
            drug_interest,
        )

        # Group data together for each group of patient
        folder_models_res_temp = f"analysis/{drug_interest}/results/{patient_categ}"
        collect_group_data(folder_models_res_temp)

    folder_result_resistant = f"analysis/{drug_interest}/results/resistant"
    folder_result_sensitive = f"analysis/{drug_interest}/results/sensitive"

    df_res_combined = pd.read_csv(
        f"{folder_result_resistant}/combined_results.csv", index_col=0
    )  # index: inputs, columns: phenotypes
    df_sens_combined = pd.read_csv(
        f"{folder_result_sensitive}/combined_results.csv", index_col=0
    )

    # Compute statistics for each condition-phenotype pair
    folder_result_temp = f"analysis/{drug_interest}/results/statistics"
    if not os.path.exists(folder_result_temp):
        os.makedirs(folder_result_temp)

    compute_mannwhitneyu_test_means(
        folder_result_temp, df_res_combined, df_sens_combined, drug_interest
    )

    data_greater_sign = pd.read_csv(
        f"{folder_result_temp}/p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv",
        index_col=0,
    )

    # Create boxplot for each phenotype
    create_boxplot(folder_result, df_res_combined, df_sens_combined, data_greater_sign)

    folder_result_resistant = f"analysis/{drug_interest}/results/resistant"
    folder_result_sensitive = f"analysis/{drug_interest}/results/sensitive"

    df_res_combined = pd.read_csv(
        f"{folder_result_resistant}/combined_results.csv", index_col=0
    )  # index: inputs, columns: phenotypes
    df_sens_combined = pd.read_csv(
        f"{folder_result_sensitive}/combined_results.csv", index_col=0
    )

    patient_resistant_mean = compute_mean_phenotype_values(df_res_combined)
    patient_sensitive_mean = compute_mean_phenotype_values(df_sens_combined)

    # Create a heatmap of phenotype distribution under diff growth condition between resistant and sensitive
    plot_side_by_side_heatmaps(
        patient_resistant_mean, patient_sensitive_mean, folder_result
    )

    rna_seq_data_filtered = pd.read_csv(
        f"analysis/{drug_interest}/data_filtered/rna_seq_data_filtered.csv"
    )

    res_tables_path = f"{folder_result_resistant}"
    sens_tables_path = f"{folder_result_sensitive}"

    greater_sign_results = pd.read_csv(
        f"{folder_result}/statistics/p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv"
    )
    conditions_phenotypes_df = greater_sign_results[["Condition", "Phenotype"]]

    # combine all the patients together (resistant and sensitive)
    patients_phenot_table = create_combined_table_patients(
        res_tables_path, sens_tables_path, folder_result
    )
    patients_phenot_table.to_csv(f"analysis/{drug_interest}/results/outputs")

    # gene differential expression analysis
    create_results_gene_enrichment(
        rna_seq_data_filtered,
        patients_phenot_table,
        top_resistant_ids,
        top_sensitive_ids,
        conditions_phenotypes_df,
        folder_result,
    )


def downstream_analysis_interv(
    drug_interest,
    top_resistant_ids,
    top_sensitive_ids,
    inputs_list,
    phenotype_interest,
    list_genes_interv=None,
):
    patients_categ = ["resistant", "sensitive"]
    genes_to_target = []

    for patient_categ in patients_categ:
        top_patients_ids = (
            top_resistant_ids if patient_categ == "resistant" else top_sensitive_ids
        )

        # Create a temporary directory for this gene's models intervention
        with (
            tempfile.TemporaryDirectory() as temp_dir,
        ):
            # 2. Copy the original models into the temp directory
            models_folder_temp = f"analysis/{drug_interest}/models/{patient_categ}"
            shutil.copytree(models_folder_temp, temp_dir, dirs_exist_ok=True)

        if list_genes_interv is not None:
            for gene in list_genes_interv:
                tailor_bnd_genes_intervention(
                    gene, top_patients_ids, temp_dir, drug_interest
                )

        # 3. Compute phenotype table for each patient (attractor distribution) and Boxplot to vizualise the distribution and extract the significant phenotypes
        data_greater_sign = downstream_analysis(
            drug_interest,
            top_resistant_ids,
            top_sensitive_ids,
            inputs_list,
            phenotype_interest,
        )

        if "Proliferation" not in data_greater_sign["Phenotype"].values:
            # Ensured gene is always defined when appending.
            if gene is not None:
                genes_to_target.append(gene)
        else:
            print("Gene does not block proliferation", gene)
    if not genes_to_target:
        print("No genes to target found in the analysis.")
    else:
        return genes_to_target
