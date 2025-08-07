# Identify what genes are differently expressed in the resistant with high proliferation upon EGF

import pandas as pd
import numpy as np
import scipy.stats as stats
import os
from statsmodels.stats.multitest import multipletests


# compute genes mean and variability
def compute_genes_mean_signature(
    rna_seq_data,
    folder,
    phenotype,
    condition,
    data_phenotype_patients,
    top_resistant_ids,
    top_sensitive_ids,
    annotations_models,
):
    # identify the differently genes (higher expressed) in the resistant group compared to the sensitive group upon a specific phenotype- condition
    # data_phenotype_patients["Model_ID"] = (
    #     data_phenotype_patients.index
    # )
    data_phenotype_patients = data_phenotype_patients.reset_index().rename(
        columns={"index": "Model_ID"}
    )

    conditions = [
        data_phenotype_patients["Model_ID"].isin(top_resistant_ids),
        data_phenotype_patients["Model_ID"].isin(top_sensitive_ids),

    ]

    choices = ["Resistant", "Sensitive"]
    data_phenotype_patients.loc[:, "Drug status"] = np.select(
        conditions, choices, default=""
    )
    # resistant group changes according to what is the condition and the phenotype
    # group_proliferation_resistant: group with high phenotype

    # adjust the expression values here
    group_phenotype_resistant = data_phenotype_patients[
        (data_phenotype_patients["Drug status"] == "Resistant")
        & (data_phenotype_patients[f"{condition}_{phenotype}"] >= 0.8)
    ]

    # 2 groups: resistant_proliferating_group and sensitive_group_ids
    resistant_group_ids = group_phenotype_resistant["Model_ID"].tolist()

    group_phenotype_sensitive = data_phenotype_patients[
        (data_phenotype_patients["Drug status"] == "Sensitive")
        & (data_phenotype_patients[f"{condition}_{phenotype}"] <= 0.4)
    ]

    sensitive_group_ids = group_phenotype_sensitive["Model_ID"].tolist()

    # annotations
    tissue_status_distribution_res = annotations_models[
        annotations_models["model_id"].isin(resistant_group_ids)
    ]["tissue_status"].value_counts()

    tissue_status_distribution_sens = annotations_models[
        annotations_models["model_id"].isin(sensitive_group_ids)
    ]["tissue_status"].value_counts()

    gender_distribution_res = annotations_models[
        annotations_models["model_id"].isin(resistant_group_ids)
    ]["gender"].value_counts()

    gender_distribution_sens = annotations_models[
        annotations_models["model_id"].isin(sensitive_group_ids)
    ]["gender"].value_counts()

    tissue_distribution_res = annotations_models[
        annotations_models["model_id"].isin(resistant_group_ids)
    ]["tissue"].value_counts()

    tissue_distribution_sens = annotations_models[
        annotations_models["model_id"].isin(sensitive_group_ids)
    ]["tissue"].value_counts()

    ratio_gender_res = gender_distribution_res / len(resistant_group_ids)
    ratio_gender_sens = gender_distribution_sens / len(sensitive_group_ids)

    ratio_tissue_status_res = tissue_status_distribution_res / len(resistant_group_ids)
    ratio_tissue_status_sens = tissue_status_distribution_sens / len(
        sensitive_group_ids
    )

    ratio_tissue_res = tissue_distribution_res / len(resistant_group_ids)
    ratio_tissue_sens = tissue_distribution_sens / len(sensitive_group_ids)

    # extract gene expression data
    rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]

    mean_gene_rsem = (
        rna_seq_data.groupby(["model_id", "gene_symbol"])["rsem_tpm"]
        .mean()
        .reset_index()
    )
    # compute mean and sd of each gene in the two groups
    genes_names = mean_gene_rsem["gene_symbol"].unique().tolist()
    genes_stats_results = pd.DataFrame(
        index=genes_names, columns=["Group Resistant Mean", "Group Sensitive Mean"]
    )

    for gene in genes_names:
        # Filter by gene and by patient group
        group_resistant_phenotype = mean_gene_rsem[
            (mean_gene_rsem["gene_symbol"] == gene)
            & (mean_gene_rsem["model_id"].isin(resistant_group_ids))
        ]["rsem_tpm"]

        group_sensitive = mean_gene_rsem[
            (mean_gene_rsem["gene_symbol"] == gene)
            & (mean_gene_rsem["model_id"].isin(sensitive_group_ids))
        ]["rsem_tpm"]

        # Compute the means
        group_resistant_phenotype_mean = group_resistant_phenotype.mean()
        group_sensitive_mean = group_sensitive.mean()

        # Compute the variance
        group_resistant_variance = group_resistant_phenotype.var(
            ddof=1
        )  # sample variance
        group_sensitive_variance = group_sensitive.var(ddof=1)

        # Perform statistical test (T-test or Mann-Whitney U test)
        # Check if data is normally distributed (Shapiro-Wilk test)
        _, p_normal_group_1 = stats.shapiro(group_resistant_phenotype)
        _, p_normal_group_2 = stats.shapiro(group_sensitive)

        # If both groups are normally distributed, use t-test
        if p_normal_group_1 > 0.05 and p_normal_group_2 > 0.05:
            t_stat, p_value = stats.ttest_ind(
                group_resistant_phenotype, group_sensitive
            )
        else:
            # If either group is not normally distributed, use Mann-Whitney U test
            u_stat, p_value = stats.mannwhitneyu(
                group_resistant_phenotype, group_sensitive
            )

        # Assign results to DataFrame
        genes_stats_results.at[gene, "Group Resistant Mean"] = (
            group_resistant_phenotype_mean
        )
        genes_stats_results.at[gene, "Group Sensitive Mean"] = group_sensitive_mean
        genes_stats_results.at[gene, "Group Resistant Variance"] = (
            group_resistant_variance
        )
        genes_stats_results.at[gene, "Group Sensitive Variance"] = (
            group_sensitive_variance
        )

        genes_stats_results.at[gene, "P-value"] = p_value





        ## TEST TEST TEST - p adjusted value 
        genes_stats_results["P-value"] = genes_stats_results["P-value"].astype(float)
        pvals = genes_stats_results["P-value"].values
        # Compute adjusted p-values (FDR, Benjamini-Hochberg)
        _, pvals_adj, _, _ = multipletests(pvals, method="fdr_bh")
        genes_stats_results["P-adj"] = pvals_adj



    significant_genes = genes_stats_results[genes_stats_results["P-adj"] <= 0.05]
    output_dir = f"{folder}/genes_diff_expressed"
    os.makedirs(output_dir, exist_ok=True)
    file_path = f"{output_dir}/significant_genes_{condition}_ON_{phenotype}.csv"

    # Save the main DataFrame of significant genes
    significant_genes.to_csv(file_path, index=True)

    # Append additional information
    with open(file_path, "a") as f:
        f.write("\n\n--- Summary Annotations for Resistant Group ---\n")

        f.write("\nTissue Status Distribution:\n")
        ratio_tissue_status_res.to_csv(f, header=False)

        f.write(f"\nResistant Group Size,{len(resistant_group_ids)}\n")
        f.write("\nGender Ratio (Resistant Group / Total):\n")
        ratio_gender_res.to_csv(f, header=False)

        f.write("\nTissue Ratio (Resistant Group / Total):\n")
        ratio_tissue_res.to_csv(f, header=False)

        # Sensitive Group Annotations
        f.write("\n\n--- Summary Annotations for Sensitive Group ---\n")

        f.write(f"Sensitive Group Size,{len(sensitive_group_ids)}\n")

        f.write("\nTissue Status Distribution:\n")
        ratio_tissue_status_sens.to_csv(f, header=False)

        f.write(f"\nSensitive Group Size,{len(sensitive_group_ids)}\n")

        f.write("\nGender Ratio (Sensitive Group / Total):\n")
        ratio_gender_sens.to_csv(f, header=False)

        f.write("\nTissue Ratio (Sensitive Group / Total):\n")
        ratio_tissue_sens.to_csv(f, header=False)


def create_results_gene_enrichment(
    rna_seq_data_filtered,
    patients_phenot_table,
    top_resistant_ids,
    top_sensitive_ids,
    conditions_phenotypes_df,
    folder_result,
    annotations_models,
):
    # save every gene enrichment dataframe in a dictionary of index condtion_phenotype

    # indexes_list = []
    # for condition, phenotype in conditions_phenotypes_df.values:
    #     indexes_list.append(f"{condition}_{phenotype}")
    # results_genes_enrich = {key: None for key in indexes_list}

    for condition, phenotype in conditions_phenotypes_df.values:
        compute_genes_mean_signature(
            rna_seq_data_filtered,
            folder_result,
            phenotype,
            condition,
            patients_phenot_table,
            top_resistant_ids,
            top_sensitive_ids,
            annotations_models,
        )
        # if isinstance(genes_stats_results, pd.DataFrame):
        #     results_genes_enrich[f"{condition}_{phenotype}"] = genes_stats_results
