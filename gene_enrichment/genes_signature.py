# Identify what genes are differently expressed in the resistant with high proliferation upon EGF

import pandas as pd
import numpy as np
import scipy.stats as stats


# compute genes mean and variability
def compute_genes_mean_signature(
    rna_seq_data,
    folder,
    phenotype,
    condition,
    data_phenotype_patients,
    top_resistant_ids,
    top_sensitive_ids,
):
    # identify the differently genes (higher expressed) in the resistant group compared to the sensitive group upon a specific phenotype- condition
    data_phenotype_patients["Model_ID"] = (
        data_phenotype_patients["Unnamed: 0"].astype(str).str.split("_").str[0]
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

    group_phenotype_resistant = data_phenotype_patients[
        (data_phenotype_patients["Drug status"] == "Resistant")
        & (data_phenotype_patients[f"{condition}_{phenotype}"] >= 0.1)
    ]

    # 2 groups: resistant_proliferating_group and sensitive_group_ids
    resistant_group_ids = group_phenotype_resistant["Model_ID"].tolist()
    sensitive_group = data_phenotype_patients[
        data_phenotype_patients["Drug status"] == "Sensitive"
    ]
    sensitive_group_ids = sensitive_group["Model_ID"].tolist()

    # len_resistant_group_ids = len(resistant_group_ids)
    # len_sensitive_group_ids = len(sensitive_group_ids)

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

    significant_genes = genes_stats_results[genes_stats_results["P-value"] <= 0.05]
    significant_genes.to_csv(
        f"{folder}/sensitive_resistant_results/genes_diff_expressed/significant_genes_{condition}_ON_{phenotype}.csv",
        index=True,
    )
    return significant_genes, resistant_group_ids, sensitive_group_ids
