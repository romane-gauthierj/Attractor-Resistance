
# Identify what genes are differently expressed in the resistant with high proliferation upon EGF

import pandas as pd 
from get_patients_sens_res import get_patients_top_10
import numpy as np
import scipy.stats as stats





def compute_genes_mean_signature(folder, montagud_nodes, phenotype,condition, data_phenotype_patients,top_resistant_ids, top_sensitive_ids):
# identify the differently genes (higher expressed) in the resistant group compared to the sensitive group upon a specific phenotype- condition
    data_phenotype_patients['Model_ID'] = data_phenotype_patients['Unnamed: 0'].astype(str).str.split('_').str[0]
    print(data_phenotype_patients)

    conditions = [
        data_phenotype_patients['Model_ID'].isin(top_resistant_ids),
        data_phenotype_patients['Model_ID'].isin(top_sensitive_ids)
        ]

    choices = ['Resistant', 'Sensitive']
    data_phenotype_patients.loc[:,'Drug status'] = np.select(conditions, choices, default = '')
    print(data_phenotype_patients)


# resistant group changes according to what is the condition and the phenotype
# group_proliferation_resistant: group with high phenotype 

    group_phenotype_resistant = data_phenotype_patients[
        (data_phenotype_patients['Drug status'] == 'Resistant') & 
        (data_phenotype_patients[f'{condition}_ON_{phenotype}'] >= 0.1)
    ]
    
    # 2 groups: resistant_proliferating_group and sensitive_group_ids
    resistant_group_ids = group_phenotype_resistant['Model_ID'].tolist()
    sensitive_group = data_phenotype_patients[data_phenotype_patients['Drug status'] == 'Sensitive']
    sensitive_group_ids = sensitive_group['Model_ID'].tolist()

    # extract gene expression data 
    patients_ids = top_resistant_ids + top_sensitive_ids
    genes_data = (
    pd.read_csv('../cancer_data/rnaseq_merged/rnaseq_merged_20250117.csv')
    .query('`model_id` in @patients_ids and `gene_symbol` in @montagud_nodes')
    .loc[:, ['model_id', 'gene_symbol', 'rsem_tpm']]
)
    
    mean_gene_rsem = (
        genes_data
        .groupby(['model_id', 'gene_symbol'])['rsem_tpm']
        .mean()
        .reset_index()
    )
    # compute mean and sd of each gene in the two groups
    genes_names = mean_gene_rsem['gene_symbol'].unique().tolist()
    genes_stats_results = pd.DataFrame(index=genes_names, columns=['Group Resistant', 'Group Sensitive'])

    for gene in genes_names:
        # Filter by gene and by patient group
        group_resistant_phenotype = mean_gene_rsem[
            (mean_gene_rsem['gene_symbol'] == gene) &
            (mean_gene_rsem['model_id'].isin(resistant_group_ids))
        ]['rsem_tpm']
        
        group_sensitive = mean_gene_rsem[
            (mean_gene_rsem['gene_symbol'] == gene) &
            (mean_gene_rsem['model_id'].isin(sensitive_group_ids))
        ]['rsem_tpm']


        # Compute the means
        group_resistant_phenotype_mean = group_resistant_phenotype.mean()
        group_sensitive_mean = group_sensitive.mean()


          # Perform statistical test (T-test or Mann-Whitney U test)
        # Check if data is normally distributed (Shapiro-Wilk test)
        _, p_normal_group_1 = stats.shapiro(group_resistant_phenotype)
        _, p_normal_group_2 = stats.shapiro(group_sensitive)

        # If both groups are normally distributed, use t-test
        if p_normal_group_1 > 0.05 and p_normal_group_2 > 0.05:
            t_stat, p_value = stats.ttest_ind(group_resistant_phenotype, group_sensitive)
        else:
            # If either group is not normally distributed, use Mann-Whitney U test
            u_stat, p_value = stats.mannwhitneyu(group_resistant_phenotype, group_sensitive)



        # Assign results to DataFrame
        genes_stats_results.at[gene, 'Group Resistant'] = group_resistant_phenotype_mean
        genes_stats_results.at[gene, 'Group Sensitive'] = group_sensitive_mean
        genes_stats_results.at[gene, 'P-value'] = p_value

    significant_genes = genes_stats_results[genes_stats_results['P-value'] <= 0.05]
    significant_genes.to_csv(f'{folder}/sensitive_resistant_results/genes_diff_expressed/significant_genes_{condition}_ON_{phenotype}.csv', index=True)
    print(significant_genes)
    return significant_genes













# def compute_genes_std_signature(data_phenotype_patients, top_resistant_ids, top_sensitive_ids):

#     data_phenotype_resistant = data_phenotype_patients[data_phenotype_patients['Unnamed: 0'].isin(top_resistant_ids)]
#     data_phenotype_resistant_filtered = data_phenotype_resistant[['Unnamed: 0', 'EGF_ON_Proliferation']]
#     data_phenotype_resistant_filtered = data_phenotype_resistant_filtered.sort_values(by="EGF_ON_Proliferation", ascending=False)
#     data_phenotype_sensitive = data_phenotype_patients[data_phenotype_patients['Unnamed: 0'].isin(top_sensitive_ids)]
#     data_phenotype_sensitive_filtered = data_phenotype_sensitive[['Unnamed: 0', 'EGF_ON_Proliferation']]
#     data_phenotype_sensitive_filtered = data_phenotype_sensitive_filtered.sort_values(by="EGF_ON_Proliferation", ascending=False)

#     proliferating_patients = ['SIDM00719', 'SIDM00119']

#     other_patients = [x for x in patients_ids if x not in ['SIDM00719', 'SIDM00119']]

#     # extract gene expression data 
#     genes_data = (
#         pd.read_csv('../cancer_data/rnaseq_merged/rnaseq_merged_20250117.csv')
#         .query('model_id in @patients_ids')
#         .loc[:, ['model_id', 'gene_symbol', 'rsem_tpm']]
#     )

#     mean_gene_rsem = (
#         genes_data
#         .groupby(['model_id', 'gene_symbol'])['rsem_tpm']
#         .mean()
#         .reset_index()
#     )
#     # compute mean and sd of each gene in the two groups
#     genes_names = mean_gene_rsem['gene_symbol'].unique().tolist()
#     genes_std_results = pd.DataFrame(index=genes_names, columns=['Group Proliferation', 'Group Not Proliferation'])

#     for gene in genes_names:
#         # Filter by gene and by patient group
#         group_proliferation = mean_gene_rsem[
#             (mean_gene_rsem['gene_symbol'] == gene) &
#             (mean_gene_rsem['model_id'].isin(proliferating_patients))
#         ]['rsem_tpm'].std()
        
#         group_not_proliferation = mean_gene_rsem[
#             (mean_gene_rsem['gene_symbol'] == gene) &
#             (mean_gene_rsem['model_id'].isin(other_patients))
#         ]['rsem_tpm'].std()
        
#         # Assign to the results DataFrame
#         genes_std_results.at[gene, 'Group Proliferation'] = group_proliferation
#         genes_std_results.at[gene, 'Group Not Proliferation'] = group_not_proliferation

#     return(genes_std_results)




# data_phenotype_patients = pd.read_csv('personalized_boolean/patients_phenot_table.csv')

# tissue_interest = 'Lung'
# drug_interest = 'Avagacestat'
# top_sensitive_ids, top_resistant_ids = get_patients_top_10(tissue_interest, drug_interest)
# patients_ids = top_sensitive_ids + top_resistant_ids
# #genes_mean_results = compute_genes_mean_signature(data_phenotype_patients,top_resistant_ids, top_sensitive_ids)
# #genes_mean_results.to_csv('personalized_boolean/genes_signature_mean_results.csv', index=True)



# genes_std_results = compute_genes_std_signature(data_phenotype_patients, top_resistant_ids, top_sensitive_ids)
# genes_std_results.to_csv('personalized_boolean/genes_signature_std_results.csv', index=True)




