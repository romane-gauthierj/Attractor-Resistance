# Pipeline to compute patients phenotypes 

import pandas as pd 
from get_patients_sens_res import get_patients_top_10
from create_generic_patients_cfgs import create_generic_patients_cfgs_bnds
from pre_process_genes import process_genes

from tailor_cfgs_patients_gene import personalized_patients_genes_cfgs
from MaBoSS_phenotype_distribution import compute_phenotypes_distribution, compute_mean_patients
# from pre_process_profiles_table_data_lung import create_genes_patients
# from identify_mutations_patients import identif_mutations_kras_egfr
from stats_proba import compute_mannwhitneyu_test_means
from tailor_bnd_mutations import personalized_patients_mutations_bnds
from boxplot_phenotype_V2 import create_boxplot
from create_phenotypes_patients_table import vizualise_table_phenotype_condition
from patients_ids_phenotype_table import create_table_patients_phenotypes
from genes_signature import compute_genes_mean_signature
import numpy as np



# Import data
annotations_models = pd.read_csv('../data/model_list_20250407.csv')
mutations_data = pd.read_csv('../data/mutations_all_20250318.csv')
drug_data = pd.read_csv('../data/drug_sensitivity.csv')

montagud_data = (
    pd.read_csv('../data/Montagud_inter_nodes_data.csv', header=1)
    .loc[:, ['Target node', 'Interaction type', 'Source']])
montagud_nodes = list(set(montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()))
rna_seq_data = pd.read_csv('../data/rnaseq_merged/rnaseq_merged_20250117.csv')
#genes_data_filtered = pd.read_csv('filtered_data/rna_seq_lung_clean.csv')




# chose folder where we want all the models and results saved 
folder='personalized_boolean_large_groups'



# Output directories
output_dir_resistant = f'{folder}/resistant_patient/generic_models'
output_dir_sec_resistant = f'{folder}/resistant_patient/personalized_boolean_modified'
output_dir_sensitive = f'{folder}/sensitive_patient/generic_models'
output_dir_sec_sensitive = f'{folder}/sensitive_patient/personalized_boolean_modified'
bnd_dir_res = f'{folder}/resistant_patient/personalized_boolean_modified/models_gene_expression'
bnd_dir_sens = f'{folder}/sensitive_patient/personalized_boolean_modified/models_gene_expression'





# ---------------- Step 1: Select cancer and drug of interest (tissue_interest, drug_interest)------------------------
# Pre-process genes data 
# tissue_interest = 'Lung'
# top_resistant_ids, top_sensitive_ids= get_patients_top_10(drug_data, annotations_models, drug_interest, tissue_interest)
drug_interest = 'AZD8931' #'Avagacestat' AZD8931
tissue_remove = 'Haematopoietic and Lymphoid'
top_resistant_ids, top_sensitive_ids, drug_tissue_data= get_patients_top_10(drug_data, annotations_models, drug_interest)
patients_ids = top_sensitive_ids + top_resistant_ids

# check if KRAS is also in the montagud_data
rna_seq_data_filtered = process_genes(patients_ids, montagud_data, rna_seq_data)



# ---------------- Step 2: Create generic boolean networks with the sensitive and resistant ID names ---------------------
create_generic_patients_cfgs_bnds(folder, top_resistant_ids, top_sensitive_ids, drug_interest)



# ----------------Step 3: Personalize the cfg files with genes/ proteins ---------------------------------
personalized_patients_genes_cfgs(montagud_data, output_dir_resistant, output_dir_sec_resistant, patients_ids, rna_seq_data_filtered, drug_interest)
personalized_patients_genes_cfgs(montagud_data, output_dir_sensitive, output_dir_sec_sensitive, patients_ids, rna_seq_data_filtered, drug_interest)




# ---------------- Step 4: Identification of which patients id have KRAS or EGFR mutation  ---------------------------------
# #         personalize the bnd files with the mutations common to Lung (KRAS/ EGFR): 
personalized_patients_mutations_bnds(mutations_data,patients_ids,bnd_dir_res, drug_interest)
personalized_patients_mutations_bnds(mutations_data,patients_ids,bnd_dir_sens, drug_interest)



# ---------------- Step 5: compute the phenotype distribution 
dic_patient_resistant =f'{folder}/resistant_patient/personalized_boolean_modified/models_gene_expression'
dic_patient_sensitive =f'{folder}/sensitive_patient/personalized_boolean_modified/models_gene_expression'
inputs_list = ['EGF', 'FGF', 'TGFb', 'Nutrients', 'Hypoxia', 'Acidosis', 'Androgen', 'TNFalpha', 'Carcinogen']

patient_res_data_dict = compute_phenotypes_distribution(folder, dic_patient_resistant, inputs_list, 'resistant', drug_interest)
patient_sens_data_dict = compute_phenotypes_distribution(folder, dic_patient_sensitive, inputs_list, 'sensitive', drug_interest)
patients_res_df_mean, patients_res_df_std, stats_results_data_res_df = compute_mean_patients(patient_res_data_dict)
patients_sens_df_mean, patients_sens_df_std, stats_results_data_sens_df = compute_mean_patients(patient_sens_data_dict)



patients_res_df_mean.to_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_df_mean_{drug_interest}.csv', index=True)
patients_res_df_std.to_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_df_std_{drug_interest}.csv', index=True)
stats_results_data_res_df.to_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_values_stats_{drug_interest}.csv', index=True)
stats_results_data_sens_df.to_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_values_stats_{drug_interest}.csv', index=True)
patients_sens_df_mean.to_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_df_mean_{drug_interest}.csv', index=True)
patients_sens_df_std.to_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_df_std_{drug_interest}.csv', index=True)



# ----------------  Step 6: Compute stats test between two mean datasets---------------------- 
patient_res_stats_values = pd.read_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_values_stats_{drug_interest}.csv')
patient_sens_stats_values = pd.read_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_values_stats_{drug_interest}.csv')
compute_mannwhitneyu_test_means(folder,patient_res_stats_values, patient_sens_stats_values, drug_interest)




# ----------------  Step 7: Vizualise the boxplot of phenotype distribution output--------------------
patient_res_values = pd.read_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_values_stats_{drug_interest}.csv')
patient_sens_values = pd.read_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_values_stats_{drug_interest}.csv')
data_greater_side = pd.read_csv(f'{folder}/sensitive_resistant_results/p_values_df_mannwhitneyu_greater_sign_{drug_interest}.csv')
create_boxplot(folder, patient_res_values, patient_sens_values, data_greater_side)


# ----------------- Step 8: create table of patients with conditions- phenotype
dir_res_data = f'{folder}/resistant_results/only_gene_expression/single_input_on/phenotype_distribution_patients'
dir_sens_data = f'{folder}/sensitive_results/only_gene_expression/single_input_on/phenotype_distribution_patients'
patients_phenot_table = create_table_patients_phenotypes(folder, dir_res_data, dir_sens_data)






# -----------------  Step 9: Create heatmap figure 

patient_resistant_mean = pd.read_csv(f'{folder}/resistant_results/only_gene_expression/single_input_on/patients_resistant_df_mean_{drug_interest}.csv')
patient_sensitive_mean = pd.read_csv(f'{folder}/sensitive_results/only_gene_expression/single_input_on/patients_sensitive_df_mean_{drug_interest}.csv')
vizualise_table_phenotype_condition(folder, patient_resistant_mean, patient_sensitive_mean)






# -----------------  Step 10: Identify genes differently expressed in the patients with high 

patients_phenot_table = pd.read_csv(f'{folder}/sensitive_resistant_results/patients_phenot_table.csv')
genes_stats_results_metast_TGFb = compute_genes_mean_signature(folder, montagud_nodes, 'Metastasis', 'TGFb', patients_phenot_table, top_resistant_ids, top_sensitive_ids)
genes_stats_results_prolif_egf = compute_genes_mean_signature(folder, montagud_nodes, 'Proliferation', 'EGF', patients_phenot_table, top_resistant_ids, top_sensitive_ids)






# -----------------  Step 11: check there is not correlation between phenotype distribution and cancer type ----------------------------
patients_phenot_table['SANGER_MODEL_ID'] = patients_phenot_table['Unnamed: 0'].str.split('_').str[0]
conditions = [
    patients_phenot_table['SANGER_MODEL_ID'].isin(top_resistant_ids),
    patients_phenot_table['SANGER_MODEL_ID'].isin(top_sensitive_ids)
    ]
choices = ['Resistant', 'Sensitive']
patients_phenot_table.loc[:,'Drug status'] = np.select(conditions, choices, default = '')



ids_tissue_data = drug_tissue_data[['SANGER_MODEL_ID', 'tissue']]
ids_tissue_data = ids_tissue_data.drop_duplicates(subset='SANGER_MODEL_ID')


# merge tissues and model id 
patients_phenot_table = pd.merge(patients_phenot_table, ids_tissue_data, on = 'SANGER_MODEL_ID')
print(patients_phenot_table)



# look the number of each cancer for the condition-phenotype of interest
condition = 'TGFb'
phenotype = 'Metastasis'



# resistant group changes according to what is the condition and the phenotype
# group_proliferation_resistant: group with high phenotype 

group_phenotype_resistant = patients_phenot_table[
    (patients_phenot_table['Drug status'] == 'Resistant') & 
    (patients_phenot_table[f'{condition}_ON_{phenotype}'] >= 0.1)
]




print(group_phenotype_resistant['tissue'].value_counts()) # EGF- proliferation: 4 lung, 1 breast, 1 haematopoetic
                                                          # TGFb- Metastasis: 21 haemato, 2 skin, 2 breast, 1 lung, 1 large intestine, 1 endom, 1 liver








# ???????
# how to check that specific cancers are more present in one category???




# check what drug is the best to keep (the one with most resistant and sensitive)
# results = {}
# drug_interests = drug_data['DRUG_NAME'].unique().tolist()
# #print(drug_interests)
# for drug_interest in drug_interests:
#     results[drug_interest] = get_patients_top_10(drug_data, annotations_models, tissue_interest, drug_interest)

# drug_interest: {
#     "name": drug_interest,
#     "<-1.5": float((df['Z_SCORE'] < -1.5).sum()),
#     ">1.5": float((df['Z_SCORE'] > 1.5).sum()),
#     "mean": float(df['Z_SCORE'].mean()),
#     "std": float(df['Z_SCORE'].std()),
#     "abs_zscore": float(df['Z_SCORE'].abs().mean()),
# }

# list_results = results.values()
# pd_results = pd.DataFrame(list_results)
# pd_results.to_csv(f'{folder}/drug_analysis.csv")
