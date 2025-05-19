
import pandas as pd
import os
import re
import numpy as np




def get_patients_top_10(drug_data, annotations_models, drug_interest, tissue_interest=None):
    # extract the top 100 and bottom 100 patients resitant to the drug
    if tissue_interest is not None:
        tissue_interest = tissue_interest.upper()
        annotations_models['tissue'] = annotations_models['tissue'].str.upper()
        annotations_models = annotations_models[annotations_models['tissue'] == tissue_interest]
    
    models_id = annotations_models['model_id'].tolist()
    drug_data_filtered = drug_data[drug_data['SANGER_MODEL_ID'].isin(models_id)]
    drug_data_filtered = drug_data_filtered[['SANGER_MODEL_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'AUC', 'LN_IC50', 'Z_SCORE']]
    drug_data_filtered = drug_data_filtered[drug_data_filtered['DRUG_NAME'] == drug_interest]
    grouped_drug_data_filtered = drug_data_filtered.groupby('SANGER_MODEL_ID')['Z_SCORE'].mean().reset_index()

    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(by='Z_SCORE')
    # print(drug_data_filtered_ranked)
    
    top_sensitive = drug_data_filtered_ranked.iloc[0:100, :]

    top_sensitive_ids = top_sensitive['SANGER_MODEL_ID'].tolist()
    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(by='Z_SCORE', ascending=False)
    top_resistant = drug_data_filtered_ranked.iloc[0:100, :]

    top_resistant_ids = top_resistant['SANGER_MODEL_ID'].tolist()
    top_resistant_ids = list(set(top_resistant_ids))
    top_sensitive_ids = list(set(top_sensitive_ids))

    return top_resistant_ids, top_sensitive_ids






# def get_patients_top_10(tissue_interest, drug_interest):
#     drug_data = pd.read_csv('../cancer_data/drug_sensitivity.csv')
#     annotations_models = pd.read_csv('../cancer_data/model_list_20250407.csv')
#     annotations_models['tissue'] = annotations_models['tissue'].str.upper()

#     #annotations_models = annotations_models[annotations_models['cancer_type_detail'] == 'Lung Adenocarcinoma']
#     tissue_interest = tissue_interest.upper()
#     annotations_models = annotations_models[annotations_models['tissue'] == tissue_interest]
#     models_id = annotations_models['model_id'].tolist()
#     drug_data_filtered = drug_data[drug_data['SANGER_MODEL_ID'].isin(models_id)]
#     drug_data_filtered = drug_data_filtered[['SANGER_MODEL_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'AUC', 'LN_IC50', 'Z_SCORE']]
#     drug_data_filtered = drug_data_filtered[drug_data_filtered['DRUG_NAME'] == drug_interest]
#     grouped_drug_data_filtered = drug_data_filtered.groupby('SANGER_MODEL_ID')['Z_SCORE'].mean().reset_index()
#     drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(by='Z_SCORE')
#     top_sensitive = drug_data_filtered_ranked.iloc[1:11, :]
#     top_sensitive_ids = top_sensitive['SANGER_MODEL_ID'].tolist()
#     drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(by='Z_SCORE', ascending=False)
#     top_resistant = drug_data_filtered_ranked.iloc[1:11, :]
#     top_resistant_ids = top_resistant['SANGER_MODEL_ID'].tolist()
#     return top_sensitive_ids, top_resistant_ids






# def get_patients_zscore_strict(tissue_interest, drug_interest):
#     drug_data = pd.read_csv('../cancer_data/drug_sensitivity.csv')
#     # identify the drug data for LUNG cancer
#     annotations_models = pd.read_csv('../cancer_data/model_list_20250407.csv')
#     annotations_models['tissue'] = annotations_models['tissue'].str.upper()

#     #annotations_models = annotations_models[annotations_models['cancer_type_detail'] == 'Lung Adenocarcinoma']
#     annotations_models = annotations_models[annotations_models['tissue'] == tissue_interest]
#     models_id = annotations_models['model_id'].tolist()
#     drug_data_filtered = drug_data[drug_data['SANGER_MODEL_ID'].isin(models_id)]
#     drug_data_filtered = drug_data_filtered[['SANGER_MODEL_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'AUC', 'LN_IC50', 'Z_SCORE']]

#     # Identify the drug for which we explore drivers of resistance
#     drug_counts = drug_data_filtered['DRUG_NAME'].value_counts()
#     # print(drug_counts)
#     # select the drug with the most data available
#     #print(drug_data_filtered[drug_data_filtered['DRUG_NAME'] == 'Erlotinib']) # 66
#     #print(drug_data_filtered[drug_data_filtered['DRUG_NAME'] == 'Avagacestat']) #357
#     #print(drug_data_filtered[drug_data_filtered['DRUG_NAME'] == 'Gefitinib']) #181

#     # select Avagacestat
#     drug_data_filtered = drug_data_filtered[drug_data_filtered['DRUG_NAME'] == drug_interest]
#     # create new column for drug sensitivity based on LN IC50 mean 
#     #ic50_mean_avagacestat = drug_data_filtered['LN_IC50'].mean()
#     #drug_data_filtered['drug_sensitivity'] = drug_data_filtered['LN_IC50'].apply(
#     # lambda x: 'Sensitivity' if x < ic50_mean_avagacestat else 'Resistant')

#     drug_data_filtered['drug_sensitivity'] = drug_data_filtered['Z_SCORE'].apply(
#         lambda x: 'Sensitivity' if x < -2 else 'Resistant' if x > 2 else 'Intermediate'
#     )

#     # print(drug_data_filtered['drug_sensitivity'].value_counts())
#     drug_data_filtered = drug_data_filtered[drug_data_filtered['drug_sensitivity']!='Intermediate']

#     patients_drug_ids = drug_data_filtered['SANGER_MODEL_ID'].tolist()

#     # among the sanger model ID, select the one present in the three datasets
#     patients_ids = table_proteins_patients.index.tolist() + table_cnv_patients.index.tolist() + table_genes_patients.index.tolist()
#     patients_ids = list(set(patients_ids))

#     # Common id bw id from drug and patients
#     common_elements = list(set(patients_drug_ids) & set(patients_ids))
#     common_elements = list(set(common_elements))

#     # keep only the id in common 
#     drug_data_filtered = drug_data_filtered[drug_data_filtered['SANGER_MODEL_ID'].isin(common_elements)]

#     drug_data_filtered.to_csv('drug_sensitivity_clinical.csv', index=True)
#     resistant_patients = drug_data_filtered[drug_data_filtered['drug_sensitivity'] == 'Resistant']
#     sensitive_patients = drug_data_filtered[drug_data_filtered['drug_sensitivity'] == 'Sensitivity']
    
#     sensitive_patients_ids = list(set(sensitive_patients['SANGER_MODEL_ID'].tolist()))
#     resistant_patients_ids = list(set(resistant_patients['SANGER_MODEL_ID'].tolist()))
#     return sensitive_patients_ids, resistant_patients_ids





    ## TEST 
    # print(grouped_drug_data_filtered)
    
    # Convert Z_SCORE column to float, keeping the whole DataFrame
    # grouped_drug_data_filtered['Z_SCORE'] = pd.to_numeric(grouped_drug_data_filtered['Z_SCORE'], errors='coerce')

    # Now safely compute from the full DataFrame
    # df = grouped_drug_data_filtered

    # return {
    #     "name": drug_interest,
    #     "<-1.8": float((df['Z_SCORE'] <= -1.8).sum()),
    #     ">1.8": float((df['Z_SCORE'] >= 1.8).sum()),
    #     "mean": float(df['Z_SCORE'].mean()),
    #     "std": float(df['Z_SCORE'].std()),
    #     "abs_zscore": float(df['Z_SCORE'].abs().mean())
    # }

    ## END