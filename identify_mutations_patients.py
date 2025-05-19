import pandas as pd
from get_patients_sens_res import get_patients_top_10



def identif_mutations_kras_egfr(mutations_data, patients_ids):
    mutations_kras_ids_list = []
    mutations_egfr_ids_list = []

    mutations_lung = mutations_data[mutations_data['gene_symbol'].isin(['KRAS', 'EGFR'])]
    mutations_lung = mutations_lung[mutations_lung['model_id'].isin(patients_ids)]
    
    # select only the missense mutation as we know these mutations drive overactivation
    mutations_lung = mutations_lung[mutations_lung['effect']== 'missense']
    mutations_kras_ids = mutations_lung[mutations_lung['gene_symbol'] == 'KRAS']['model_id'].tolist()
    mutations_kras_ids_list.extend(mutations_kras_ids)


    mutations_egfr_ids = mutations_lung[mutations_lung['gene_symbol'] == 'EGFR']['model_id'].tolist()
    mutations_egfr_ids_list.extend(mutations_egfr_ids)

    return mutations_kras_ids_list, mutations_egfr_ids_list


# tissue_interest = 'Lung'
# drug_interest = 'AZD8931' #'Avagacestat'
# annotations_models = pd.read_csv('../data/model_list_20250407.csv')
# mutations_data = pd.read_csv('../data/mutations_all_20250318.csv')
# drug_data = pd.read_csv('../data/drug_sensitivity.csv')
# mutations_data = pd.read_csv('../data/mutations_all_20250318.csv')

# top_resistant_ids, top_sensitive_ids= get_patients_top_10(drug_data, annotations_models, tissue_interest, drug_interest)
# patients_ids = top_sensitive_ids + top_resistant_ids



# mutations_kras_ids_list, mutations_egfr_ids_list = identif_mutations_kras_egfr(mutations_data, patients_ids)
# print('patients with mutations KRAS:', mutations_kras_ids_list)
# print('patients with mutations EGFR:', mutations_egfr_ids_list)