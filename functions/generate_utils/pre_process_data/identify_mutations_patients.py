import pandas as pd



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
