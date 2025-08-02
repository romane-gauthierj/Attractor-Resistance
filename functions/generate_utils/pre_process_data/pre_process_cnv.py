import pandas as pd
import sys
import os


def preprocess_cnv(patients_ids, cnv_data, all_montagud_nodes, synonyms_to_nodes_dict):
    # filter out Neural Category
    cnv_data_filtered = cnv_data[cnv_data["cn_category"] != "Neutral"]

    # keep only the patients id and montagud nodes
    cnv_data_filtered = cnv_data_filtered[
        cnv_data_filtered["symbol"].isin(all_montagud_nodes)
    ]
    cnv_data_filtered = cnv_data_filtered[
        cnv_data_filtered["model_id"].isin(patients_ids)
    ]
    cnv_data_filtered = cnv_data_filtered[
        ["model_id", "symbol", "total_copy_number", "cn_category"]
    ]
    cnv_data_filtered.rename(columns={"symbol": "gene_symbol"}, inplace=True)

    cnv_data_filtered = cnv_data_filtered[
        ~cnv_data_filtered["total_copy_number"].isna()
    ]

    # remplace the cnv gene symbol column names by its synonyms in the proteins_synonyms_maps dictionary

    # apply the synonym mapping 
    for gene in cnv_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict:
            cnv_data_filtered.loc[cnv_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]

    return cnv_data_filtered



# # Add the parent directory to sys.path
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# from identification_patients.get_patients_sens_res import get_patients


# annotations_models = pd.read_csv('../data/model_list_20250407.csv')
# mutations_data = pd.read_csv('../data/mutations_all_20250318.csv')
# drug_data = pd.read_csv('../data/drug_sensitivity.csv')

# montagud_data = (
#     pd.read_csv('../data/Montagud_inter_nodes_data.csv', header=1)
#     .loc[:, ['Target node', 'Interaction type', 'Source']])
# montagud_nodes = list(set(montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()))
# montagud_nodes = [node for node in montagud_nodes if node != '0/1']
# montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]
# montagud_nodes.append('KRAS')
# to_remove = ['RAS', 'FUSED_EVENT', 'NKX3_1', 'SPOP', 'AR_ERG']

# montagud_nodes = [node for node in montagud_nodes if node not in to_remove]
# montagud_nodes = list(set(montagud_nodes))
# # chose folder where we want all the personalized boolean models and associated results saved
# drug_interest = 'Refametinib' #Pictilisib, 'Avagacestat' AZD8931
# tissue_interest = 'Lung'
# tissue_remove = 'Haematopoietic and Lymphoid'


# top_resistant_ids, top_sensitive_ids, drug_tissue_data= get_patients(20, drug_data, annotations_models, drug_interest, tissue_interest, tissue_remove)
# patients_ids = top_sensitive_ids + top_resistant_ids
# cnv_data = pd.read_csv('../data/cellmodel_data/cnv_summary_20250207.csv')


# def preprocess_cnv(cnv_data, montagud_nodes, patients_ids, synonyms_maps):
#     # filter out Neural Category
#     cnv_data_filtered = cnv_data[cnv_data["cn_category"] != "Neutral"]

#     # keep only the patients id and montagud nodes
#     cnv_data_filtered = cnv_data_filtered[
#         cnv_data_filtered["symbol"].isin(montagud_nodes)
#     ]
#     cnv_data_filtered = cnv_data_filtered[
#         cnv_data_filtered["model_id"].isin(patients_ids)
#     ]
#     cnv_data_filtered = cnv_data_filtered[
#         ["model_id", "symbol", "total_copy_number", "cn_category"]
#     ]
#     cnv_data_filtered.rename(columns={"symbol": "gene_symbol"}, inplace=True)

#     cnv_data_filtered = cnv_data_filtered[
#         ~cnv_data_filtered["total_copy_number"].isna()
#     ]

#     # remplace the cnv gene symbol column names by its synonyms in the proteins_synonyms_maps dictionary
#     cnv_data_filtered["gene_symbol"] = cnv_data_filtered["gene_symbol"].apply(
#         lambda x: synonyms_maps.get(x, x)
#     )

#     return cnv_data_filtered


