import pandas as pd 
import sys
import os
import numpy as np

# DeSeq mutations preprocess pipeline


def pre_process_mutations(mutations_data, patients_ids, all_montagud_nodes, onco_tsg_data,synonyms_to_nodes_dict):
    mutations_data_filtered = mutations_data[mutations_data['model_id'].isin(patients_ids)]
    mutations_data_filtered = mutations_data_filtered[['gene_symbol', 'effect', 'model_id']]
    mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['gene_symbol'].isin(all_montagud_nodes)]

    onco_tsg_data_filtered = onco_tsg_data[['Hugo Symbol', 'Is Oncogene', 'Is Tumor Suppressor Gene']]
    onco_tsg_data_filtered = onco_tsg_data_filtered[onco_tsg_data_filtered['Hugo Symbol'].isin(all_montagud_nodes)]


    oncogenes_data = onco_tsg_data_filtered[(onco_tsg_data_filtered['Is Oncogene'] == 'Yes') & (onco_tsg_data_filtered['Is Tumor Suppressor Gene'] == 'No')]
    oncogenes_list = list(set(oncogenes_data['Hugo Symbol']))

    tsg_data = onco_tsg_data_filtered[(onco_tsg_data_filtered['Is Oncogene'] == 'No') & (onco_tsg_data_filtered['Is Tumor Suppressor Gene'] == 'Yes')]
    tsg_list = list(set(tsg_data['Hugo Symbol']))

    mutations_data_filtered = mutations_data_filtered[['gene_symbol', 'effect', 'model_id']]

    mutations_data_filtered['oncogene_tsg'] = np.where(
    mutations_data_filtered['gene_symbol'].isin(oncogenes_list), 'onco',
    np.where(mutations_data_filtered['gene_symbol'].isin(tsg_list), 'tsg', 'unknown')
    )

    inactivating_mutations= ['nonsense', 
    'frameshift',            
    'ess_splice',          
    'start_lost',           
    'cds_disrupted',        
    '5prime_UTR_ess_splice', 
    '3prime_UTR_ess_splice',
    'stop_lost']  

    mutations_data_filtered['inactivating_status'] = mutations_data_filtered['effect'].isin(inactivating_mutations).astype(int)



    mutations_data_filtered['effect'] = np.where(
        mutations_data_filtered['inactivating_status'] == 1, 'deletion',
        np.where(
            (mutations_data_filtered['inactivating_status'] == 0) & 
            (mutations_data_filtered['oncogene_tsg'] == 'onco'), 'amplification', 
            'unknown'
        )
    )

    mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['effect']!='unknown']

    mutations_data_filtered = mutations_data_filtered[['gene_symbol', 'model_id', 'effect']]
    mutations_data_filtered


    # apply the synonym mapping 
    for gene in mutations_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict:
            mutations_data_filtered.loc[mutations_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]


    return mutations_data_filtered