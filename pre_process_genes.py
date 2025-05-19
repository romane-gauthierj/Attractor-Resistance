# Cell model passport - RNA seq expression

import numpy as np
import pandas as pd
from get_patients_sens_res import get_patients_top_10
import math


import os


def classify_expression(z):
        if z > 1:
            return 'high'
        elif z < -1:
            return 'low'
        else:
            return 'normal' 




def process_genes(patients_ids, montagud_data, rna_seq_data):
    # filter only data with patient id and the nodes of the montagud_data
    montagud_data['Target node'] = montagud_data['Target node'].astype(str)
    montagud_data['Source'] = montagud_data['Source'].astype(str)
    montagud_nodes = montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()
    montagud_nodes = list(set(montagud_nodes))
    montagud_nodes_upper = [node.upper() for node in montagud_nodes]

    rna_seq_data['gene_symbol_upper'] = rna_seq_data['gene_symbol'].str.upper()
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['gene_symbol_upper'].isin(montagud_nodes_upper)]

    rna_seq_data_filtered = rna_seq_data_filtered[rna_seq_data_filtered['model_id'].isin(patients_ids)]
    rna_seq_data_filtered['z_score'] = rna_seq_data_filtered.groupby('gene_symbol')['rsem_tpm'].transform(
            lambda x: (x - x.mean()) / x.std()
        )

    rna_seq_data_filtered['gene_expression_level'] = rna_seq_data_filtered['z_score'].apply(classify_expression)

    rna_seq_data_filtered = rna_seq_data_filtered[['model_id', 'gene_symbol', 'gene_expression_level']]
    rna_seq_data_filtered.rename(columns={'gene_symbol': 'gene_name'}, inplace=True)
    rna_seq_data_filtered = rna_seq_data_filtered[rna_seq_data_filtered['gene_expression_level'].isin(['low', 'high'])]
    
    table_rna_seq_patients = rna_seq_data_filtered.pivot_table(
            index='model_id',
            columns='gene_expression_level',
            values='gene_name',
            aggfunc=lambda x: ', '.join(sorted(set(x)))
        ).fillna('-')


    table_rna_seq_patients = table_rna_seq_patients.rename(columns={
        'low': 'Low Gene Expression',
        'high': 'High Gene Expression'
    })

    return(table_rna_seq_patients)












