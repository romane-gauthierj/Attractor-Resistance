import pandas as pd 
import sys
import os

# DeSeq mutations preprocess pipeline

def preprocess_mutations_ds(patients_ids, models_depmap_annotation,mutations_data_depseq, montagud_nodes):
    # preprocess mutations from depseq (ids - can retrieve ids of cell model passport with cell model annotation of cell passport)
    # keep only montagud nodes, keep only the TSG/ Oncogene with impact 

    models_annotation = models_depmap_annotation[['model_id', 'BROAD_ID']]
    models_annotation.rename(columns={'model_id': 'cell model id'}, inplace=True)
    models_annotation.rename(columns={'BROAD_ID': 'depmap id'}, inplace=True)
    models_annotation = models_annotation.dropna(axis=0)
    models_annotation = models_annotation[models_annotation['cell model id'].isin(patients_ids)]

    models_annot_ids = list(models_annotation['depmap id'])
    mutations_data_filtered = mutations_data_depseq[['VariantType', 'VariantInfo', 'HugoSymbol', 'OncogeneHighImpact', 'TumorSuppressorHighImpact', 'ModelID']]
    mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['ModelID'].isin(models_annot_ids)]
    mutations_data_filtered = mutations_data_filtered[(mutations_data_filtered['OncogeneHighImpact'] == True) | (mutations_data_filtered['TumorSuppressorHighImpact'] == True)]
    mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['HugoSymbol'].isin(montagud_nodes)]
    mutations_data_filtered.rename(columns={'ModelID': 'depmap id'}, inplace=True)
    mutations_data_filtered_combined = pd.merge(mutations_data_filtered,models_annotation, on='depmap id')
    mutations_data_filtered_combined = mutations_data_filtered_combined.set_index('cell model id')


    return mutations_data_filtered_combined
