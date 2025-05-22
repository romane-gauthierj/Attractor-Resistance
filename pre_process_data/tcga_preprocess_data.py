import pandas as pd 
import numpy as np


#phenotype_data = pd.read_csv('../data/TCGA_data/tcga_gdc_pancan_phenotype.csv')
tumor_stage_data = pd.read_csv('../data/TCGA_data/TCGA_GDC-PANCAN_tumor_stage_phenotype.csv')
genes_data = pd.read_csv('../data/TCGA_data/TCGA_GDC-PANCAN_genes.csv')




# print(tumor_stage_data.head())
# print(tumor_stage_data['diagnoses.tumor_stage'].value_counts())

patients_ids = list(set(tumor_stage_data['sample']))

patients_group_df = pd.DataFrame(index=patients_ids)

group_0 = ['stage 0']
group_1 = ['stage i', 'stage ia', 'stage ib']
group_2 = ['stage ii', 'stage iia', 'stage iib', 'stage iic']
group_3 = ['stage iii', 'stage iiia', 'stage iiib', 'stage iiic']
group_4 = ['stage iv', 'stage iva', 'stage ivb', 'stage ivc']

conditions = [
    tumor_stage_data['diagnoses.tumor_stage'].isin(group_0),
    tumor_stage_data['diagnoses.tumor_stage'].isin(group_1),
    tumor_stage_data['diagnoses.tumor_stage'].isin(group_2),
    tumor_stage_data['diagnoses.tumor_stage'].isin(group_3),
    tumor_stage_data['diagnoses.tumor_stage'].isin(group_4) 
    ]
choices = ['Local Cancer', 'Early Stage', 'Larger Tumor', 'Advanced Local Speed', 'Metastatic']

patients_group_df.loc[:,'Tumor Stage'] = np.select(conditions, choices, default = '')
patients_group_df = patients_group_df[patients_group_df['Tumor Stage'].notna() & (patients_group_df['Tumor Stage'] != '')]
# print(patients_group_df['Tumor Stage'].value_counts())
patients_ids = patients_group_df.index.tolist()
print(patients_ids)

print(patients_group_df.head())
print(genes_data)
# filter genes to only keep the patients_ids info
genes_data_filtered = genes_data[patients_ids]
print(genes_data_filtered)
