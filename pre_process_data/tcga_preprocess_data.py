import pandas as pd 
import numpy as np


#phenotype_data = pd.read_csv('../data/TCGA_data/tcga_gdc_pancan_phenotype.csv')
tumor_stage_data = pd.read_csv('../data/TCGA_data/TCGA_GDC-PANCAN_tumor_stage_phenotype.csv')


print(tumor_stage_data.head())
print(tumor_stage_data['diagnoses.tumor_stage'].value_counts())

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
print(patients_group_df)


patients_group_df = patients_group_df[patients_group_df['Tumor Stage'].notna() & (patients_group_df['Tumor Stage'] != '')]
print(patients_group_df)
print(patients_group_df['Tumor Stage'].value_counts())