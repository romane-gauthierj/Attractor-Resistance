import pandas as pd
import os
import re

import pandas as pd
import sys
import os

# # Add the parent directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from identification_patients.get_patients_sens_res import get_patients
from pre_process_data.pre_process_cnv import preprocess_cnv


annotations_models = pd.read_csv('../data/model_list_20250407.csv')
mutations_data = pd.read_csv('../data/mutations_all_20250318.csv')
drug_data = pd.read_csv('../data/drug_sensitivity.csv')

montagud_data = (
    pd.read_csv('../data/Montagud_inter_nodes_data.csv', header=1)
    .loc[:, ['Target node', 'Interaction type', 'Source']])
montagud_nodes = list(set(montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()))
montagud_nodes = [node for node in montagud_nodes if node != '0/1']
montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]
montagud_nodes.append('KRAS')
to_remove = ['RAS', 'FUSED_EVENT', 'NKX3_1', 'SPOP', 'AR_ERG']

montagud_nodes = [node for node in montagud_nodes if node not in to_remove]
montagud_nodes = list(set(montagud_nodes))
# chose folder where we want all the personalized boolean models and associated results saved 
drug_interest = 'Refametinib' #Pictilisib, 'Avagacestat' AZD8931 
tissue_interest = 'Lung'
tissue_remove = 'Haematopoietic and Lymphoid'
top_resistant_ids, top_sensitive_ids, drug_tissue_data= get_patients(20, drug_data, annotations_models, drug_interest, tissue_interest, tissue_remove)
patients_ids = top_sensitive_ids + top_resistant_ids
cnv_data = pd.read_csv('../data/cellmodel_data/cnv_summary_20250207.csv')

cnv_data_filtered = preprocess_cnv(cnv_data, montagud_nodes,patients_ids)
folder_pers_models=f'../models/personalized_boolean_{drug_interest}_{tissue_interest}'

bnd_dir_sens = f'{folder_pers_models}/sensitive_patient/personalized_boolean_modified/models_gene_expression'
bnd_dir_res = f'{folder_pers_models}/resistant_patient/personalized_boolean_modified/models_gene_expression'





gain_group = cnv_data_filtered[
    cnv_data_filtered['cn_category'].isin(['Amplification', 'Gain']) &
    (cnv_data_filtered['total_copy_number'] > 2.0)
]

loss_group = cnv_data_filtered[
    cnv_data_filtered['cn_category'].isin(['Deletion', 'Loss']) &
    (cnv_data_filtered['total_copy_number'] < 2.0)
]

gain_ids = set(gain_group['model_id'])
loss_ids = set(loss_group['model_id'])
patients_all = list(gain_ids | loss_ids)

for patient in patients_all:
    bnd_file = None
    files_res = [file for file in os.listdir(bnd_dir_res) if file.endswith('.bnd') and patient in file]
    if files_res:
        bnd_file = os.path.join(bnd_dir_res, files_res[0])
    else:
        files_sens = [file for file in os.listdir(bnd_dir_sens) if file.endswith('.bnd') and patient in file]
        if files_sens:
            bnd_file = os.path.join(bnd_dir_sens, files_sens[0])
    if bnd_file is None:
        print(f"No .bnd file found for patient: {patient}")
        continue

    genes = cnv_data_filtered[cnv_data_filtered['model_id'] == patient]['symbol']
    gene_list = genes.tolist()

    patient_id = os.path.splitext(os.path.basename(bnd_file))[0].replace(f'_{drug_interest}', '')
    with open(bnd_file, 'r') as file:
        content = file.read()

    modified_any = False
    for gene in gene_list:
        print(f"🔍 Processing patient {patient}, gene: {gene}")
        pattern = re.compile(
            rf'Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}',
            re.DOTALL
        )

        if patient in gain_ids and patient in loss_ids:
            print(f"Patient {patient} is in both gain and loss groups. Please review.")
            new_gene_block = f"""Node {gene} {{
    logic = 1;
    rate_up = 1;
    rate_down = 0;
}}"""
        elif patient in gain_ids:
            new_gene_block = f"""Node {gene} {{
    logic = 1;
    rate_up = 1;
    rate_down = 0;
}}"""
        elif patient in loss_ids:
            new_gene_block = f"""Node {gene} {{
    logic = 0;
    rate_up = 0;
    rate_down = 1;
}}"""
        else:
            print(f"Patient {patient} not found in gain or loss CNV groups.")
            continue

        gene_match = pattern.search(content)
        if gene_match:
            print(f"{gene} node found. Replacing...")
            content, n_subs = re.subn(pattern, new_gene_block, content)
            if n_subs > 0:
                modified_any = True
        else:
            print(f"No {gene} node found in file for patient {patient_id}")

    if modified_any:
        print(f"{patient_id}: CNV — nodes modified")
        with open(bnd_file, 'w') as file:
            file.write(content)
    else:
        print(f"{patient_id}: CNV — no nodes modified")



# check that also works with loss (worked for gain!!)