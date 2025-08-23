import pandas as pd
import numpy as np
import os
import re
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import spearmanr
import gseapy as gp

from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import maboss
from pathlib import Path

import seaborn as sns
import logging

logger = logging.getLogger(__name__)

from functions.generate_utils.pre_process_data.pre_process_genes import (
 identify_genes_distribution, compute_multi_distrib_normalization,
)


def process_montagud_nodes_synonyms(montagud_synonyms_data):

    montagud_node_synonyms = montagud_synonyms_data.dropna(subset=['Node'])

    #Strip whitespace from columns BEFORE processing
    montagud_node_synonyms['Node'] = montagud_node_synonyms['Node'].str.strip()
    if 'HGNC.symbols' in montagud_node_synonyms.columns:
        montagud_node_synonyms['HGNC.symbols'] = montagud_node_synonyms['HGNC.symbols'].str.strip()


     # Create multiple rows at once
    new_rows = pd.DataFrame({
        'Node': ['cFLAR', 'eEF2', 'eEF2K', 'Rheb'],  
        'HGNC.symbols': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'], 
        'unique': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'] 
    })

    montagud_node_synonyms = pd.concat([montagud_node_synonyms, new_rows], ignore_index=True)


    montagud_node_synonyms = montagud_node_synonyms.rename(columns={'HGNC.symbols': 'Node_synonyms'})
    
    montagud_node_synonyms['Node_synonyms'] = montagud_node_synonyms['Node_synonyms'].str.split(',')
    montagud_node_synonyms = montagud_node_synonyms.explode('Node_synonyms')

     # Strip whitespace from Node column
    montagud_node_synonyms['Node_synonyms'] = montagud_node_synonyms['Node_synonyms'].str.strip()
    montagud_node_synonyms['Node'] = montagud_node_synonyms['Node'].str.strip()


    montagud_node_synonyms = montagud_node_synonyms[montagud_node_synonyms['Node_synonyms'] != '']
    montagud_node_synonyms = montagud_node_synonyms[montagud_node_synonyms['Node'] != '']

    

    synonyms_to_nodes_dict = montagud_node_synonyms.set_index('Node_synonyms')['Node'].to_dict()

    return montagud_node_synonyms, synonyms_to_nodes_dict





def process_montagud_nodes(
    montagud_original_data_df, montagud_node_synonyms,
):

    # Create list of genes of interest (in Montagud data)
    # montagud_node_model are the nodes of the model 
    montagud_node_model = list(
        set(montagud_original_data_df["Target node"].tolist() + montagud_original_data_df["Source"].tolist())
    )
    montagud_node_model = [node for node in montagud_node_model if node != "0/1"]
    all_montagud_nodes = list(set(montagud_node_synonyms['Node_synonyms'])) + list(set(montagud_node_synonyms['Node'])) + list(set(montagud_node_model))

    return montagud_node_model, all_montagud_nodes






def process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict):
    gene_info_cols = ['Hugo_Symbol', 'Entrez_Gene_Id'] 
    patient_cols = [col for col in rna_seq_data.columns if col in patients_ids]
    
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['Hugo_Symbol'].isin(all_montagud_nodes)]


    # Strip whitespace from Hugo_Symbol BEFORE synonym mapping
    rna_seq_data_filtered['Hugo_Symbol'] = rna_seq_data_filtered['Hugo_Symbol'].str.strip()


  # Special cases handling: Create both mTORC1 and mTORC2 from MTOR data
    syn_dict = {'MTOR': ['mTORC1', 'mTORC2'], 'MYC': ['MYC', 'MYC_MAX'], 'PIK3CA': ['PI3K', 'PIP3'], 'LDHA': ['LDHA', 'Lactic_acid'], 'ERG': ['AR_ERG', 'ERG']}
    # Get all keys
    list_genes_duplicates = syn_dict.keys()

    # apply the synonym mapping 
    for gene in rna_seq_data_filtered['Hugo_Symbol'].unique():
        if gene in synonyms_to_nodes_dict and gene not in list_genes_duplicates:
            rna_seq_data_filtered.loc[rna_seq_data_filtered['Hugo_Symbol'] == gene, 'Hugo_Symbol'] = synonyms_to_nodes_dict[gene]


    for duplicate_gene in list_genes_duplicates:
        gene_duplicate_data = rna_seq_data_filtered[rna_seq_data_filtered['Hugo_Symbol'] == duplicate_gene].copy()
        if not gene_duplicate_data.empty:
            # Create mTORC1 data
            gene_duplicate_1_data = gene_duplicate_data.copy()
            gene_duplicate_1_data['Hugo_Symbol'] = syn_dict[duplicate_gene][0]
            
            # Create mTORC2 data  
            gene_duplicate_2_data = gene_duplicate_data.copy()
            gene_duplicate_2_data['Hugo_Symbol'] = syn_dict[duplicate_gene][1]
            
            # Remove original MTOR and add both complexes
            rna_seq_data_filtered = rna_seq_data_filtered[rna_seq_data_filtered['Hugo_Symbol'] != duplicate_gene]
            rna_seq_data_filtered = pd.concat([rna_seq_data_filtered, gene_duplicate_1_data, gene_duplicate_2_data], ignore_index=True)

        

    rna_seq_data_filtered['Hugo_Symbol'] = rna_seq_data_filtered['Hugo_Symbol'].str.strip()

        # Transform from wide to long format
    rna_seq_data_long = rna_seq_data_filtered.melt(
        id_vars=['Hugo_Symbol', 'Entrez_Gene_Id'],  # Keep these columns as identifiers
        value_vars=rna_seq_data_filtered.columns[2:],  # Patient columns to melt
        var_name='patient_id',  # Name for the new column containing patient IDs
        value_name='expression_value'  # Name for the new column containing expression values
    )

    rna_seq_data_final = rna_seq_data_long[['patient_id', 'Hugo_Symbol', 'expression_value']]
    rna_seq_data_final = rna_seq_data_final.rename(columns={'Hugo_Symbol': 'gene_symbol', 'patient_id':'model_id', 'expression_value':'rsem_tpm'})
    
     # Strip whitespace from gene_symbol column
    rna_seq_data_final['gene_symbol'] = rna_seq_data_final['gene_symbol'].str.strip()


    return rna_seq_data_final





def preprocess_cnv(patients_ids, cnv_data, all_montagud_nodes, synonyms_to_nodes_dict):
    
    cnv_data_filtered = cnv_data.melt(
    id_vars=['Hugo_Symbol', 'Entrez_Gene_Id'],
    value_vars=cnv_data.columns[2:],
    var_name='model_id',
    value_name='total_copy_number'
)
    
    cnv_data_filtered = cnv_data_filtered.rename(columns={'Hugo_Symbol': 'gene_symbol', 'patient_id':'model_id'})


    cnv_data_filtered = cnv_data_filtered[
        cnv_data_filtered["model_id"].isin(patients_ids)
    ]


    conditions = [
        cnv_data_filtered['total_copy_number'] == -2,
        cnv_data_filtered['total_copy_number'] == -1,  
        cnv_data_filtered['total_copy_number'] == 0,
        cnv_data_filtered['total_copy_number'] == 1,
        cnv_data_filtered['total_copy_number'] == 2
    ]

    choices = ['deletion', 'loss', 'normal', 'gain', 'amplification']
    cnv_data_filtered['effect'] = np.select(conditions, choices, default='other')

    # filter out Neural Category
    cnv_data_filtered = cnv_data_filtered[~cnv_data_filtered["effect"].isin(["normal", 'loss', 'gain'])]




    # keep only the patients id and montagud nodes
    cnv_data_filtered = cnv_data_filtered[
        cnv_data_filtered["gene_symbol"].isin(all_montagud_nodes)
    ]
    
    cnv_data_filtered = cnv_data_filtered[
        ["model_id", "gene_symbol", "total_copy_number", "effect"]
    ]

    cnv_data_filtered = cnv_data_filtered[
        ~cnv_data_filtered["total_copy_number"].isna()
    ]

    # remplace the cnv gene symbol column names by its synonyms in the proteins_synonyms_maps dictionary

    # apply the synonym mapping 
    for gene in cnv_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict:
            cnv_data_filtered.loc[cnv_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]

    return cnv_data_filtered





def pre_process_mutations(patients_ids, mutations_data, onco_tsg_data, all_montagud_nodes, synonyms_to_nodes_dict):
    mutations_data_filtered = mutations_data[mutations_data['Tumor_Sample_Barcode'].isin(patients_ids)]
    inactivating_mutations = ['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Translation_Start_Site', 'Splice_Site']


    mutations_data_filtered = mutations_data_filtered[['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode']]
    mutations_data_filtered['inactivating_status'] = mutations_data_filtered['Variant_Classification'].isin(inactivating_mutations).astype(int)


    onco_tsg_data_filtered = onco_tsg_data[['Hugo Symbol', 'Is Oncogene', 'Is Tumor Suppressor Gene']]


    onco_tsg_data_filtered = onco_tsg_data_filtered[onco_tsg_data_filtered['Hugo Symbol'].isin(all_montagud_nodes)]


    oncogenes_data = onco_tsg_data_filtered[(onco_tsg_data_filtered['Is Oncogene'] == 'Yes') & (onco_tsg_data_filtered['Is Tumor Suppressor Gene'] == 'No')]
    oncogenes_list = list(set(oncogenes_data['Hugo Symbol']))
    
    tsg_data = onco_tsg_data_filtered[(onco_tsg_data_filtered['Is Oncogene'] == 'No') & (onco_tsg_data_filtered['Is Tumor Suppressor Gene'] == 'Yes')]
    tsg_list = list(set(tsg_data['Hugo Symbol']))




    mutations_data_filtered['oncogene_tsg'] = np.where(
        mutations_data_filtered['Hugo_Symbol'].isin(oncogenes_list), 'onco',
        np.where(mutations_data_filtered['Hugo_Symbol'].isin(tsg_list), 'tsg', 'unknown')
    )


    mutations_data_filtered['mutations_effect'] = np.where(
        mutations_data_filtered['inactivating_status'] == 1, 'KO',
        np.where(
            (mutations_data_filtered['inactivating_status'] == 0) & 
            (mutations_data_filtered['oncogene_tsg'] == 'onco'), 'KI', 
            'unknown'
        )
    )

    mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['mutations_effect']!='unknown']

    mutations_data_filtered = mutations_data_filtered[['Hugo_Symbol', 'Tumor_Sample_Barcode', 'mutations_effect']]


    # Replace KI with amplification and KO with loss
    mutations_data_filtered['effect'] = mutations_data_filtered['mutations_effect'].replace({
        'KI': 'amplification',
        'KO': 'deletion'
    })

    mutations_data_filtered = mutations_data_filtered.rename(columns={'Hugo_Symbol': 'gene_symbol', 'Tumor_Sample_Barcode': 'model_id'})


    # apply the synonym mapping 
    for gene in mutations_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict:
            mutations_data_filtered.loc[mutations_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]

    return mutations_data_filtered
        






# ========= Create personalized boolean networks ========= 


def create_generic_patients_cfgs_bnds(
    folder_generic_models,
    folder_models,
    patients_ids,
):
    # --- Templates ---
    cfg_template_path = folder_generic_models + "/Fumia2013.cfg"
    bnd_template_path = folder_generic_models + "/Fumia2013.bnd"



    # Read templates once
    with open(cfg_template_path, "r") as file:
        cfg_template_content = file.read()
    with open(bnd_template_path, "r") as file:
        bnd_template_content = file.read()


    for patient_id in patients_ids:
        # patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
        # patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

        patient_cfg = cfg_template_content
        patient_bnd = bnd_template_content


        cfg_output_path = os.path.join(
            folder_models, f"{patient_id}.cfg"
        )
        bnd_output_path = os.path.join(
            folder_models, f"{patient_id}.bnd"
        )

        with open(cfg_output_path, "w") as file:
            file.write(patient_cfg)
        with open(bnd_output_path, "w") as file:
            file.write(patient_bnd)
    logger.debug(
        "All .cfg and .bnd files created for sensitive, resistant and healthy patients."
    )





def generic_models_update_phenotypes(phenotype_interest, folder_models):
    for filename in os.listdir(folder_models):
        file_path = os.path.join(folder_models, filename)
        # Only process files ending with .cfg
        if os.path.isfile(file_path) and filename.endswith(".cfg"):
            with open(file_path, "r") as file:
                content = file.read()
            modified = False
            for node in phenotype_interest:
                # Regex pattern for is_internal assignment
                pattern = re.compile(rf"({re.escape(node)}\.is_internal\s*=\s*)(TRUE|FALSE)\s*;", re.IGNORECASE)
                
                if re.search(pattern, content):
                    content, n_subs = pattern.subn(rf"\g<1>FALSE;", content)
                    if n_subs > 0:
                        modified = True
                else:
                    # Add new assignment at the end if not present
                    content += f"\n{node}.is_internal=FALSE;"
                    modified = True

            
            # Save only if modified
            if modified:
                modified_file_path = os.path.join(folder_models, filename)
                with open(modified_file_path, "w") as file:
                    file.write(content)






#     return rna_normalized

def normalize_rna_seq_efficient(rna_seq_data_models_filtered, normalize_technique):
    rna_normalized = rna_seq_data_models_filtered.copy()

    def local_sigmoid_normalize_group(group):
        """
        Apply sigmoid normalization to a group (gene) -> when data is unimodal
        1. Subtract median from all values
        2. Calculate lambda using MAD (Median Absolute Deviation)
        3. Apply sigmoid function: 1 / (1 + exp(-λ * x))
        """
        # Step 1: Calculate median and subtract from all values
        median_value = group['rsem_tpm'].median()
        centered_values = group['rsem_tpm'] - median_value
        
        # Step 2: Calculate MAD (Median Absolute Deviation)
        # MAD = median(|xi - median(X)|)
        mad_value = np.median(np.abs(centered_values))
        
        # Handle edge case where MAD is 0 (all values are identical)
        if mad_value == 0:
            # If all values are the same, set them all to 0.5 (median mapping)
            sigmoid_values = np.full(len(group), 0.5)
        else:
            # Step 3: Calculate lambda (λ) = log(3) / MAD
            lambda_param = np.log(3) / mad_value
            
            # Step 4: Apply sigmoid function with lambda scaling
            # sigmoid(x) = 1 / (1 + exp(-λ * x))
            sigmoid_values = 1 / (1 + np.exp(-lambda_param * centered_values))
        
        # Store the normalized values
        group['rsem_tpm_normalized'] = sigmoid_values
        
        return group
    
    def local_min_max_normalize_group(group):
        """
        Apply min-max normalization to a group (gene)
        Formula: (x - min) / (max - min)
        Maps values to [0, 1] range
        """

        min_value = group['rsem_tpm'].min()
        max_value = group['rsem_tpm'].max()
        value_range = max_value - min_value

        if value_range == 0:
        # If all values are the same, set them all to 0.5
            group['rsem_tpm_normalized'] = 0.5
        else:
            # Apply min-max normalization: (x - min) / (max - min)
            group['rsem_tpm_normalized'] = (group['rsem_tpm'] - min_value) / value_range
        
        return group
    
    def local_log_transform_group(group):
        """
        Apply log transformation followed by min-max normalization to a group (gene)
        Steps:
        1. Apply log2 transformation with pseudocount: log2(x + 1)
        2. Apply min-max normalization to map to [0, 1] range
        """

        # Step 1: Log2 transformation with pseudocount
        log_values = np.log2(group['rsem_tpm'] + 1)
        
        # Step 2: Min-max normalization of log values
        min_log = log_values.min()
        max_log = log_values.max()
        log_range = max_log - min_log
        
        # Handle edge case where all values are identical after log transformation
        if log_range == 0:
            # If all values are the same, set them all to 0.5
            group['rsem_tpm_normalized'] = 0.5
        else:
            # Apply min-max normalization: (log_x - log_min) / (log_max - log_min)
            group['rsem_tpm_normalized'] = (log_values - min_log) / log_range

        return group
    

    def global_minmax_normalize():
        """Global min-max normalization"""
        global_min = rna_normalized['rsem_tpm'].min()
        global_max = rna_normalized['rsem_tpm'].max()
        global_range = global_max - global_min
        
        if global_range > 0:
            rna_normalized['rsem_tpm_normalized'] = (rna_normalized['rsem_tpm'] - global_min) / global_range
        else:
            rna_normalized['rsem_tpm_normalized'] = 0.5
        
        logger.debug(f"Global normalization: min={global_min:.2f}, max={global_max:.2f}")
        return rna_normalized
    

    def global_log_normalize():
        """Log + global normalization"""
        log_values = np.log2(rna_normalized['rsem_tpm'] + 1)
        log_min = log_values.min()
        log_max = log_values.max()
        log_range = log_max - log_min
        
        if log_range > 0:
            rna_normalized['rsem_tpm_normalized'] = (log_values - log_min) / log_range
        else:
            rna_normalized['rsem_tpm_normalized'] = 0.5
        
        logger.debug(f"Log+Global normalization: log_min={log_min:.2f}, log_max={log_max:.2f}")
        return rna_normalized
        
   # Group by gene and apply normalization
    if normalize_technique == 'sigmoid':
        logger.debug("Applying global sigmoid normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_sigmoid_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'min-max':
        logger.debug("Applying global min-max normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_min_max_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'log_transf':
        logger.debug("Applying global log_transf normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_log_transform_group).reset_index(drop=True)


    elif normalize_technique == 'global_minmax':
        logger.debug("Applying global min-max normalization...")
        rna_normalized = global_minmax_normalize()
    
    elif normalize_technique == 'global_log':
        logger.debug("Applying log + global normalization...")
        rna_normalized = global_log_normalize()
    


    elif normalize_technique == 'distribution_normalization':
        logger.debug("Applying distribution normalization (paper)...")

        rna_normalized_filt = rna_seq_data_models_filtered[['model_id', 'gene_symbol', 'rsem_tpm']]

        rna_normalized_filt = rna_normalized_filt.rename(columns={'gene_symbol': 'Hugo_Symbol'})

        rna_normalized_filt = rna_normalized_filt.groupby(['model_id', 'Hugo_Symbol'], as_index=False).agg({'rsem_tpm': 'mean'})

        rna_normalized_filt = rna_normalized_filt.pivot(
        index='Hugo_Symbol',      # Genes as index
        columns='model_id',       # Patients as columns
        values='rsem_tpm'         # Expression values
        )

        rna_normalized_filt.index.name = 'Hugo_Symbol'
        rna_normalized_filt.columns.name = None


        rna_seq_data_filtered_analysis_pivoted_distribution = identify_genes_distribution(rna_normalized_filt)
        rna_normalized = compute_multi_distrib_normalization(rna_seq_data_filtered_analysis_pivoted_distribution)


    return rna_normalized



def personalized_patients_genes_cfgs(
    rna_seq_data_models_filtered,
    all_montagud_nodes,
    folder_models,
    amplif_factor,
    normalization_method
):
    # Apply the normalization
   
#     rna_seq_data_models_filtered = rna_seq_data_models_filtered.loc[
#     rna_seq_data_models_filtered.index.isin(all_montagud_nodes)
# ]
    
#     print('rna_seq_data_models_filtered:', rna_seq_data_models_filtered)

    rna_seq_data_normalized = normalize_rna_seq_efficient(rna_seq_data_models_filtered, normalization_method)


    # Loop through each file in the directory
    for filename in os.listdir(folder_models):
        if not filename.endswith('.cfg'):
            continue

        file_path = os.path.join(folder_models, filename)
        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0]
            
            with open(file_path, "r") as file:
                content = file.read()
      

            for gene in all_montagud_nodes:
                gene_data_normalized = rna_seq_data_normalized[
                        (rna_seq_data_normalized["gene_symbol"] == gene) &
                        (rna_seq_data_normalized["model_id"] == model_id_in_file)
                    ]['rsem_tpm_normalized']
                
                if gene_data_normalized.empty:
                        logger.debug(f"  No data for {gene}")
                        continue

                # Calculate the transition up and down values

                k_up = amplif_factor**(2 * (gene_data_normalized.iloc[0] - 0.5))
                if k_up != 0:
                    k_down = 1/ k_up
                else:
                    # Handle edge case
                    k_up = 0.001
                    k_down = 1000

                # Apply modifications
                u_pattern = re.compile(
                    rf"\$u_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?\s*;", 
                    re.DOTALL
                )
                d_pattern = re.compile(
                    rf"\$d_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?\s*;", 
                    re.DOTALL
                )

                u_line = f"$u_{gene} = {k_up:.4f};"
                d_line = f"$d_{gene} = {k_down:.4f};"
                
                content = re.sub(u_pattern, u_line, content)
                content = re.sub(d_pattern, d_line, content)



        # Save modified file
        modified_file_path = os.path.join(folder_models, filename)
        with open(modified_file_path, "w") as file:
            file.write(content)


def tailor_bnd_cnv_cm(cnv_data_filtered, folder_models):
    gain_group = cnv_data_filtered[
        cnv_data_filtered["effect"]=='amplification'
    ]

    loss_group = cnv_data_filtered[
        cnv_data_filtered["effect"]=='deletion'
    ]

    gain_ids = set(gain_group["model_id"])
    loss_ids = set(loss_group["model_id"])
    patients_all = list(gain_ids | loss_ids)

    for patient in patients_all:
        files_categ = [
            file
            for file in os.listdir(folder_models)
            if file.endswith(".bnd") and patient in file
        ]
        if not files_categ:
            logger.debug(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(folder_models, files_categ[0])

        if bnd_file is None:
            logger.debug(f"No .bnd file found for patient: {patient}")
            continue

        genes = cnv_data_filtered[cnv_data_filtered["model_id"] == patient][
            "gene_symbol"
        ]
        gene_list = genes.tolist()

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0]
        with open(bnd_file, "r") as file:
            content = file.read()

        modified_any = False
        for gene in gene_list:            
            # Check if this specific patient-gene combination is in both gain and loss
            patient_gene_in_gain = len(gain_group[
                (gain_group["model_id"] == patient) & 
                (gain_group["gene_symbol"] == gene)
            ]) > 0
            
            patient_gene_in_loss = len(loss_group[
                (loss_group["model_id"] == patient) & 
                (loss_group["gene_symbol"] == gene)
            ]) > 0
            
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            if patient_gene_in_gain and patient_gene_in_loss:
                logger.warning(
                    f"Patient {patient} with gene {gene} is in both gain and loss groups. Please review."
                )
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient_gene_in_gain:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient_gene_in_loss:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                logger.debug(f"Patient {patient} with gene {gene} not found in gain or loss CNV groups.")
                continue

            gene_match = pattern.search(content)
            if gene_match:
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    modified_any = True
            else:
                logger.debug(f"No {gene} node found in file for patient {patient_id}")

        if modified_any:
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            logger.debug(f"{patient_id}: CNV — no nodes modified")



def compute_phenotype_table(
    folder_save_results,
    folder_models,
    patient_id,
    input_nodes,
    phenotypes_interest,
):


    model_pers_bnd = f"{folder_models}/{patient_id}.bnd"
    model_pers_cfg = f"{folder_models}/{patient_id}.cfg"

    if not os.path.exists(model_pers_bnd) or not os.path.exists(model_pers_cfg):
        logger.warning(f"Missing model files for {patient_id}")
        return None

    model_pers_lung = maboss.load(model_pers_bnd, model_pers_cfg)

    results = pd.DataFrame(index=input_nodes, columns=phenotypes_interest)

    threshold_diff = 0.01
    initial_max_time = 5
    max_time_limit = 30

    for active_node in input_nodes:
        model_pers_lung.network.set_istate(active_node, [0, 1])
        for inactive_node in input_nodes:
            if inactive_node != active_node:
                model_pers_lung.network.set_istate(inactive_node, [1, 0])

        converged = False
        current_max_time = initial_max_time
        while not converged and current_max_time < max_time_limit:
            model_pers_lung.update_parameters(
                time_tick=0.1, max_time=current_max_time, sample_count=500
            )
            res_lung_pers = model_pers_lung.run()

            last_state = res_lung_pers.get_nodes_probtraj()
            if len(last_state) < 10:
                logger.warning(
                    f"[{active_node}] Warning: Not enough data for convergence check."
                )
                break
            diff = last_state.diff().abs()
            max_change_recent = diff.tail(10).max()
            if (max_change_recent > threshold_diff).any():
                current_max_time += 1
            else:
                converged = True
        if not converged:
            logger.warning(
                f"[{active_node}] Did not fully converge by max_time = {max_time_limit}."
            )
        final_probs = last_state.iloc[-1]
        for phenotype in phenotypes_interest:
            if phenotype in final_probs.index:
                results.loc[active_node, phenotype] = final_probs[phenotype]
    
    os.makedirs(folder_save_results, exist_ok=True)  # Ensure the folder exists
    results.to_csv(f"{folder_save_results}/{patient_id}.csv")
    return results




def combine_patient_results(results_folder):
    """
    Combine individual patient CSV files into one dataframe
    """
    all_data = []
    
    # Get all CSV files in the results folder
    csv_files = [f for f in os.listdir(results_folder) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        # Extract patient ID from filename (remove .csv extension)
        patient_id = Path(csv_file).stem
        
        # Read the CSV file
        df = pd.read_csv(os.path.join(results_folder, csv_file), index_col=0)
        
        # Calculate mean values across all inputs for each phenotype
        patient_data = {
            'patient_id': patient_id,
            'Proliferation': df['Proliferation'].mean(),
            'Apoptosis': df['Apoptosis'].mean()
        }
        
        all_data.append(patient_data)
    
    # Create combined dataframe
    combined_df = pd.DataFrame(all_data)
    combined_df.set_index('patient_id', inplace=True)
    combined_df.to_csv(f'{results_folder}/combined_patients_ids.csv')
    
    return combined_df





def correlate_boolean_predictions_with_gene_signatures(validation, proba_phenotype, hallmark, phenotype, rna_seq_data, patients_ids):

    # Melt the dataframe so that patient columns become rows with a 'model_id' column
    if validation:
        rna_seq_breast_long = pd.melt(
            rna_seq_data,
            id_vars=['Hugo_Symbol', 'Entrez_Gene_Id'],  # Keep these columns as identifiers
            value_vars=rna_seq_data.columns[2:],      # All patient/model columns
            var_name='model_id',                        # New column for patient/model ID
            value_name='expression_value'               # New column for expression value
        )

        rna_seq_breast_long = rna_seq_breast_long.rename(columns={'Hugo_Symbol': 'gene_symbol', 'expression_value':'rsem_tpm'})
    else:
        rna_seq_breast_long = rna_seq_data.copy()
        
    # Preprocessing -> get gene expression matrix
    rna_seq_data_filtered = rna_seq_breast_long[rna_seq_breast_long['model_id'].isin(patients_ids)]


    rna_seq_data_filtered = rna_seq_data_filtered[['model_id', 'gene_symbol', 'rsem_tpm']]


    # remove duplicates (compute mean)
    rna_seq_data_deduplicated = rna_seq_data_filtered.groupby(['model_id', 'gene_symbol']).agg({
        'rsem_tpm': 'mean'
    }).reset_index()


    # Reshape the DataFrame: genes as rows (index), samples as columns, expression as values
    gene_expression_matrix = rna_seq_data_deduplicated.pivot(
        index='gene_symbol',     # Genes become row index
        columns='model_id',      # Model IDs become columns
        values='rsem_tpm'        # Expression values fill the matrix
    )



    # extract genes signature for the specific hallmark given
    msigdb_hallmark = gp.get_library('MSigDB_Hallmark_2020')
    hallmark_phenotype = msigdb_hallmark[hallmark]


    # get phenotype values 

    phenotype_score_df = []
    for patient_id in patients_ids:
       
        sample_expr = gene_expression_matrix[patient_id].dropna()
        # Run ssGSEA for this sample
        gene_sets = {hallmark: hallmark_phenotype}

        ss_result = gp.ssgsea(
            data=sample_expr,
            gene_sets=gene_sets,
            outdir=None,  # Don't save files
            no_plot=True,
            processes=1
        )

        phenotype_score = ss_result.res2d.loc[0, 'ES']

        phenotype_score_df.append({
            'patient_id': patient_id,
            f'{hallmark}_gene_signature_score': phenotype_score
        })
    
     

    # filtering proba phenotype df

    phenotype_cols = [col for col in proba_phenotype.columns if phenotype in col]  
    proba_phenotype_filtered = proba_phenotype[phenotype_cols]
    
    proba_phenotype_filtered = proba_phenotype_filtered[proba_phenotype_filtered.index.isin(patients_ids)]
    
    proba_phenotype_filtered = proba_phenotype_filtered.reset_index()
    proba_phenotype_filtered = proba_phenotype_filtered.rename(columns={proba_phenotype_filtered.columns[0]: 'patient_id'})

    phenotype_scores_df = pd.DataFrame(phenotype_score_df) 

    # create correlation data (combination of the proba phenotype and score calculated before)
    correlation_data = phenotype_scores_df.merge(proba_phenotype_filtered, on='patient_id', how='inner')
    results_corr_df = pd.DataFrame(index = phenotype_cols, columns = ['correlation_value', 'p_val'])


    for col in phenotype_cols:
        corr, pval = spearmanr(correlation_data[col], correlation_data[f'{hallmark}_gene_signature_score'])
        results_corr_df.loc[col, 'correlation_value'] = corr
        results_corr_df.loc[col, 'p_val'] = pval


    # Add mean row
    mean_corr = results_corr_df['correlation_value'].astype(float).mean()
    mean_pval = results_corr_df['p_val'].astype(float).mean()

    results_corr_df.loc['Mean', 'correlation_value'] = mean_corr
    results_corr_df.loc['Mean', 'p_val'] = mean_pval


    return results_corr_df, correlation_data

# vizualise the results 
def visualize_correlation_results(results_corr_df, correlation_data, hallmark):
    """
    Create comprehensive visualization of Boolean network validation results
    """
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Get gene signature column name
    gene_signature_col = f'{hallmark}_gene_signature_score'
    gene_scores = correlation_data[gene_signature_col]
    
    # Get Boolean network columns (excluding patient_id and gene signature)
    boolean_cols = [col for col in correlation_data.columns 
                   if col not in ['patient_id', gene_signature_col]]
    
    # Calculate number of subplots
    n_scenarios = len(boolean_cols)
    n_cols = min(3, n_scenarios)
    n_rows = max(2, (n_scenarios + n_cols - 1) // n_cols)
    
    fig = plt.figure(figsize=(5*n_cols, 4*n_rows))
    
    # Main title
    fig.suptitle(f'Boolean Network Validation: {hallmark}', fontsize=16, fontweight='bold')
    
    # Individual scatter plots
    for i, col in enumerate(boolean_cols):
        ax = plt.subplot(n_rows, n_cols, i+1)
        
        # Get data for this scenario
        boolean_probs = correlation_data[col]
        
        # Get correlation from results
        corr = float(results_corr_df.loc[col, 'correlation_value'])
        pval = float(results_corr_df.loc[col, 'p_val'])
        
        # Create scatter plot with patient labels
        scatter = ax.scatter(gene_scores, boolean_probs, 
                           alpha=0.8, s=120, 
                           c=range(len(gene_scores)), 
                           cmap='viridis',
                           edgecolors='black', linewidth=0.5)
        
        # Add patient labels
        for j, patient_id in enumerate(correlation_data['patient_id']):
            ax.annotate(patient_id, 
                       (gene_scores.iloc[j], boolean_probs.iloc[j]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=9, alpha=0.8, fontweight='bold')
        
        # Add correlation line
        if len(gene_scores) > 1:  # Only if we have multiple points
            z = np.polyfit(gene_scores, boolean_probs, 1)
            p = np.poly1d(z)
            x_line = np.linspace(gene_scores.min(), gene_scores.max(), 100)
            ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2)
        
        # Formatting
        ax.set_xlabel(f'{hallmark}\nGene Signature Score', fontsize=10)
        ax.set_ylabel(f'Boolean Network\n{col.replace("_", " ")}', fontsize=10)
        
        # Title with correlation info
        if pval < 0.001:
            sig_text = "***"
        elif pval < 0.01:
            sig_text = "**"
        elif pval < 0.05:
            sig_text = "*"
        else:
            sig_text = "ns"
        
        title_text = f'{col.replace("_", " ")}\nr = {corr:.3f} ({sig_text})'
        ax.set_title(title_text, fontsize=11, fontweight='bold')
        
        # Color code the title based on significance
        if pval < 0.05:
            ax.title.set_color('green')
        else:
            ax.title.set_color('red')
        
        ax.grid(True, alpha=0.3)
        
        # Add correlation strength text
        abs_corr = abs(corr)
        if abs_corr >= 0.7:
            strength = "Strong"
            strength_color = 'darkgreen'
        elif abs_corr >= 0.5:
            strength = "Moderate"
            strength_color = 'orange'
        elif abs_corr >= 0.3:
            strength = "Weak"
            strength_color = 'red'
        else:
            strength = "Very Weak"
            strength_color = 'darkred'
        
        ax.text(0.05, 0.95, f'{strength}', transform=ax.transAxes,
                fontsize=9, fontweight='bold', color=strength_color,
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
    
    # Summary statistics subplot (if space allows)
    if n_scenarios < n_rows * n_cols:
        ax_summary = plt.subplot(n_rows, n_cols, n_scenarios + 1)
        
        # Create summary statistics
        correlations = results_corr_df['correlation_value'].astype(float)
        p_values = results_corr_df['p_val'].astype(float)
        
        # Bar plot of correlations
        scenario_names = [name.split('_')[0] for name in results_corr_df.index]
        bars = ax_summary.bar(range(len(correlations)), correlations, 
                             color=['green' if p < 0.05 else 'red' for p in p_values],
                             alpha=0.7, edgecolor='black')
        
        ax_summary.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax_summary.axhline(y=0.7, color='green', linestyle='--', alpha=0.5, label='Strong (0.7)')
        ax_summary.axhline(y=-0.7, color='green', linestyle='--', alpha=0.5)
        ax_summary.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='Moderate (0.5)')
        ax_summary.axhline(y=-0.5, color='orange', linestyle='--', alpha=0.5)
        
        ax_summary.set_xlabel('Scenarios')
        ax_summary.set_ylabel('Correlation Coefficient')
        ax_summary.set_title('Correlation Summary', fontweight='bold')
        ax_summary.set_xticks(range(len(correlations)))
        ax_summary.set_xticklabels(scenario_names, rotation=45, ha='right')
        ax_summary.legend()
        ax_summary.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (bar, corr, pval) in enumerate(zip(bars, correlations, p_values)):
            height = bar.get_height()
            ax_summary.text(bar.get_x() + bar.get_width()/2., height + 0.02*np.sign(height),
                           f'{corr:.2f}', ha='center', va='bottom' if height > 0 else 'top',
                           fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    return fig






# Survival Analysis 
def survival_analysis_comparison(normalization_method, discrete_variable, continuous_variable, data, group_col, time_col='OS_MONTHS', event_col='OS_STATUS', save_plots=True):
    """
    Perform complete survival analysis comparison between groups
    """
    
    # Get unique groups
    groups = data[group_col].unique()
    
    if len(groups) != 2:
        logger.warning(f"Warning: Expected 2 groups, found {len(groups)}")
        return None
    
    # Separate groups
    group1_data = data[data[group_col] == groups[0]]
    group2_data = data[data[group_col] == groups[1]]
    
    
    # Log-rank test
    results = logrank_test(
        group1_data[time_col], 
        group2_data[time_col],
        group1_data[event_col], 
        group2_data[event_col]
    )
    
    
    # Plot Kaplan-Meier curves
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(10, 6))
    
    # Plot each group
    colors = ['red', 'blue']
    for i, group in enumerate(groups):
        group_data = data[data[group_col] == group]
        kmf.fit(group_data[time_col], group_data[event_col], label=f'{group_col}: {group}')
        kmf.plot_survival_function(color=colors[i])
    
    plt.title(f'Survival Analysis by {group_col}\nLog-rank p-value: {results.p_value:.4f}')
    plt.xlabel('Time (months)')
    plt.ylabel('Survival Probability')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save plot if requested
    if save_plots:
        # Create output directory
    
        output_dir = f'analysis/validation_Breast/{discrete_variable}_{continuous_variable}/{normalization_method}/results/outputs'
        os.makedirs(output_dir, exist_ok=True)
        
        # Create filename
        filename = f"survival_analysis_{group_col.lower().replace('_', '-')}.png"
        filepath = os.path.join(output_dir, filename)
        
        # Save the plot
        plt.savefig(filepath, dpi=300, bbox_inches='tight', facecolor='white')    
    # plt.show()
    


    