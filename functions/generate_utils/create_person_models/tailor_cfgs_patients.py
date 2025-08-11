# pipeline to tailor the generic lung adenocarcinoma boolean network to patient specific boolean

import pandas as pd
import re
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler


from functions.generate_utils.pre_process_data.pre_process_genes import (
    process_genes, identify_genes_distribution, compute_multi_distrib_normalization,
)
import logging

logger = logging.getLogger(__name__)

def normalize_rna_seq_efficient(data_filtered, normalize_technique):
    data_normalized = data_filtered.copy()
    symbol_col = 'gene_symbol' if 'gene_symbol' in data_filtered.columns else 'protein_symbol'

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
        global_min = data_normalized['rsem_tpm'].min()
        global_max = data_normalized['rsem_tpm'].max()
        global_range = global_max - global_min
        
        if global_range > 0:
            data_normalized['rsem_tpm_normalized'] = (data_normalized['rsem_tpm'] - global_min) / global_range
        else:
            data_normalized['rsem_tpm_normalized'] = 0.5
        
        logger.debug(f"Global normalization: min={global_min:.2f}, max={global_max:.2f}")
        return data_normalized


    def global_log_normalize():
        """Log + global normalization"""
        log_values = np.log2(data_normalized['rsem_tpm'] + 1)
        log_min = log_values.min()
        log_max = log_values.max()
        log_range = log_max - log_min
        
        if log_range > 0:
            data_normalized['rsem_tpm_normalized'] = (log_values - log_min) / log_range
        else:
            data_normalized['rsem_tpm_normalized'] = 0.5
        
        logger.debug(f"Log+Global normalization: log_min={log_min:.2f}, log_max={log_max:.2f}")
        return data_normalized
        

   # Group by gene and apply normalization
    if normalize_technique == 'sigmoid':
        logger.debug("Applying global sigmoid normalization...")
        data_normalized = data_normalized.groupby(symbol_col).apply(local_sigmoid_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'min-max':
        logger.debug("Applying global min-max normalization...")
        data_normalized = data_normalized.groupby(symbol_col).apply(local_min_max_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'log_transf':
        logger.debug("Applying global log_transf normalization...")
        data_normalized = data_normalized.groupby(symbol_col).apply(local_log_transform_group).reset_index(drop=True)


    elif normalize_technique == 'global_minmax':
        logger.debug("Applying global min-max normalization...")
        data_normalized = global_minmax_normalize()
    
    elif normalize_technique == 'global_log':
        logger.debug("Applying log + global normalization...")
        data_normalized = global_log_normalize()



# test
    elif normalize_technique == 'distribution_normalization':
        logger.debug("Applying distribution normalization (paper)...")



        data_normalized_filt = data_filtered[['model_id', symbol_col, 'rsem_tpm']]

        data_normalized_filt = data_normalized_filt.rename(columns={symbol_col: 'Hugo_Symbol'})

        data_normalized_filt = data_normalized_filt.groupby(['model_id', 'Hugo_Symbol'], as_index=False).agg({'rsem_tpm': 'mean'})

        data_normalized_filt = data_normalized_filt.pivot(
        index='Hugo_Symbol',      # Genes as index
        columns='model_id',       # Patients as columns
        values='rsem_tpm'         # Expression values
        )

        data_normalized_filt.index.name = 'Hugo_Symbol'
        data_normalized_filt.columns.name = None



        data_filtered_analysis_pivoted_distribution = identify_genes_distribution(data_normalized_filt)
        data_normalized = compute_multi_distrib_normalization(data_filtered_analysis_pivoted_distribution)

    return data_normalized



def personalized_patients_genes_cfgs(
    rna_seq_data_models_filtered,
    montagud_node_model,
    folder_models,
    amplif_factor,
    context_label,
    normalization_method,
):
    # Apply the normalization
    rna_seq_data_normalized = normalize_rna_seq_efficient(rna_seq_data_models_filtered, normalization_method)


    # Loop through each file in the directory
    for filename in os.listdir(folder_models):
        if not filename.endswith('.cfg'):
            continue

        file_path = os.path.join(folder_models, filename)

        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0].replace(
                f"_{context_label}", ""
            )
            with open(file_path, "r") as file:
                content = file.read()
    
            for gene in montagud_node_model:
                gene_data = rna_seq_data_normalized[
                    (rna_seq_data_normalized["gene_symbol"] == gene) &
                    (rna_seq_data_normalized["model_id"] == model_id_in_file)
                ]['rsem_tpm_normalized']
                
                if gene_data.empty:
                    logger.debug(f"  No data for {gene}")
                    continue
            
                # Get the expression value
                expression_value = gene_data.iloc[0]
                # Calculate the transition up and down values
                
                k_up = amplif_factor**(2 * (expression_value - 0.5))
               
                # Handle edge cases
                if k_up <= 0 or not np.isfinite(k_up):
                    k_up = 0.001
                    k_down = 1000
                else:
                    k_down = 1 / k_up

                # Apply modifications
                u_pattern = re.compile(rf"\$u_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                d_pattern = re.compile(rf"\$d_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                
                u_line = f"$u_{gene} = {k_up:.4f};"
                d_line = f"$d_{gene} = {k_down:.4f};"
                
                content = re.sub(u_pattern, u_line, content)
                content = re.sub(d_pattern, d_line, content)

        # Save modified file
        modified_file_path = os.path.join(folder_models, filename)
        with open(modified_file_path, "w") as file:
            file.write(content)






# TEST 

def personalized_patients_proteins_cfgs(
    protein_models_filtered,
    montagud_node_model,
    folder_models,
    amplif_factor,
    context_label,
    normalization_method,
):
    # Apply the normalization
    rna_seq_data_normalized = normalize_rna_seq_efficient(protein_models_filtered, normalization_method)


    # Loop through each file in the directory
    for filename in os.listdir(folder_models):
        if not filename.endswith('.cfg'):
            continue

        file_path = os.path.join(folder_models, filename)

        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0].replace(
                f"_{context_label}", ""
            )
            with open(file_path, "r") as file:
                content = file.read()
    
            for protein in montagud_node_model:
                protein_data = rna_seq_data_normalized[
                    (rna_seq_data_normalized["protein_symbol"] == protein) &
                    (rna_seq_data_normalized["model_id"] == model_id_in_file)
                ]['rsem_tpm_normalized']
                
                if protein_data.empty:
                    logger.debug(f"  No data for {protein}")
                    continue
            
                # Get the expression value
                expression_value = protein_data.iloc[0]
                # Calculate the transition up and down values
                
                k_up = amplif_factor**(2 * (expression_value - 0.5))
               
                # Handle edge cases
                if k_up <= 0 or not np.isfinite(k_up):
                    k_up = 0.001
                    k_down = 1000
                else:
                    k_down = 1 / k_up

                # Apply modifications
                u_pattern = re.compile(rf"\$u_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                d_pattern = re.compile(rf"\$d_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                
                u_line = f"$u_{protein} = {k_up:.4f};"
                d_line = f"$d_{protein} = {k_down:.4f};"
                
                content = re.sub(u_pattern, u_line, content)
                content = re.sub(d_pattern, d_line, content)

        # Save modified file
        modified_file_path = os.path.join(folder_models, filename)
        with open(modified_file_path, "w") as file:
            file.write(content)






# def personalized_patients_proteins_cfgs(
#     proteins_models_filtered,
#     montagud_node_model,
#     folder_models,
#     context_label,
#     normalization_method,
# ):


#     df_melted_proteins = df_melted_proteins[["model_id", "protein_symbol", "rsem_tpm"]]


#     protein_data_max = df_melted_proteins.copy()
#     protein_data_max["protein_max"] = protein_data_max.groupby("protein_symbol")[
#         "rsem_tpm"
#     ].transform("max")

#     # Loop through each file in the directory
#     for filename in os.listdir(folder_models):
#         if not filename.endswith(".cfg"):
#             continue
#         file_path = os.path.join(folder_models, filename)
#         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory

#             model_id_in_file = os.path.splitext(filename)[0].replace(
#                 f"_{context_label}", ""
#             )
            
#             # Only proceed if model_id is in table_proteins_patients
#             with open(file_path, "r") as file:
#                 content = file.read()
            
#             for protein in montagud_node_model:

#                 protein_max_data = protein_data_max[
#                     protein_data_max["protein_symbol"] == protein
#                 ]
                
#                 if protein_max_data.empty:
#                     logger.debug(f"  No max data for protein: {protein}")
#                     continue
                
#                 # Check if this patient has data for this protein
#                 patient_protein_data = df_melted_proteins[
#                     (df_melted_proteins["protein_symbol"] == protein) &
#                     (df_melted_proteins["model_id"] == model_id_in_file)
#                 ]
                
#                 if patient_protein_data.empty:
#                     logger.debug(f"  No data for protein {protein} in patient {model_id_in_file}")
#                     continue
                                

#                 expr_max = protein_max_data["protein_max"].iloc[0]
#                 prot_expr = patient_protein_data["rsem_tpm"].iloc[0]


#                 prob_1 = min(max(prot_expr / expr_max, 0), 1)

#                 u_pattern = re.compile(rf"\$u_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
#                 d_pattern = re.compile(rf"\$d_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)

#                 u_line = f"$u_{protein} = {prob_1:.4f};"
#                 d_line = f"$d_{protein} = {1 - prob_1:.4f};"
                    
                        
#                 content = re.sub(u_pattern, u_line, content)
#                 content = re.sub(d_pattern, d_line, content)


#                 # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
#         modified_file_path = os.path.join(folder_models, filename)
#         with open(modified_file_path, "w") as file:
#             file.write(content)

