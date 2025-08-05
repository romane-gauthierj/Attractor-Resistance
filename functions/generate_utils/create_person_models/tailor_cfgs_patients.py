# pipeline to tailor the generic lung adenocarcinoma boolean network to patient specific boolean

import pandas as pd
import re
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler

from functions.generate_utils.pre_process_data.pre_process_genes import (
    process_genes, identify_genes_distribution, compute_multi_distrib_normalization,
)

# Step 2 - personalization (with their method)

# from sklearn.preprocessing import MinMaxScaler

# def normalize_rna_seq_efficient(rna_seq_data_models_filtered):
#     """
#     Efficient min-max normalization using pandas groupby
#     """
#     rna_normalized = rna_seq_data_models_filtered.copy()
    
#     def minmax_normalize_group(group):
#         """Apply min-max normalization to a group (gene)"""
#         scaler = MinMaxScaler()
#         group['rsem_tpm_normalized'] = scaler.fit_transform(group[['rsem_tpm']]).flatten()
#         return group
    
#     print("Applying min-max normalization per gene...")
    
#     # Group by gene and apply normalization
#     rna_normalized = rna_normalized.groupby('gene_symbol').apply(minmax_normalize_group).reset_index(drop=True)
    
#     print("Normalization completed!")
#     print(f"Original column: 'rsem_tpm', Normalized column: 'rsem_tpm_normalized'")
#     return rna_normalized




# def normalize_rna_seq_efficient(rna_seq_data_models_filtered):
#     rna_normalized = rna_seq_data_models_filtered.copy()

#     def sigmoid_normalize_group(group):
#         """
#         Apply sigmoid normalization to a group (gene) -> when data is unimodal
#         1. Subtract median from all values
#         2. Calculate lambda using MAD (Median Absolute Deviation)
#         3. Apply sigmoid function: 1 / (1 + exp(-λ * x))
#         """
#         # Step 1: Calculate median and subtract from all values
#         median_value = group['rsem_tpm'].median()
#         centered_values = group['rsem_tpm'] - median_value
        
#         # Step 2: Calculate MAD (Median Absolute Deviation)
#         # MAD = median(|xi - median(X)|)
#         mad_value = np.median(np.abs(centered_values))
        
#         # Handle edge case where MAD is 0 (all values are identical)
#         if mad_value == 0:
#             # If all values are the same, set them all to 0.5 (median mapping)
#             sigmoid_values = np.full(len(group), 0.5)
#         else:
#             # Step 3: Calculate lambda (λ) = log(3) / MAD
#             lambda_param = np.log(3) / mad_value
            
#             # Step 4: Apply sigmoid function with lambda scaling
#             # sigmoid(x) = 1 / (1 + exp(-λ * x))
#             sigmoid_values = 1 / (1 + np.exp(-lambda_param * centered_values))
        
#         # Store the normalized values
#         group['rsem_tpm_normalized'] = sigmoid_values
        
#         # Store additional metrics for debugging/validation
#         group['median_value'] = median_value
#         group['mad_value'] = mad_value
#         group['lambda_param'] = lambda_param if mad_value != 0 else np.nan
        
#         return group
#    # Group by gene and apply normalization
#     rna_normalized = rna_normalized.groupby('gene_symbol').apply(sigmoid_normalize_group).reset_index(drop=True)

    
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
        
        print(f"Global normalization: min={global_min:.2f}, max={global_max:.2f}")
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
        
        print(f"Log+Global normalization: log_min={log_min:.2f}, log_max={log_max:.2f}")
        return rna_normalized
        

   # Group by gene and apply normalization
    if normalize_technique == 'sigmoid':
        print("Applying global sigmoid normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_sigmoid_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'min-max':
        print("Applying global min-max normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_min_max_normalize_group).reset_index(drop=True)

    elif normalize_technique == 'log_transf':
        print("Applying global log_transf normalization...")
        rna_normalized = rna_normalized.groupby('gene_symbol').apply(local_log_transform_group).reset_index(drop=True)


    elif normalize_technique == 'global_minmax':
        print("Applying global min-max normalization...")
        rna_normalized = global_minmax_normalize()
    
    elif normalize_technique == 'global_log':
        print("Applying log + global normalization...")
        rna_normalized = global_log_normalize()



# test
    elif normalize_technique == 'distribution_normalization':
        print("Applying distribution normalization (paper)...")

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



# checkkkk
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
                    print(f"  No data for {gene}")
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

                # if k_up != 0:
                #     k_down = 1/ k_up
                # else:
                #     # Handle edge case
                #     k_up = 0.001
                #     k_down = 1000

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




# def personalized_patients_genes_cfgs(
#     rna_seq_data_models_filtered,
#     montagud_node_model,
#     folder_models,
#     amplif_factor,
#     context_label,
# ):
#     # Apply the normalization
#     rna_seq_data_normalized = normalize_rna_seq_efficient(rna_seq_data_models_filtered)


#     print("==== Processing all genes per patient ===")
#     # Loop through each file in the directory
#     for filename in os.listdir(folder_models):
#         if not filename.endswith('.cfg'):
#             continue

#         file_path = os.path.join(folder_models, filename)
#         print(f"\nProcessing file: {filename}")

#         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
#             # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
#             # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
#             model_id_in_file = os.path.splitext(filename)[0].replace(
#                 f"_{context_label}", ""
#             )
            
#             with open(file_path, "r") as file:
#                 content = file.read()
      

#             for gene in montagud_node_model:
#                 gene_data_normalized = rna_seq_data_normalized[
#                         (rna_seq_data_normalized["gene_symbol"] == gene) &
#                         (rna_seq_data_normalized["model_id"] == model_id_in_file)
#                     ]['rsem_tpm_normalized']
                
#                 if gene_data_normalized.empty:
#                         print(f"  No data for {gene}")
#                         continue

#                 # Calculate the transition up and down values

#                 k_up = amplif_factor**(2 * (gene_data_normalized.iloc[0] - 0.5))
#                 if k_up != 0:
#                     k_down = 1/ k_up
#                 else:
#                     # Handle edge case
#                     k_up = 0.001
#                     k_down = 1000

#                 # Apply modifications
#                 u_pattern = re.compile(rf"\$u_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
#                 d_pattern = re.compile(rf"\$d_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                
#                 u_line = f"$u_{gene} = {k_up:.4f};"
#                 d_line = f"$d_{gene} = {k_down:.4f};"
                
#                 content = re.sub(u_pattern, u_line, content)
#                 content = re.sub(d_pattern, d_line, content)



#         # Save modified file
#         modified_file_path = os.path.join(folder_models, filename)
#         with open(modified_file_path, "w") as file:
#             file.write(content)
#         print(f"Modified and saved: {modified_file_path}")





# # Step 2 - personalization (with my method)
# # def personalized_patients_genes_cfgs(
# #     rna_seq_data_models_filtered,
# #     montagud_node_model,
# #     folder_models,
# #     context_label,
# # ):
# #     rna_seq_data_max = rna_seq_data_models_filtered.copy()
# #     rna_seq_data_max["gene_max"] = rna_seq_data_max.groupby("gene_symbol")[
# #         "rsem_tpm"
# #     ].transform("max")


# #     print("==== Processing all genes per patient ===")
# #     # Loop through each file in the directory
# #     for filename in os.listdir(folder_models):
# #         if not filename.endswith('.cfg'):
# #             continue

# #         file_path = os.path.join(folder_models, filename)
# #         print(f"\nProcessing file: {filename}")

# #         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
# #             # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
# #             # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
# #             model_id_in_file = os.path.splitext(filename)[0].replace(
# #                 f"_{context_label}", ""
# #             )
            
# #             with open(file_path, "r") as file:
# #                 content = file.read()
      

# #             for gene in montagud_node_model:
# #                 gene_data = rna_seq_data_models_filtered[
# #                         (rna_seq_data_models_filtered["gene_symbol"] == gene) &
# #                         (rna_seq_data_models_filtered["model_id"] == model_id_in_file)
# #                     ]
# #                 if gene_data.empty:
#                         print(f"  No data for {gene}")
#                         continue

#                 # Calculate expression probability
#                 expr_max = rna_seq_data_max.loc[
#                     rna_seq_data_max["gene_symbol"] == gene, "gene_max"
#                 ].iloc[0]
#                 gene_expr = gene_data["rsem_tpm"].iloc[0]
#                 prob_1 = min(max(gene_expr / expr_max, 0), 1)


#                 # Apply modifications
#                 u_pattern = re.compile(rf"\$u_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
#                 d_pattern = re.compile(rf"\$d_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                
#                 u_line = f"$u_{gene} = {prob_1:.4f};"
#                 d_line = f"$d_{gene} = {1 - prob_1:.4f};"
                
#                 content = re.sub(u_pattern, u_line, content)
#                 content = re.sub(d_pattern, d_line, content)



#         # Save modified file
#         modified_file_path = os.path.join(folder_models, filename)
#         with open(modified_file_path, "w") as file:
#             file.write(content)
#         print(f"Modified and saved: {modified_file_path}")



def personalized_patients_proteins_cfgs(
    df_melted_proteins,
    montagud_node_model,
    folder_models,
    context_label,
):

    df_melted_proteins = df_melted_proteins[["model_id", "protein_symbol", "rsem_tpm"]]


    protein_data_max = df_melted_proteins.copy()
    protein_data_max["protein_max"] = protein_data_max.groupby("protein_symbol")[
        "rsem_tpm"
    ].transform("max")

    # Loop through each file in the directory
    for filename in os.listdir(folder_models):
        if not filename.endswith(".cfg"):
            continue
        file_path = os.path.join(folder_models, filename)
        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory

            model_id_in_file = os.path.splitext(filename)[0].replace(
                f"_{context_label}", ""
            )
            
            # Only proceed if model_id is in table_proteins_patients
            with open(file_path, "r") as file:
                content = file.read()
            
            for protein in montagud_node_model:

                protein_max_data = protein_data_max[
                    protein_data_max["protein_symbol"] == protein
                ]
                
                if protein_max_data.empty:
                    print(f"  No max data for protein: {protein}")
                    continue
                
                # Check if this patient has data for this protein
                patient_protein_data = df_melted_proteins[
                    (df_melted_proteins["protein_symbol"] == protein) &
                    (df_melted_proteins["model_id"] == model_id_in_file)
                ]
                
                if patient_protein_data.empty:
                    print(f"  No data for protein {protein} in patient {model_id_in_file}")
                    continue
                                

                expr_max = protein_max_data["protein_max"].iloc[0]
                prot_expr = patient_protein_data["rsem_tpm"].iloc[0]


                prob_1 = min(max(prot_expr / expr_max, 0), 1)

                print(f"  Protein: {protein}, Expression: {prot_expr:.4f}, Max: {expr_max:.4f}, Prob: {prob_1:.4f}")


                u_pattern = re.compile(rf"\$u_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
                d_pattern = re.compile(rf"\$d_{re.escape(protein)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)

                u_line = f"$u_{protein} = {prob_1:.4f};"
                d_line = f"$d_{protein} = {1 - prob_1:.4f};"
                    
                        
                content = re.sub(u_pattern, u_line, content)
                content = re.sub(d_pattern, d_line, content)


                # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
        modified_file_path = os.path.join(folder_models, filename)
        with open(modified_file_path, "w") as file:
            file.write(content)
        print(f"Modified and saved: {modified_file_path}")



# def personalized_patients_genes_cfgs(
#     rna_seq_data_models_filtered,
#     montagud_node_model,
#     folder_models,
#     table_rna_seq_patients,
#     context_label,
# ):
   

#     rna_seq_data_max = rna_seq_data_models_filtered.copy()
#     rna_seq_data_max["gene_max"] = rna_seq_data_max.groupby("gene_symbol")[
#         "rsem_tpm"
#     ].transform("max")


#     print("==== Processing all genes per patient ===")
#     # Loop through each file in the directory
#     for filename in os.listdir(folder_models):
#         file_path = os.path.join(folder_models, filename)
#         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
#             # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
#             # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
#             model_id_in_file = os.path.splitext(filename)[0].replace(
#                 f"_{context_label}", ""
#             )
#             # Only proceed if model_id is in table_genes_patients
#             if model_id_in_file in table_rna_seq_patients.index:
#                 with open(file_path, "r") as file:
#                     content = file.read()
#                     # rna_seq_data_filtered_filtered = rna_seq_data_filtered[rna_seq_data_filtered.index == model_id_in_file]
#                 if model_id_in_file in table_rna_seq_patients.index:
#                     row = table_rna_seq_patients.loc[model_id_in_file]
#                     # Only modify if this row corresponds to the file being processed
#                     high_genes = row["High Gene Expression"]
#                     low_genes = row["Low Gene Expression"]
#                     gene_high_list = [
#                         gene.strip() for gene in high_genes.split(",") if gene.strip()
#                     ]
#                     gene_low_list = [
#                         gene.strip() for gene in low_genes.split(",") if gene.strip()
#                     ]

#                     gene_list = list(set(gene_high_list + gene_low_list))
#                     for gene in gene_list:  # add condition if gene is in the genes of the boolean generic network
#                         if gene in montagud_node_model:
#                             expr_max = rna_seq_data_max.loc[
#                                 rna_seq_data_max["gene_symbol"] == gene, "gene_max"
#                             ].iloc[0]
#                             gene_expr = rna_seq_data_models_filtered.loc[
#                                 (rna_seq_data_models_filtered["gene_symbol"] == gene)
#                                 & (rna_seq_data_models_filtered["model_id"] == model_id_in_file),
#                                 "rsem_tpm",
#                             ].iloc[0]
#                             prob_1 = min(max(gene_expr / expr_max, 0), 1)

                           
#                             u_pattern = re.compile(rf"\$u_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)
#                             d_pattern = re.compile(rf"\$d_{re.escape(gene)}\s*=\s*[0-9]*\.?[0-9]+;", re.DOTALL)


#                             u_line = f"$u_{gene} = {prob_1:.4f};"
#                             d_line = f"$d_{gene} = {1 - prob_1:.4f};"
#                             content = re.sub(u_pattern, u_line, content)
#                             content = re.sub(d_pattern, d_line, content)
#                 # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
#                 modified_file_path = os.path.join(folder_models, filename)

#                 with open(modified_file_path, "w") as file:
#                     file.write(content)
#                 print(f"Modified and saved: {modified_file_path}")


# function used for validation and general pipeline
# def personalized_patients_genes_cfgs(
#     rna_seq_data,
#     montagud_node_model,
#     folder_models,
#     patients_ids,
#     rna_seq_data_filtered,
#     context_label,
# ):
#     rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]

#     rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]

#     print('rna_seq_data', rna_seq_data)


#     rna_seq_data_max = rna_seq_data
#     rna_seq_data_max["gene_max"] = rna_seq_data_max.groupby("gene_symbol")[
#         "rsem_tpm"
#     ].transform("max")



#     # Loop through each file in the directory
#     for filename in os.listdir(folder_models):

#         file_path = os.path.join(folder_models, filename)
#         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
#             # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
#             # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
#             model_id_in_file = os.path.splitext(filename)[0].replace(
#                 f"_{context_label}", ""
#             )
#             # Only proceed if model_id is in table_genes_patients
#             if model_id_in_file in rna_seq_data_filtered.index:
#                 with open(file_path, "r") as file:
#                     content = file.read()
#                     # rna_seq_data_filtered_filtered = rna_seq_data_filtered[rna_seq_data_filtered.index == model_id_in_file]
#                 if model_id_in_file in rna_seq_data_filtered.index:
#                     row = rna_seq_data_filtered.loc[model_id_in_file]
#                     # Only modify if this row corresponds to the file being processed
#                     high_genes = row["High Gene Expression"]
#                     low_genes = row["Low Gene Expression"]
#                     gene_high_list = [
#                         gene.strip() for gene in high_genes.split(",") if gene.strip()
#                     ]
#                     gene_low_list = [
#                         gene.strip() for gene in low_genes.split(",") if gene.strip()
#                     ]
#                     gene_list = list(set(gene_high_list + gene_low_list))
#                     for gene in gene_list:  # add condition if gene is in the genes of the boolean generic network
#                         if gene in montagud_node_model:
#                             expr_max = rna_seq_data_max.loc[
#                                 rna_seq_data_max["gene_symbol"] == gene, "gene_max"
#                             ].iloc[0]
#                             gene_expr = rna_seq_data.loc[
#                                 (rna_seq_data["gene_symbol"] == gene)
#                                 & (rna_seq_data["model_id"] == model_id_in_file),
#                                 "rsem_tpm",
#                             ].iloc[0]
#                             prob_1 = min(max(gene_expr / expr_max, 0), 1)

#                             u_pattern = re.compile(
#                                 rf"\$u_{gene}\s*=\s*(0|1);", re.DOTALL
#                             )
#                             d_pattern = re.compile(
#                                 rf"\$d_{gene}\s*=\s*(0|1);", re.DOTALL
#                             )
#                             u_line = f"$u_{gene} = {prob_1:.4f};"
#                             d_line = f"$d_{gene} = {1 - prob_1:.4f};"
#                             content = re.sub(u_pattern, u_line, content)
#                             content = re.sub(d_pattern, d_line, content)
#                 # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
#                 modified_file_path = os.path.join(folder_models, filename)

#                 with open(modified_file_path, "w") as file:
#                     file.write(content)
#                 print(f"Modified and saved: {modified_file_path}")


# function used for validation and general pipeline
# def personalized_patients_proteins_cfgs(
#     df_melted_proteins,
#     montagud_node_model,
#     folder_models,
#     patients_ids,
#     table_proteins_patients,
#     drug_interest,
# ):
#     df_melted_proteins = df_melted_proteins[
#         df_melted_proteins["model_id"].isin(patients_ids)
#     ]

#     df_melted_proteins = df_melted_proteins[["model_id", "protein_symbol", "rsem_tpm"]]

#     # Open each file in the directory
#     # os.makedirs(original_data_dir, exist_ok=True)

#     # Directory where modified files will be saved
#     # os.makedirs(results_dir, exist_ok=True)

#     protein_data_max = df_melted_proteins
#     protein_data_max["protein_max"] = protein_data_max.groupby("protein_symbol")[
#         "rsem_tpm"
#     ].transform("max")

#     # Loop through each file in the directory
#     for filename in os.listdir(folder_models):
#         if not filename.endswith(".cfg"):
#             continue
#         file_path = os.path.join(folder_models, filename)
#         if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
#             # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
#             # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
#             model_id_in_file = os.path.splitext(filename)[0].replace(
#                 f"_{drug_interest}", ""
#             )
#             # Only proceed if model_id is in table_proteins_patients
#             if model_id_in_file in table_proteins_patients.index:
#                 with open(file_path, "r") as file:
#                     content = file.read()
#                     # rna_seq_data_filtered_filtered = rna_seq_data_filtered[rna_seq_data_filtered.index == model_id_in_file]
#                 if model_id_in_file in table_proteins_patients.index:
#                     row = table_proteins_patients.loc[model_id_in_file]
#                     # Only modify if this row corresponds to the file being processed
#                     high_proteins = row["High Protein Abundance"]
#                     low_proteins = row["Low Protein Abundance"]
#                     protein_high_list = [
#                         protein.strip()
#                         for protein in high_proteins.split(",")
#                         if protein.strip()
#                     ]
#                     protein_low_list = [
#                         protein.strip()
#                         for protein in low_proteins.split(",")
#                         if protein.strip()
#                     ]
#                     protein_list = list(set(protein_high_list + protein_low_list))
#                     for protein in protein_list:  # add condition if protein is in the proteins of the boolean generic network
#                         if protein in montagud_node_model:
#                             expr_max = protein_data_max.loc[
#                                 protein_data_max["protein_symbol"] == protein,
#                                 "protein_max",
#                             ].iloc[0]

#                             prot_expr = df_melted_proteins.loc[
#                                 (df_melted_proteins["protein_symbol"] == protein)
#                                 & (df_melted_proteins["model_id"] == model_id_in_file),
#                                 "rsem_tpm",
#                             ].iloc[0]

#                             prob_1 = min(max(prot_expr / expr_max, 0), 1)

#                             u_pattern = re.compile(
#                                 rf"\$u_{protein}\s*=\s*(0|1);", re.DOTALL
#                             )

#                             d_pattern = re.compile(
#                                 rf"\$d_{protein}\s*=\s*(0|1);", re.DOTALL
#                             )
#                             u_line = f"$u_{protein} = {prob_1:.4f};"
#                             d_line = f"$d_{protein} = {1 - prob_1:.4f};"
#                             content = re.sub(u_pattern, u_line, content)
#                             content = re.sub(d_pattern, d_line, content)

                    
#                             u_pattern = re.compile(
#                                 rf"\$u_{protein}\s*=\s*\d+\.?\d*;", re.DOTALL
#                             )
#                             d_pattern = re.compile(
#                                 rf"\$d_{protein}\s*=\s*\d+\.?\d*;", re.DOTALL
#                             )

#                             u_line = f"$u_{protein} = {prob_1:.4f};"
#                             d_line = f"$d_{protein} = {1 - prob_1:.4f};"

#                             content, num_u = re.subn(u_pattern, u_line, content)
#                             content, num_d = re.subn(d_pattern, d_line, content)

#                 # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
#                 modified_file_path = os.path.join(folder_models, filename)

#                 with open(modified_file_path, "w") as file:
#                     file.write(content)
#                 print(f"Modified and saved: {modified_file_path}")
