import pandas as pd
import numpy as np
import os
import re
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import spearmanr
import gseapy as gp


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

            print(f" Duplicated {duplicate_gene}: {syn_dict[duplicate_gene][0]} ({len(gene_duplicate_1_data)} rows) + {syn_dict[duplicate_gene][1]} ({len(gene_duplicate_2_data)} rows)")

        

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
    cnv_data_filtered['cn_category'] = np.select(conditions, choices, default='other')

    # filter out Neural Category
    cnv_data_filtered = cnv_data_filtered[~cnv_data_filtered["cn_category"].isin(["normal", 'loss', 'gain'])]




    # keep only the patients id and montagud nodes
    cnv_data_filtered = cnv_data_filtered[
        cnv_data_filtered["gene_symbol"].isin(all_montagud_nodes)
    ]
    
    cnv_data_filtered = cnv_data_filtered[
        ["model_id", "gene_symbol", "total_copy_number", "cn_category"]
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






# Create personalized boolean networks


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
        patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
        patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

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
    print(
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
                        print(f"Updated {node}.is_internal=FALSE in {filename}")
                        modified = True
                else:
                    # Add new assignment at the end if not present
                    content += f"\n{node}.is_internal=FALSE;"
                    print(f"Added {node}.is_internal=FALSE to {filename}")
                    modified = True

            
            # Save only if modified
            if modified:
                modified_file_path = os.path.join(folder_models, filename)
                with open(modified_file_path, "w") as file:
                    file.write(content)
                    print(f"Modified and saved: {modified_file_path}")





def normalize_rna_seq_efficient(rna_seq_data_models_filtered):
    """
    Efficient min-max normalization using pandas groupby
    """
    rna_normalized = rna_seq_data_models_filtered.copy()
    
    def minmax_normalize_group(group):
        """Apply min-max normalization to a group (gene)"""
        scaler = MinMaxScaler()
        group['rsem_tpm_normalized'] = scaler.fit_transform(group[['rsem_tpm']]).flatten()
        return group
    
    print("Applying min-max normalization per gene...")
    
    # Group by gene and apply normalization
    rna_normalized = rna_normalized.groupby('gene_symbol').apply(minmax_normalize_group).reset_index(drop=True)
    
    print("Normalization completed!")
    print(f"Original column: 'rsem_tpm', Normalized column: 'rsem_tpm_normalized'")
    return rna_normalized



def personalized_patients_genes_cfgs(
    rna_seq_data_models_filtered,
    montagud_node_model,
    folder_models,
    amplif_factor,
):
    # Apply the normalization
    rna_seq_data_normalized = normalize_rna_seq_efficient(rna_seq_data_models_filtered)


    print("==== Processing all genes per patient ===")
    # Loop through each file in the directory
    for filename in os.listdir(folder_models):
        if not filename.endswith('.cfg'):
            continue

        file_path = os.path.join(folder_models, filename)
        print(f"\nProcessing file: {filename}")

        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0]
            
            with open(file_path, "r") as file:
                content = file.read()
      

            for gene in montagud_node_model:
                gene_data_normalized = rna_seq_data_normalized[
                        (rna_seq_data_normalized["gene_symbol"] == gene) &
                        (rna_seq_data_normalized["model_id"] == model_id_in_file)
                    ]['rsem_tpm_normalized']
                
                if gene_data_normalized.empty:
                        print(f"  No data for {gene}")
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
        print(f"Modified and saved: {modified_file_path}")






def tailor_bnd_cnv_cm(cnv_data_filtered, folder_models):
    gain_group = cnv_data_filtered[
        cnv_data_filtered["cn_category"]=='amplification'
    ]

    loss_group = cnv_data_filtered[
        cnv_data_filtered["cn_category"]=='deletion'
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
            print(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(folder_models, files_categ[0])

        if bnd_file is None:
            print(f"No .bnd file found for patient: {patient}")
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
            print(f"🔍 Processing patient {patient}, gene: {gene}")
            
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
                print(
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
                print(f"Patient {patient} with gene {gene} not found in gain or loss CNV groups.")
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
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            print(f"{patient_id}: CNV — no nodes modified")







def correlate_boolean_predictions_with_gene_signatures(proba_phenotype, hallmark, phenotype, rna_seq_data, patients_ids):
    # Preprocessing -> get gene expression matrix
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['model_id'].isin(patients_ids)]


    rna_seq_data_filtered = rna_seq_data_filtered[['model_id', 'gene_symbol', 'rsem_tpm']]


    # remove duplicates (compute mean)
    rna_seq_data_deduplicated = rna_seq_data_filtered.groupby(['model_id', 'gene_symbol']).agg({
        'rsem_tpm': 'mean'
    }).reset_index()

    deduplicated_gene_sample = rna_seq_data_deduplicated.duplicated(subset=['model_id', 'gene_symbol'])

    # Reshape the DataFrame: genes as rows (index), samples as columns, expression as values
    gene_expression_matrix = rna_seq_data_deduplicated.pivot(
        index='gene_symbol',     # Genes become row index
        columns='model_id',      # Model IDs become columns
        values='rsem_tpm'        # Expression values fill the matrix
    )

    patients_in_matrix = gene_expression_matrix.columns.tolist()



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
        print(f"  ✅ {hallmark} score: {phenotype_score:.3f}")

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
    print('correlation_data', correlation_data)
    results_corr_df = pd.DataFrame(index = phenotype_cols, columns = ['correlation_value', 'p_val'])


    for col in phenotype_cols:
        corr, pval = spearmanr(correlation_data[col], correlation_data[f'{hallmark}_gene_signature_score'])
        results_corr_df.loc[col, 'correlation_value'] = corr
        results_corr_df.loc[col, 'p_val'] = pval


    return results_corr_df, correlation_data
    