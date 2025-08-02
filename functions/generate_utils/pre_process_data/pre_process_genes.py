# Cell model passport - RNA seq expression

import numpy as np
import pandas as pd
import os


def identify_genes_synonyms(data_synonyms, uniprot_data, montagud_nodes):
    data_synonyms = data_synonyms[
        ~(data_synonyms["Gene name"].isna() & data_synonyms["Gene Synonym"].isna())
    ]

    # added - check ??
    data_synonyms["Gene name"] = data_synonyms["Gene name"].astype(str)
    data_synonyms["Gene Synonym"] = data_synonyms["Gene Synonym"].astype(str)


    data_synonyms["Gene name"] = data_synonyms["Gene name"].str.upper()
    data_synonyms["Gene Synonym"] = data_synonyms["Gene Synonym"].str.upper()

    data_synonyms = data_synonyms[
        data_synonyms["Gene name"].isin(montagud_nodes)
        | data_synonyms["Gene Synonym"].isin(montagud_nodes)
    ]

    def match_montagud_node(row, montagud_genes):
        if row["Gene name"] in montagud_genes:
            return row["Gene name"]
        elif row["Gene Synonym"] in montagud_genes:
            return row["Gene Synonym"]
        else:
            return None

    # Apply the function to each row
    data_synonyms["montagud_node"] = data_synonyms.apply(
        lambda row: match_montagud_node(row, montagud_nodes), axis=1
    )

    merged_df = pd.merge(data_synonyms, uniprot_data, on="Gene stable ID")
    merged_df["Gene Synonym"] = merged_df["Gene Synonym"].str.replace(
        r"[_-]", "", regex=True
    )
    merged_df = merged_df[["Gene name_x", "Gene Synonym", "montagud_node"]]

    return merged_df




# def process_genes(patients_ids, montagud_nodes, rna_seq_data, synonyms_maps):
#     rna_seq_data["gene_symbol_upper"] = rna_seq_data["gene_symbol"].str.upper()
#     rna_seq_data = rna_seq_data[rna_seq_data["gene_symbol_upper"].isin(montagud_nodes)]
#     rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]
#     rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]

#     # remplace the gene symbol column names by its synonyms in the synonyms_maps dictionary
#     rna_seq_data["gene_symbol"] = rna_seq_data["gene_symbol"].apply(
#         lambda x: synonyms_maps.get(x, x)
#     )
#     return rna_seq_data

    
# def process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict):

#     rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]


#     # first filter the large datasets with the nodes of the montagud as well as their synonyms
#     rna_seq_data_filtered = rna_seq_data[rna_seq_data['gene_symbol'].isin(all_montagud_nodes)]

#     # apply the synonym mapping 
#     for gene in rna_seq_data_filtered['gene_symbol'].unique():
#         if gene in synonyms_to_nodes_dict:
#             rna_seq_data_filtered.loc[rna_seq_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]

#     rna_seq_data_models_filtered = rna_seq_data_filtered[["model_id", "gene_symbol", "rsem_tpm"]]

#     return rna_seq_data_models_filtered

    



def process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict):
    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['gene_symbol'].isin(all_montagud_nodes)]


  # Special cases handling: Create both mTORC1 and mTORC2 from MTOR data
    syn_dict = {'MTOR': ['mTORC1', 'mTORC2'], 'MYC': ['MYC', 'MYC_MAX'], 'PIK3CA': ['PI3K', 'PIP3'], 'LDHA': ['LDHA', 'Lactic_acid'], 'ERG': ['AR_ERG', 'ERG']}
    # Get all keys
    list_genes_duplicates = syn_dict.keys()

    # apply the synonym mapping 
    for gene in rna_seq_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict and gene not in list_genes_duplicates:
            rna_seq_data_filtered.loc[rna_seq_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]


    for duplicate_gene in list_genes_duplicates:
        gene_duplicate_data = rna_seq_data_filtered[rna_seq_data_filtered['gene_symbol'] == duplicate_gene].copy()
        if not gene_duplicate_data.empty:
            # Create mTORC1 data
            gene_duplicate_1_data = gene_duplicate_data.copy()
            gene_duplicate_1_data['gene_symbol'] = syn_dict[duplicate_gene][0]
            
            # Create mTORC2 data  
            gene_duplicate_2_data = gene_duplicate_data.copy()
            gene_duplicate_2_data['gene_symbol'] = syn_dict[duplicate_gene][1]
            
            # Remove original MTOR and add both complexes
            rna_seq_data_filtered = rna_seq_data_filtered[rna_seq_data_filtered['gene_symbol'] != duplicate_gene]
            rna_seq_data_filtered = pd.concat([rna_seq_data_filtered, gene_duplicate_1_data, gene_duplicate_2_data], ignore_index=True)

            print(f" Duplicated {duplicate_gene}: {syn_dict[duplicate_gene][0]} ({len(gene_duplicate_1_data)} rows) + {syn_dict[duplicate_gene][1]} ({len(gene_duplicate_2_data)} rows)")

        

    rna_seq_data_models_filtered = rna_seq_data_filtered[["model_id", "gene_symbol", "rsem_tpm"]]
    return rna_seq_data_models_filtered





# def classify_expression(z):
#     if z > 1:
#         return "high"
#     elif z < -1:
#         return "low"
#     else:
#         return "normal"


# def create_table_rna_seq_patients(rna_seq_data_filtered):
#     # filter only data with patient id and the nodes of the montagud_data

#     rna_seq_data_filtered["z_score"] = rna_seq_data_filtered.groupby("gene_symbol")["rsem_tpm"].transform(
#         lambda x: (x - x.mean()) / x.std()
#     )

#     rna_seq_data_filtered["gene_expression_level"] = rna_seq_data_filtered["z_score"].apply(
#         classify_expression
#     )

#     rna_seq_data_filtered = rna_seq_data_filtered[["model_id", "gene_symbol", "gene_expression_level"]]
#     rna_seq_data_filtered.rename(columns={"gene_symbol": "gene_name"}, inplace=True)
#     rna_seq_data_filtered = rna_seq_data_filtered[
#         rna_seq_data_filtered["gene_expression_level"].isin(["low", "high"])
#     ]

#     table_rna_seq_patients = rna_seq_data_filtered.pivot_table(
#         index="model_id",
#         columns="gene_expression_level",
#         values="gene_name",
#         aggfunc=lambda x: ", ".join(sorted(set(x))),
#     ).fillna("-")

#     table_rna_seq_patients = table_rna_seq_patients.rename(
#         columns={"low": "Low Gene Expression", "high": "High Gene Expression"}
#     )

#     return table_rna_seq_patients
