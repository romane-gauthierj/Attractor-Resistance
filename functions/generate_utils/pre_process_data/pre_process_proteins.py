# Pre-process proteins data (only id in the patients id list  and proteins in montagud list)



# def classify_expression(z):
#     if z > 1:
#         return "high"
#     elif z < -1:
#         return "low"
#     else:
#         return "normal"


# def create_table_proteins_patients(proteins_data):
#     # filter only data with patient id and the nodes of the montagud_data

#     proteins_data["z_score"] = proteins_data.groupby("protein_symbol")[
#         "rsem_tpm"
#     ].transform(lambda x: (x - x.mean()) / x.std())

#     proteins_data["protein_expression_level"] = proteins_data["z_score"].apply(
#         classify_expression
#     )

#     proteins_data = proteins_data[
#         ["model_id", "protein_symbol", "protein_expression_level"]
#     ]
#     proteins_data.rename(columns={"protein_symbol": "protein_name"}, inplace=True)
#     proteins_data = proteins_data[
#         proteins_data["protein_expression_level"].isin(["low", "high"])
#     ]

#     table_proteins_patients = proteins_data.pivot_table(
#         index="model_id",
#         columns="protein_expression_level",
#         values="protein_name",
#         aggfunc=lambda x: ", ".join(sorted(set(x))),
#     ).fillna("-")

#     table_proteins_patients = table_proteins_patients.rename(
#         columns={"low": "Low Protein Abundance", "high": "High Protein Abundance"}
#     )

#     return table_proteins_patients




def process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict):

    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]


    # first filter the large datasets with the nodes of the montagud as well as their synonyms
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['gene_symbol'].isin(all_montagud_nodes)]

   
    rna_seq_data_filtered = rna_seq_data_filtered[["model_id", "gene_symbol", "rsem_tpm"]]

    return rna_seq_data_filtered




def process_proteins(patients_ids, proteins_data, all_montagud_nodes, synonyms_to_nodes_dict):
    # explode the 'symbol' column if it contains multiple values separated by commas or semicolons

    proteins_data["symbol"] = proteins_data["symbol"].astype(str)
    proteins_data["symbol"] = proteins_data["symbol"].str.split(r"[;,]")
    proteins_data = proteins_data.explode("symbol")

    proteins_data["symbol"] = proteins_data["symbol"].str.strip()

    proteins_data_patient_id = list(set(proteins_data["model_id"]))

    common_col = list(set(proteins_data_patient_id) & set(patients_ids))
    col_keep = ["symbol"] + common_col

    proteins_data_patient_id_filtered = proteins_data[
        proteins_data["model_id"].isin(col_keep)
    ]


    proteins_names_list = list(set(proteins_data_patient_id_filtered["symbol"]))

    common_proteins = list(set(proteins_names_list) & set(all_montagud_nodes))

    # filter the big datasets with the all montagud nodes proteins 
    proteins_data_patient_id_filtered = proteins_data_patient_id_filtered[
        proteins_data_patient_id_filtered["symbol"].isin(common_proteins)
    ]



    # remplace all the proteins synonyms to the unique protein name of the models 
    for protein in proteins_data_patient_id_filtered["symbol"]:
        if protein in synonyms_to_nodes_dict:
            proteins_data_patient_id_filtered.loc[proteins_data_patient_id_filtered['symbol'] == protein, 'symbol'] = synonyms_to_nodes_dict[protein]


    #  uniprot_id	model_id	model_name	protein_intensity	zscore	symbol
    proteins_data_patient_id_filtered = proteins_data_patient_id_filtered[
        ["model_id", "protein_intensity", "zscore", "symbol"]
    ]

    proteins_data_patient_id_filtered = proteins_data_patient_id_filtered.rename(
        columns={
            "symbol": "protein_symbol",
            "protein_intensity": "rsem_tpm",
            "zscore": "z_score",
        }
    )

    return proteins_data_patient_id_filtered














# def process_proteins(patients_ids, proteins_data, montagud_nodes, synonyms_maps):
#     # explode the 'symbol' column if it contains multiple values separated by commas or semicolons

#     proteins_data["symbol"] = proteins_data["symbol"].astype(str)
#     proteins_data["symbol"] = proteins_data["symbol"].str.split(r"[;,]")
#     proteins_data = proteins_data.explode("symbol")

#     proteins_data["symbol"] = proteins_data["symbol"].str.strip()

#     proteins_data_patient_id = list(set(proteins_data["model_id"]))

#     common_col = list(set(proteins_data_patient_id) & set(patients_ids))
#     col_keep = ["symbol"] + common_col

#     proteins_data_patient_id_filtered = proteins_data[
#         proteins_data["model_id"].isin(col_keep)
#     ]

#     # remplace the protein symbol column names by its synonyms in the proteins_synonyms_maps dictionary
#     proteins_data_patient_id_filtered["symbol"] = proteins_data_patient_id_filtered[
#         "symbol"
#     ].apply(lambda x: synonyms_maps.get(x, x))

#     proteins_names_list = list(set(proteins_data_patient_id_filtered["symbol"]))

#     common_proteins = list(set(proteins_names_list) & set(montagud_nodes))

#     proteins_data_patient_id_filtered = proteins_data_patient_id_filtered[
#         proteins_data_patient_id_filtered["symbol"].isin(common_proteins)
#     ]

#     #  uniprot_id	model_id	model_name	protein_intensity	zscore	symbol
#     proteins_data_patient_id_filtered = proteins_data_patient_id_filtered[
#         ["model_id", "protein_intensity", "zscore", "symbol"]
#     ]

#     proteins_data_patient_id_filtered = proteins_data_patient_id_filtered.rename(
#         columns={
#             "symbol": "protein_symbol",
#             "protein_intensity": "rsem_tpm",
#             "zscore": "z_score",
#         }
#     )

#     return proteins_data_patient_id_filtered
