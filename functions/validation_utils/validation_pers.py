import maboss
import ginsim
import pandas as pd
import numpy as np
import mygene
import os
import shutil
import ast
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from pathlib import Path

from functions.generate_utils.pre_process_data.pre_process_montagud_nodes import (
    process_montagud_nodes,process_montagud_nodes_synonyms
)

from functions.generate_utils.create_generic_models.create_generic_patients_cfgs import (
    create_generic_patients_cfg_bnd_validation,
)
from functions.generate_utils.create_generic_models.update_phenotypes_generic_models import (
    generic_models_update_phenotypes,
)



from functions.generate_utils.create_person_models.tailor_cfgs_patients import (
    personalized_patients_genes_cfgs,
    personalized_patients_proteins_cfgs,
)
from functions.generate_utils.create_person_models.tailor_bnd_cnv import (
    tailor_bnd_cnv_validation,
)


def create_gleason_score_groups(phenotype_data, size_group):
    # stratify by gleason score-> create groups
    # create 3 groups: gleason score of 6, gleason score of 7, and of gleason score of > 8

    phenotype_data_filtered = phenotype_data[
        phenotype_data["sample_type"] != "Solid Tissue Normal"
    ][["sampleID", "gleason_score"]]

    conditions = [
        phenotype_data_filtered["gleason_score"].isin([6]),
        phenotype_data_filtered["gleason_score"].isin([7]),
        phenotype_data_filtered["gleason_score"].isin([8, 9, 10]),
    ]
    choices = ["low_aggressive", "middle_aggressive", "high_aggressive"]

    phenotype_data_filtered.loc[:, "Gleason_group"] = np.select(
        conditions, choices, default=""
    )

    sampled_df = phenotype_data_filtered.groupby(
        "Gleason_group", group_keys=False
    ).apply(lambda x: x.sample(n=min(len(x), size_group), random_state=42))
    patients_id = list(sampled_df["sampleID"])
    return sampled_df, patients_id




def pre_process_cnv(cnv_data, patients_id, all_montagud_nodes, synonyms_to_nodes_dict):
    # pre-process TCGA CNV data

    cnv_data_col = list(cnv_data.columns)
    common_col = list(set(cnv_data_col) & set(patients_id))
    col_keep = ["Gene Symbol"] + common_col
    cnv_data_filtered = cnv_data[col_keep]
    cnv_data_filtered = cnv_data_filtered[cnv_data_filtered["Gene Symbol"].isin(all_montagud_nodes)]

    for gene in cnv_data_filtered['Gene Symbol'].unique():
        if gene in synonyms_to_nodes_dict:
            cnv_data_filtered.loc[cnv_data_filtered['Gene Symbol'] == gene, 'Gene Symbol'] = synonyms_to_nodes_dict[gene]

    # remplace the cnv symbol column names by its synonyms in the synonyms_maps dictionary
    # cnv_data_filtered["Gene Symbol"] = cnv_data_filtered["Gene Symbol"].apply(
    #     lambda x: synonyms_maps.get(x, x)
    # )


    df_melted_cnv = cnv_data_filtered.melt(
        id_vars=["Gene Symbol"], var_name="samples_id", value_name="expression_value"
    )

    df_melted_cnv["Gene Symbol"] = df_melted_cnv["Gene Symbol"].str.split("|").str[0]

    df_melted_cnv = df_melted_cnv.rename(
        columns={
            "samples_id": "model_id",
            "Gene Symbol": "gene_symbol",
            "expression_value": "rsem_tpm",
        }
    )

    conditions = [
        df_melted_cnv["rsem_tpm"].isin([-1, -2]),
        df_melted_cnv["rsem_tpm"].isin([0]),
        df_melted_cnv["rsem_tpm"].isin([1, 2]),
    ]
    choices = ["Loss", "Normal", "Gain"]
    df_melted_cnv.loc[:, "effect"] = np.select(conditions, choices, default="")

    
    df_melted_cnv.to_csv("data/TCGA_data/prostate/filtered_data/cnv_samples_table.csv")
    return df_melted_cnv




def pre_process_genes(genes_data, patients_id, all_montagud_nodes, synonyms_to_nodes_dict):
    # pre-process genes data

    genes_data_col = list(genes_data.columns)
    common_col = list(set(genes_data_col) & set(patients_id))
    col_keep = ["sample"] + common_col
    genes_data_filtered = genes_data[col_keep]
    genes_data_filtered = genes_data_filtered[genes_data_filtered["sample"].isin(all_montagud_nodes)]

    # Special cases handling: Create both mTORC1 and mTORC2 from MTOR data
    syn_dict = {'MTOR': ['mTORC1', 'mTORC2'], 'MYC': ['MYC', 'MYC_MAX'], 'PIK3CA': ['PI3K', 'PIP3'], 'LDHA': ['LDHA', 'Lactic_acid'], 'ERG': ['AR_ERG', 'ERG']}
    # Get all keys
    list_genes_duplicates = syn_dict.keys()


    # remplace the gene symbol column names by its synonyms in the synonyms_maps dictionary
    for gene in genes_data_filtered['sample'].unique():
        if gene in synonyms_to_nodes_dict and gene not in list_genes_duplicates:
            genes_data_filtered.loc[genes_data_filtered['sample'] == gene, 'sample'] = synonyms_to_nodes_dict[gene]

    for duplicate_gene in list_genes_duplicates:
        gene_duplicate_data = genes_data_filtered[genes_data_filtered['sample'] == duplicate_gene].copy()
        if not gene_duplicate_data.empty:
            # Create mTORC1 data
            gene_duplicate_1_data = gene_duplicate_data.copy()
            gene_duplicate_1_data['sample'] = syn_dict[duplicate_gene][0]
            
            # Create mTORC2 data  
            gene_duplicate_2_data = gene_duplicate_data.copy()
            gene_duplicate_2_data['sample'] = syn_dict[duplicate_gene][1]
            
            # Remove original MTOR and add both complexes
            genes_data_filtered = genes_data_filtered[genes_data_filtered['sample'] != duplicate_gene]
            genes_data_filtered = pd.concat([genes_data_filtered, gene_duplicate_1_data, gene_duplicate_2_data], ignore_index=True)

            print(f" Duplicated {duplicate_gene}: {syn_dict[duplicate_gene][0]} ({len(gene_duplicate_1_data)} rows) + {syn_dict[duplicate_gene][1]} ({len(gene_duplicate_2_data)} rows)")


    df_melted_gene = genes_data_filtered.melt(
        id_vars=["sample"],  # columns to keep fixed
        var_name="samples_id",  # name for the variable column (sample IDs)
        value_name="expression_value",  # name for the values
    )

    df_melted_gene["sample"] = df_melted_gene["sample"].str.split("|").str[0]

    df_melted_gene = df_melted_gene.rename(
        columns={
            "samples_id": "model_id",
            "sample": "gene_symbol",
            "expression_value": "rsem_tpm",
        }
    )
 
    df_melted_gene.to_csv(
        "data/TCGA_data/prostate/filtered_data/genes_samples_table.csv"
    )

    return df_melted_gene


def process_proteins_validation(
    proteins_data, patients_id, montagud_nodes, synonyms_to_nodes_dict
):
    proteins_data["sample"] = proteins_data["sample"].str.rsplit("-", n=2).str[0]
    proteins_data["sample"] = proteins_data["sample"].str.replace("_", "", regex=False)
    proteins_data["sample"] = proteins_data["sample"].str.replace("-", "", regex=False)


    montagud_nodes = [str(node) for node in montagud_nodes if pd.notna(node) and node != '']
    montagud_nodes = [node for node in montagud_nodes if node != 'nan']

    pattern = "|".join(montagud_nodes)

    # Filter rows where 'sample' contains any name from the list
    # check??
    proteins_data["sample"] = proteins_data["sample"].astype(str)

    proteins_data = proteins_data[
        proteins_data["sample"].str.contains(pattern, case=False, na=False)
    ]

    proteins_data = proteins_data.dropna(how="all", subset=proteins_data.columns[1:])

    mods = ["_PS", "_PT", "_PY"]
    proteins_data = proteins_data[
        ~proteins_data["sample"].apply(lambda p: any(mod in p for mod in mods))
    ]

    # remplace the protein symbol column names by its synonyms in the proteins_synonyms_maps dictionary
    for protein in proteins_data['sample'].unique():
        if protein in synonyms_to_nodes_dict:
            proteins_data.loc[proteins_data['sample'] == protein, 'sample'] = synonyms_to_nodes_dict[protein]



    # Keep only the proteins present in the montagud list

    proteins_data_col = list(proteins_data.columns)
    common_col = list(set(proteins_data_col) & set(patients_id))
    col_keep = ["sample"] + common_col
    proteins_data_filtered = proteins_data[col_keep]

    df_melted_protein = proteins_data_filtered.melt(
        id_vars=["sample"],  # columns to keep fixed
        var_name="samples_id",  # name for the variable column (sample IDs)
        value_name="expression_value",  # name for the values
    )

    df_melted_protein["sample"] = df_melted_protein["sample"].str.split("|").str[0]

    df_melted_protein = df_melted_protein.rename(
        columns={
            "samples_id": "model_id",
            "sample": "protein_symbol",
            "expression_value": "rsem_tpm",
        }
    )
  
    def replace_with_base_name(protein_name):
        for base in montagud_nodes:
            if protein_name.startswith(base):
                return base
        return protein_name  # if no match found, keep original

    # Assuming your dataframe is df and column to replace is 'protein_symbol'
    df_melted_protein["protein_symbol"] = df_melted_protein["protein_symbol"].apply(
        replace_with_base_name
    )

    df_melted_protein = df_melted_protein[
        df_melted_protein["protein_symbol"].isin(montagud_nodes)
    ]
    df_melted_protein["protein_symbol"] = df_melted_protein[
        "protein_symbol"
    ].str.replace("_", "", regex=False)
    df_melted_protein = df_melted_protein[df_melted_protein["rsem_tpm"].notna()]
    # df_melted_protein.to_csv(
    #     "data/TCGA_data/prostate/filtered_data/proteins_samples_table.csv"
    # )
    return df_melted_protein




#   MAIN function



def pre_process_data_validation(
    montagud_original_data_df,
    nodes_montagud_synonyms,
    genes_data,
    cnv_data,
    size_group,
    phenotype_data,
    proteins_data,
    
):
    # get the patients
    patients_groups, patients_id = create_gleason_score_groups(
        phenotype_data, size_group
    )

    # pre process montagud nodes

    montagud_node_synonyms, synonyms_to_nodes_dict = process_montagud_nodes_synonyms(nodes_montagud_synonyms)


    montagud_node_model, all_montagud_nodes = process_montagud_nodes(
        montagud_original_data_df, montagud_node_synonyms
    )






    # pre process CNV data
    df_melted_cnv = pre_process_cnv(
        cnv_data, patients_id, all_montagud_nodes, synonyms_to_nodes_dict
    )


    # pre process genes data
    df_melted_gene = pre_process_genes(
        genes_data, patients_id, all_montagud_nodes, synonyms_to_nodes_dict
    )



    # pre process proteins data

    df_melted_protein = process_proteins_validation(
        proteins_data, patients_id, all_montagud_nodes, synonyms_to_nodes_dict
    )



    return (
        patients_groups,
        patients_id,
        montagud_node_model,
        df_melted_gene,
        df_melted_cnv,
        df_melted_protein,
    )



def create_pers_models_generic(
    folder_validation,
    original_model_cfg,
    original_model_bnd,
    patients_id,
    tissue,
    type_models,
    phenotype_interest,
    patients_groups,
    montagud_node_model,
    df_melted_gene,
    df_melted_cnv,
    df_melted_proteins,
    amplif_factor,
    context_label,
):
    model_generic = "analysis/validation/generic"
    os.makedirs(model_generic, exist_ok=True)

    original_model_bnd_copy = os.path.join(
        model_generic, os.path.basename(original_model_bnd)
    )
    original_model_cfg_copy = os.path.join(
        model_generic, os.path.basename(original_model_cfg)
    )

    shutil.copyfile(original_model_bnd, original_model_bnd_copy)
    shutil.copyfile(original_model_cfg, original_model_cfg_copy)

    create_generic_patients_cfg_bnd_validation(
        original_model_cfg_copy,
        original_model_bnd_copy,
        folder_validation,
        patients_id,
        patients_groups,
        tissue,
    )

    groups = list(set(patients_groups["Gleason_group"]))

    for group in groups:
        folder_validation_temp = os.path.join(folder_validation, group)
        # # update phenotypes in generic models
        generic_models_update_phenotypes(phenotype_interest, folder_validation_temp)

        # # personalize the networks with genes
        # # tissue instead of drug name

        if type_models == "genes_models":
            print("Personalizing genes models...")
            personalized_patients_genes_cfgs(
                df_melted_gene,
                montagud_node_model,
                folder_validation_temp,
                amplif_factor,
                context_label = context_label,
            )
        else:
            personalized_patients_proteins_cfgs(
                df_melted_proteins,
                montagud_node_model,
                folder_validation_temp,
                context_label = context_label,
            )

        # # personalize the networks with CNV
        tailor_bnd_cnv_validation(df_melted_cnv, folder_validation_temp, context_label)



# def pre_process_montagud_nodes(montagud_data, name_maps, nodes_to_add):
#     # transform nodes_to_add to list if single element
#     if isinstance(nodes_to_add, str):
#         nodes_to_add = [nodes_to_add]

#     # Create list of genes of interest (in Montagud data)
#     montagud_nodes = list(
#         set(montagud_data["Target node"].tolist() + montagud_data["Source"].tolist())
#     )
#     montagud_nodes = [node for node in montagud_nodes if node != "0/1"]
#     # Apply mapping
#     montagud_nodes = [name_maps.get(x, x) for x in montagud_nodes]
#     # Add any additional nodes
#     montagud_nodes = list(set(montagud_nodes + nodes_to_add))

#     return montagud_nodes