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


from functions.generate_utils.create_generic_models.create_generic_patients_cfgs import (
    create_generic_patients_cfg_bnd_validation,
)
from functions.generate_utils.create_generic_models.update_phenotypes_generic_models import (
    generic_models_update_phenotypes,
)
from functions.generate_utils.pre_process_data.pre_process_genes import (
    create_table_rna_seq_patients,
)

from functions.generate_utils.pre_process_data.pre_process_proteins import (
    process_proteins_validation,
    create_table_proteins_patients,
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


def pre_process_cnv(cnv_data, patients_id, montagud_nodes, synonyms_maps):
    # pre-process TCGA CNV data

    cnv_data_col = list(cnv_data.columns)
    common_col = list(set(cnv_data_col) & set(patients_id))
    col_keep = ["Gene Symbol"] + common_col
    cnv_data_filtered = cnv_data[col_keep]

    # remplace the cnv symbol column names by its synonyms in the synonyms_maps dictionary
    cnv_data_filtered["Gene Symbol"] = cnv_data_filtered["Gene Symbol"].apply(
        lambda x: synonyms_maps.get(x, x)
    )

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

    df_melted_cnv = df_melted_cnv[df_melted_cnv["gene_symbol"].isin(montagud_nodes)]
    df_melted_cnv.to_csv("data/TCGA_data/prostate/filtered_data/cnv_samples_table.csv")
    return df_melted_cnv


def pre_process_genes(genes_data, patients_id, synonyms_maps):
    # pre-process genes data

    genes_data_col = list(genes_data.columns)
    common_col = list(set(genes_data_col) & set(patients_id))
    col_keep = ["sample"] + common_col
    genes_data_filtered = genes_data[col_keep]

    # remplace the gene symbol column names by its synonyms in the synonyms_maps dictionary
    genes_data_filtered["sample"] = genes_data_filtered["sample"].apply(
        lambda x: synonyms_maps.get(x, x)
    )

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
    df_melted_gene["gene_symbol"] = df_melted_gene["gene_symbol"].str.upper()
    df_melted_gene["gene_symbol"] = df_melted_gene["gene_symbol"].str.replace(
        "_", "", regex=False
    )
    df_melted_gene.to_csv(
        "data/TCGA_data/prostate/filtered_data/genes_samples_table.csv"
    )
    table_rna_seq_patients = create_table_rna_seq_patients(df_melted_gene)

    return df_melted_gene, table_rna_seq_patients


def pre_process_montagud_nodes(montagud_data, name_maps, nodes_to_add):
    # transform nodes_to_add to list if single element
    if isinstance(nodes_to_add, str):
        nodes_to_add = [nodes_to_add]

    # Create list of genes of interest (in Montagud data)
    montagud_nodes = list(
        set(montagud_data["Target node"].tolist() + montagud_data["Source"].tolist())
    )
    montagud_nodes = [node for node in montagud_nodes if node != "0/1"]
    montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]
    # Apply mapping
    montagud_nodes = [name_maps.get(x, x) for x in montagud_nodes]
    # Add any additional nodes
    montagud_nodes = list(set(montagud_nodes + nodes_to_add))

    return montagud_nodes


#   MAIN function
def pre_process_data_validation(
    size_group,
    phenotype_data,
    cnv_data,
    montagud_data,
    name_maps,
    nodes_to_add,
    genes_data,
    proteins_data,
    synonyms_maps,
):
    # get the patients
    patients_groups, patients_id = create_gleason_score_groups(
        phenotype_data, size_group
    )

    # pre process montagud nodes
    montagud_nodes = pre_process_montagud_nodes(montagud_data, name_maps, nodes_to_add)

    # pre process CNV data
    df_melted_cnv = pre_process_cnv(
        cnv_data, patients_id, montagud_nodes, synonyms_maps
    )

    # pre process genes data
    df_melted_gene, table_rna_seq_patients = pre_process_genes(
        genes_data, patients_id, synonyms_maps
    )

    # pre process proteins data
    df_melted_protein = process_proteins_validation(
        proteins_data, montagud_nodes, synonyms_maps, patients_id
    )
    table_proteins_patients = create_table_proteins_patients(df_melted_protein)

    return (
        patients_groups,
        patients_id,
        montagud_nodes,
        df_melted_gene,
        table_rna_seq_patients,
        df_melted_cnv,
        df_melted_protein,
        table_proteins_patients,
    )


def create_pers_models_generic(
    folder_validation,
    original_model_cfg,
    original_model_bnd,
    patients_id,
    tissue,
    name_maps,
    type_models,
    phenotype_interest,
    patients_groups,
    montagud_nodes,
    df_melted_gene,
    df_melted_cnv,
    table_rna_seq_patients,
    df_melted_proteins,
    table_proteins_patients,
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
        name_maps,
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
                montagud_nodes,
                folder_validation_temp,
                patients_id,
                table_rna_seq_patients,
                context_label,
            )
        else:
            personalized_patients_proteins_cfgs(
                df_melted_proteins,
                montagud_nodes,
                folder_validation_temp,
                patients_id,
                table_proteins_patients,
                context_label,
            )

        # # personalize the networks with CNV
        tailor_bnd_cnv_validation(df_melted_cnv, folder_validation_temp, context_label)
