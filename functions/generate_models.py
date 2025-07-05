import os

from functions.generate_utils.identification_patients.get_patients_sens_res import (
    get_patients,
)
from functions.generate_utils.pre_process_data.pre_process_cnv import preprocess_cnv
from functions.generate_utils.pre_process_data.pre_process_genes import (
    process_genes,
    create_table_rna_seq_patients,
)
from functions.generate_utils.pre_process_data.pre_process_montagud_nodes import (
    process_montagud_nodes,
)
from functions.generate_utils.create_generic_models.update_phenotypes_generic_models import (
    generic_models_update_phenotypes,
)
from functions.generate_utils.create_generic_models.create_generic_patients_cfgs import (
    create_generic_patients_cfgs_bnds,
)
from functions.generate_utils.create_person_models.tailor_cfgs_patients import (
    personalized_patients_genes_cfgs,
    personalized_patients_proteins_cfgs,
)
from functions.generate_utils.pre_process_data.pre_process_genes import (
    create_table_rna_seq_patients,
    process_genes,
)
from functions.generate_utils.create_person_models.tailor_bnd_cnv import (
    tailor_bnd_cnv_cm,
)

from functions.analysis_utils.genes_intervention.pers_interventions import (
    tailor_bnd_genes_intervention,
)

from functions.generate_utils.pre_process_data.pre_process_proteins import (
    process_proteins,
    create_table_proteins_patients,
)


def pre_process_re(
    montagud_data,
    rna_seq_data,
    cnv_data,
    number_patients,
    drug_data,
    annotations_models,
    drug_interest,
    proteins_data,
    type_models,
    name_montagud_maps,
    nodes_to_add,
    synonyms_maps,
    tissue_interest=None,
    tissue_remove=None,
    node_to_remove=None,
):
    top_resistant_ids, top_sensitive_ids, top_healthy_ids, drug_data_filtered = (
        get_patients(
            number_patients,
            drug_data,
            annotations_models,
            drug_interest,
            tissue_interest=None,
            tissue_remove=tissue_remove,
        )
    )

    patients_ids = top_resistant_ids + top_sensitive_ids

    # preprocess montagud nodes
    montagud_nodes = process_montagud_nodes(
        montagud_data, name_montagud_maps, nodes_to_add
    )

    # preprocess rna seq data
    rna_seq_data_filtered = process_genes(
        patients_ids, montagud_nodes, rna_seq_data, synonyms_maps
    )
    table_rna_seq_patients = create_table_rna_seq_patients(rna_seq_data_filtered)

    # pre process proteins data
    df_melted_protein = process_proteins(
        proteins_data, montagud_nodes, synonyms_maps, patients_ids
    )

    table_proteins_patients = create_table_proteins_patients(df_melted_protein)

    if type_models == "genes_models":
        top_resistant_ids = list(
            set(table_rna_seq_patients.index) & set(top_resistant_ids)
        )
        top_sensitive_ids = list(
            set(table_rna_seq_patients.index) & set(top_sensitive_ids)
        )
    else:
        # if proteins models
        top_resistant_ids = list(
            set(table_proteins_patients.index) & set(top_resistant_ids)
        )
        top_sensitive_ids = list(
            set(table_proteins_patients.index) & set(top_sensitive_ids)
        )

    patients_ids = top_resistant_ids + top_sensitive_ids

    # Save IDs to txt files
    with open("analysis/{}/top_resistant_ids.txt".format(drug_interest), "w") as f:
        for pid in top_resistant_ids:
            f.write(f"{pid}\n")
    with open("analysis/{}/top_sensitive_ids.txt".format(drug_interest), "w") as f:
        for pid in top_sensitive_ids:
            f.write(f"{pid}\n")

    # preprocess cnv data
    cnv_data_filtered = preprocess_cnv(
        cnv_data, montagud_nodes, patients_ids, synonyms_maps
    )

    return (
        top_resistant_ids,
        top_sensitive_ids,
        top_healthy_ids,
        montagud_nodes,
        rna_seq_data_filtered,
        cnv_data_filtered,
        table_rna_seq_patients,
        df_melted_protein,
        table_proteins_patients,
    )


def generate_models_re(
    folder_generic_models,
    folder_models,
    top_resistant_ids,
    top_sensitive_ids,
    drug_interest,
    drug_targets,
    phenotype_interest,
    rna_seq_data,
    montagud_nodes,
    table_rna_seq_patients,
    cnv_data_filtered,
    name_maps,
    type_models,
    df_melted_proteins,
    table_proteins_patients,
):
    patients_categ = ["resistant", "sensitive"]
    patients_ids = top_resistant_ids + top_sensitive_ids

    # create generic models cfgs and bnds
    create_generic_patients_cfgs_bnds(
        folder_generic_models,
        folder_models,
        top_resistant_ids,
        top_sensitive_ids,
        drug_interest,
        name_maps,
        type_models,
    )

    models_folder_res = f"{folder_models}/resistant/pers_models"
    models_folder_sens = f"{folder_models}/sensitive/pers_models"

    # simulate drug target
    tailor_bnd_genes_intervention(
        drug_targets, top_resistant_ids, models_folder_res, drug_interest
    )

    tailor_bnd_genes_intervention(
        drug_targets, top_sensitive_ids, models_folder_sens, drug_interest
    )

    for patient_categ in patients_categ:
        folder_models_pheno = f"{folder_models}/{patient_categ}/pers_models"
        generic_models_update_phenotypes(phenotype_interest, folder_models_pheno)

        folder_models_categ = f"{folder_models}/{patient_categ}/pers_models"
        patients_ids_categ = (
            top_sensitive_ids if patient_categ == "sensitive" else top_resistant_ids
        )

        # personalization of the patients models according to the type of models
        if type_models == "genes_models":
            personalized_patients_genes_cfgs(
                rna_seq_data,
                montagud_nodes,
                folder_models_categ,
                patients_ids_categ,
                table_rna_seq_patients,
                drug_interest,
            )
        else:
            personalized_patients_proteins_cfgs(
                df_melted_proteins,
                montagud_nodes,
                folder_models_categ,
                patients_ids,
                table_proteins_patients,
                drug_interest,
            )

        # create personalized models (CNV expression)
        folder_models_cnv = f"{folder_models}/{patient_categ}/pers_models"
        tailor_bnd_cnv_cm(
            cnv_data_filtered, folder_models_cnv, drug_interest=drug_interest
        )
