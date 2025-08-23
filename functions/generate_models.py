import os
import logging

logger = logging.getLogger(__name__)

from functions.generate_utils.identification_patients.get_patients_sens_res import (
    get_patients,
)
from functions.generate_utils.pre_process_data.pre_process_cnv import preprocess_cnv
from functions.generate_utils.pre_process_data.pre_process_genes import (
    process_genes,process_genes_proteins
)
from functions.generate_utils.pre_process_data.pre_process_montagud_nodes import (
    process_montagud_nodes,process_montagud_nodes_synonyms
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

from functions.generate_utils.create_person_models.tailor_bnd_cnv import (
    tailor_bnd_cnv_cm,
)

from functions.analysis_utils.genes_intervention.pers_interventions import (
    tailor_bnd_genes_intervention,
)

from functions.generate_utils.pre_process_data.pre_process_proteins import (
    process_proteins,
)
from functions.generate_utils.pre_process_data.preprocess_mutations import pre_process_mutations


def pre_process_re(
    montagud_original_data_df,
    nodes_montagud_synonyms,
    rna_seq_data,
    cnv_data,
    mutations_data,
    number_patients,
    drug_data,
    annotations_models,
    drug_interest,
    proteins_data,
    type_models,
    onco_tsg_data,
    tissue_interest=None,
    tissue_remove=None,
):
    top_resistant_ids, top_sensitive_ids, healthy_ids = get_patients(
        number_patients,
        drug_data,
        annotations_models,
        drug_interest,
        tissue_interest=tissue_interest,
        tissue_remove=tissue_remove,
    )

    patients_ids = top_resistant_ids + top_sensitive_ids + healthy_ids

  

    montagud_node_synonyms, synonyms_to_nodes_dict = process_montagud_nodes_synonyms(nodes_montagud_synonyms)


    montagud_node_model, all_montagud_nodes = process_montagud_nodes(
        montagud_original_data_df, montagud_node_synonyms
    )

 
    rna_seq_data_models_filtered = process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict)
    

    # pre process proteins data
    df_melted_protein = process_proteins(patients_ids, proteins_data, all_montagud_nodes, synonyms_to_nodes_dict)



 #TODO - check this is correct
    mutations_data_filtered = pre_process_mutations(mutations_data, patients_ids, all_montagud_nodes, onco_tsg_data,synonyms_to_nodes_dict)




    if type_models == "genes_models":
        top_resistant_ids = list(
            set(rna_seq_data_models_filtered["model_id"]) & set(top_resistant_ids)
        )
        top_sensitive_ids = list(
            set(rna_seq_data_models_filtered["model_id"]) & set(top_sensitive_ids)
        )
        top_healthy_ids = list(set(rna_seq_data_models_filtered["model_id"]) & set(healthy_ids))

    elif type_models == "proteins_models":
        # if proteins models
        top_resistant_ids = list(
            set(df_melted_protein["model_id"]) & set(top_resistant_ids)
        )
        top_sensitive_ids = list(
            set(df_melted_protein["model_id"]) & set(top_sensitive_ids)
        )
        top_healthy_ids = list(set(df_melted_protein["model_id"]) & set(healthy_ids))

    elif type_models == "genes_proteins_models":

        # priority to proteins  -> only use genes data if not in protein
        rna_seq_data_models_filtered = process_genes_proteins(df_melted_protein,rna_seq_data_models_filtered)

        top_resistant_ids = list(
            set(df_melted_protein["model_id"])
            & set(rna_seq_data_models_filtered["model_id"])
            & set(top_resistant_ids)
        )
        top_sensitive_ids = list(
            set(df_melted_protein["model_id"])
            & set(rna_seq_data_models_filtered["model_id"])
            & set(top_sensitive_ids)
        )
        top_healthy_ids = list(
            set(df_melted_protein["model_id"])
            & set(rna_seq_data_models_filtered["model_id"])
            & set(healthy_ids)
        )
    else:
        raise ValueError(
            "Type of models not recognized. Please choose from 'genes_models', 'proteins_models', or 'genes_proteins_models'."
        )

    patients_ids = top_resistant_ids + top_sensitive_ids + top_healthy_ids

    os.makedirs(f"analysis/{drug_interest}", exist_ok=True)

    # Save IDs to txt files
    with open("analysis/{}/top_resistant_ids.txt".format(drug_interest), "w") as f:
        for pid in top_resistant_ids:
            f.write(f"{pid}\n")
    with open("analysis/{}/top_sensitive_ids.txt".format(drug_interest), "w") as f:
        for pid in top_sensitive_ids:
            f.write(f"{pid}\n")
    with open("analysis/{}/top_healthy_ids.txt".format(drug_interest), "w") as f:
        for pid in top_healthy_ids:
            f.write(f"{pid}\n")

    # preprocess cnv data
    cnv_data_filtered = preprocess_cnv(patients_ids,
        cnv_data, all_montagud_nodes, synonyms_to_nodes_dict
    )


    return (
        top_resistant_ids,
        top_sensitive_ids,
        top_healthy_ids,
        montagud_node_model,
        all_montagud_nodes,
        rna_seq_data_models_filtered,
        cnv_data_filtered,
        df_melted_protein,
        mutations_data_filtered,
    )


def generate_models_re(
    normalization_method,
    discrete_variable,
    continuous_variable,
    folder_generic_models,
    folder_models,
    top_resistant_ids,
    top_sensitive_ids,
    top_healthy_ids,
    drug_interest,
    drug_targets,
    phenotype_interest,
    rna_seq_data_models_filtered,
    montagud_node_model,
    cnv_data_filtered,
    mutations_data_filtered,
    type_models,
    df_melted_proteins,
    amplif_factor,
    intervention_gene=None,
    genetic_intervention=None,
):
    patients_categ = ["resistant", "sensitive", "healthy"]
    # patients_ids = top_resistant_ids + top_sensitive_ids + top_healthy_ids

    # create generic models cfgs and bnds

    create_generic_patients_cfgs_bnds(
        folder_generic_models,
        folder_models,
        top_resistant_ids,
        top_sensitive_ids,
        top_healthy_ids,
        drug_interest,
    )


    models_folder_res = f"{folder_models}/resistant/pers_models"
    models_folder_sens = f"{folder_models}/sensitive/pers_models"
    models_folder_healthy = f"{folder_models}/healthy/pers_models"




    for patient_categ in patients_categ:
        folder_models_pheno = f"{folder_models}/{patient_categ}/pers_models"
        generic_models_update_phenotypes(phenotype_interest, folder_models_pheno)

        folder_models_categ = f"{folder_models}/{patient_categ}/pers_models"

        if patient_categ == "sensitive":
            patients_ids_categ = top_sensitive_ids
        elif patient_categ == "resistant":
            patients_ids_categ = top_resistant_ids
        else:  # "healthy"
            patients_ids_categ = top_healthy_ids

        # personalization of the patients models according to the type of models
        if continuous_variable == "genes":
            personalized_patients_genes_cfgs(
                rna_seq_data_models_filtered,
                montagud_node_model,
                folder_models_categ,
                amplif_factor,
                drug_interest, # context label
                normalization_method,
            )


# test
        elif continuous_variable == "proteins":
            personalized_patients_proteins_cfgs(
                df_melted_proteins,
                montagud_node_model,
                folder_models_categ,
                amplif_factor,
                drug_interest,  # context label
                normalization_method,
            )


        elif continuous_variable == "genes_proteins":
            # proteins will overwrite genes
            personalized_patients_genes_cfgs(
                rna_seq_data_models_filtered,
                montagud_node_model,
                folder_models_categ,
                amplif_factor,
                drug_interest, # context label
                normalization_method,
            )
            # test
            personalized_patients_proteins_cfgs(
                df_melted_proteins,
                montagud_node_model,
                folder_models_categ,
                amplif_factor,
                drug_interest,  # context label
                normalization_method,
            )

    
        else:
            raise ValueError(
                "Type of models not recognized. Please choose from 'genes_models', 'proteins_models', or 'genes_proteins_models'."
            )

        # create personalized models (CNV expression)
        folder_models_cnv = f"{folder_models}/{patient_categ}/pers_models"


        if discrete_variable == 'mutations':
            logger.debug('mutations is used to personalize the networks')
            tailor_bnd_cnv_cm(
                mutations_data_filtered, folder_models_cnv,
            )
        elif discrete_variable == 'cnv':
            logger.debug('cnv is used to personalize the networks')
            tailor_bnd_cnv_cm(
                cnv_data_filtered, folder_models_cnv,
            )
        elif discrete_variable == 'mutations_cnv':
            tailor_bnd_cnv_cm(
                cnv_data_filtered, folder_models_cnv,
            )
            tailor_bnd_cnv_cm(
                mutations_data_filtered, folder_models_cnv,
            )


            logger.debug('mutations and cnv are used to personalize the networks')
        
        else:
            logger.debug('Please select a discrete variable between cnv, mutations or mutations_cnv')

     # simulate drug target (overwrite all)
    tailor_bnd_genes_intervention(
        drug_targets,
        genetic_intervention,
        top_resistant_ids,
        models_folder_res,
        drug_interest,
    )

    tailor_bnd_genes_intervention(
        drug_targets,
        genetic_intervention,
        top_sensitive_ids,
        models_folder_sens,
        drug_interest,
    )
    tailor_bnd_genes_intervention(
        drug_targets,
        genetic_intervention,
        top_healthy_ids,
        models_folder_healthy,
        drug_interest,
    )


    # simulate gene intervention if provided
    if intervention_gene is not None and genetic_intervention is not None:
        tailor_bnd_genes_intervention(
            intervention_gene,
            genetic_intervention,
            top_resistant_ids,
            models_folder_res,
            drug_interest,
        )
        tailor_bnd_genes_intervention(
            intervention_gene,
            genetic_intervention,
            top_sensitive_ids,
            models_folder_sens,
            drug_interest,
        )
        tailor_bnd_genes_intervention(
            intervention_gene,
            genetic_intervention,
            top_healthy_ids,
            models_folder_healthy,
            drug_interest,
        )
