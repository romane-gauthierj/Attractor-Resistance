import pandas as pd
import os
import re


from .update_nodes_names import replace_node_names_in_file


def create_generic_patients_cfg_bnd_validation(
    cfg_template_path,
    bnd_template_path,
    folder_pers_models,
    patients_ids,
    patients_groups,
    tissue,
    name_maps,
):
    """
    Create personalized .cfg and .bnd files for a list of patient IDs based on generic templates.

    Parameters:
    - cfg_template_path: str, path to the generic .cfg template file
    - bnd_template_path: str, path to the generic .bnd template file
    - folder_pers_models: str, output folder where personalized files will be saved
    - patients_ids: list of str, patient identifiers
    - tissue: str, tissue type to include in the output filenames
    - type_model: gene/ proteins model
    """

    # --- Pre process the generic model ---
    replace_node_names_in_file(cfg_template_path, name_maps)
    replace_node_names_in_file(bnd_template_path, name_maps)

    # Ensure output directory exists
    os.makedirs(folder_pers_models, exist_ok=True)

    # Read templates once
    with open(cfg_template_path, "r") as file:
        cfg_template_content = file.read()

    with open(bnd_template_path, "r") as file:
        bnd_template_content = file.read()

    # Create personalized .cfg and .bnd for each patient and saved them in their group folders
    for patient_id in patients_ids:
        group = patients_groups[patients_groups["sampleID"] == patient_id][
            "Gleason_group"
        ].iloc[0]
        group_folder = os.path.join(folder_pers_models, group)
        os.makedirs(group_folder, exist_ok=True)

        patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
        patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

        cfg_output_path = os.path.join(group_folder, f"{patient_id}_{tissue}.cfg")
        bnd_output_path = os.path.join(group_folder, f"{patient_id}_{tissue}.bnd")

        with open(cfg_output_path, "w") as file:
            file.write(patient_cfg)
        with open(bnd_output_path, "w") as file:
            file.write(patient_bnd)

    print("All .cfg and .bnd files created for the validation.")


# Load patient IDs


def create_generic_patients_cfgs_bnds(
    folder_generic_models,
    folder_models,
    resistant_patients_ids,
    sensitive_patients_ids,
    top_healthy_ids,
    drug_interest,
    name_maps,
    type_models,
):
    # --- Templates ---
    cfg_template_path = (
        folder_generic_models + "Montagud2022_Prostate_Cancer_original.cfg"
    )
    bnd_template_path = (
        folder_generic_models + "Montagud2022_Prostate_Cancer_original.bnd"
    )

    # --- Pre process the generic model (proteins or genes names) ---
    replace_node_names_in_file(cfg_template_path, name_maps)
    replace_node_names_in_file(bnd_template_path, name_maps)

    # --- Sensitive Patients ---

    patients_categs = ["sensitive", "resistant", "healthy"]
    for patient_categ in patients_categs:
        cfg_output_dir = f"{folder_models}/{patient_categ}/pers_models"
        bnd_output_dir = f"{folder_models}/{patient_categ}/pers_models"
        os.makedirs(cfg_output_dir, exist_ok=True)
        os.makedirs(bnd_output_dir, exist_ok=True)

        # Read templates once
        with open(cfg_template_path, "r") as file:
            cfg_template_content = file.read()
        with open(bnd_template_path, "r") as file:
            bnd_template_content = file.read()

        if patient_categ == "sensitive":
            patients_ids = sensitive_patients_ids
        elif patient_categ == "healthy":
            patients_ids = top_healthy_ids
        else:  # resistant patients
            patients_ids = resistant_patients_ids

        for patient_id in patients_ids:
            patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
            patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

            cfg_output_path = os.path.join(
                cfg_output_dir, f"{patient_id}_{drug_interest}.cfg"
            )
            bnd_output_path = os.path.join(
                bnd_output_dir, f"{patient_id}_{drug_interest}.bnd"
            )

            with open(cfg_output_path, "w") as file:
                file.write(patient_cfg)
            with open(bnd_output_path, "w") as file:
                file.write(patient_bnd)
    print(
        "All .cfg and .bnd files created for sensitive, resistant and healthy patients."
    )
