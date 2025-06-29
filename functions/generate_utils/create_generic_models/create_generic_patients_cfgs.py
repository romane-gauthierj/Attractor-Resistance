import pandas as pd
import os
import re


import os


def create_generic_patients_cfg_bnd_validation(
    cfg_template_path, bnd_template_path, folder_pers_models, patients_ids, tissue
):
    """
    Create personalized .cfg and .bnd files for a list of patient IDs based on generic templates.

    Parameters:
    - cfg_template_path: str, path to the generic .cfg template file
    - bnd_template_path: str, path to the generic .bnd template file
    - folder_pers_models: str, output folder where personalized files will be saved
    - patients_ids: list of str, patient identifiers
    - tissue: str, tissue type to include in the output filenames
    """
    # Ensure output directory exists
    os.makedirs(folder_pers_models, exist_ok=True)

    # Read templates once
    with open(cfg_template_path, "r") as file:
        cfg_template_content = file.read()

    with open(bnd_template_path, "r") as file:
        bnd_template_content = file.read()

    # Create personalized .cfg and .bnd for each patient
    for patient_id in patients_ids:
        patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
        patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

        cfg_output_path = os.path.join(folder_pers_models, f"{patient_id}_{tissue}.cfg")
        bnd_output_path = os.path.join(folder_pers_models, f"{patient_id}_{tissue}.bnd")

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
    drug_interest,
):
    # --- Templates ---
    cfg_template_path = folder_generic_models + "Montagud2022_Prostate_Cancer.cfg"
    bnd_template_path = folder_generic_models + "Montagud2022_Prostate_Cancer.bnd"

    # --- Sensitive Patients ---

    patients_categs = ["sensitive", "resistant"]
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
        else:
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
    print("All .cfg and .bnd files created for sensitive and resistant patients.")

    # # --- Resistant Patients ---
    # resistant_cfg_output_dir = f"{folder_pers_models}/resistant_patient/generic_models"
    # resistant_bnd_output_dir = f"{folder_pers_models}/resistant_patient/generic_models"
    # os.makedirs(resistant_cfg_output_dir, exist_ok=True)
    # os.makedirs(resistant_bnd_output_dir, exist_ok=True)

    # # Read templates again if needed (same ones reused here)
    # # You can skip this if using same content as above

    # # Create personalized .cfg and .bnd for resistant patients
    # for patient_id in resistant_patients_ids:
    #     patient_cfg = cfg_template_content.replace("{PATIENT_ID}", patient_id)
    #     patient_bnd = bnd_template_content.replace("{PATIENT_ID}", patient_id)

    #     cfg_output_path = os.path.join(
    #         resistant_cfg_output_dir, f"{patient_id}_{drug_interest}.cfg"
    #     )
    #     bnd_output_path = os.path.join(
    #         resistant_bnd_output_dir, f"{patient_id}_{drug_interest}.bnd"
    #     )

    #     with open(cfg_output_path, "w") as file:
    #         file.write(patient_cfg)
    #     with open(bnd_output_path, "w") as file:
    #         file.write(patient_bnd)
