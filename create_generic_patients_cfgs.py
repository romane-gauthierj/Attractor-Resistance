import pandas as pd
import os
import re

# Import the function to get the patients ID
from get_patients_sens_res import get_patients_top_10

# Load patient IDs

def create_generic_patients_cfgs_bnds(folder, resistant_patients_ids, sensitive_patients_ids, drug_interest):
    # --- Templates ---
    cfg_template_path = 'generic_lung_model/Adeno_lung_Cancer.cfg'
    bnd_template_path = 'generic_lung_model/Adeno_lung_Cancer.bnd'  

    # --- Sensitive Patients ---
    sensitive_cfg_output_dir = f'{folder}/sensitive_patient/generic_models'
    sensitive_bnd_output_dir = f'{folder}/sensitive_patient/generic_models'
    os.makedirs(sensitive_cfg_output_dir, exist_ok=True)
    os.makedirs(sensitive_bnd_output_dir, exist_ok=True)

    # Read templates once
    with open(cfg_template_path, 'r') as file:
        cfg_template_content = file.read()

    with open(bnd_template_path, 'r') as file:
        bnd_template_content = file.read()

    # Create personalized .cfg and .bnd for sensitive patients
    for patient_id in sensitive_patients_ids:
        patient_cfg = cfg_template_content.replace('{PATIENT_ID}', patient_id)
        patient_bnd = bnd_template_content.replace('{PATIENT_ID}', patient_id)

        cfg_output_path = os.path.join(sensitive_cfg_output_dir, f'{patient_id}_{drug_interest}.cfg')
        bnd_output_path = os.path.join(sensitive_bnd_output_dir, f'{patient_id}_{drug_interest}.bnd')

        with open(cfg_output_path, 'w') as file:
            file.write(patient_cfg)
        with open(bnd_output_path, 'w') as file:
            file.write(patient_bnd)

    # --- Resistant Patients ---
    resistant_cfg_output_dir = f'{folder}/resistant_patient/generic_models'
    resistant_bnd_output_dir = f'{folder}/resistant_patient/generic_models'
    os.makedirs(resistant_cfg_output_dir, exist_ok=True)
    os.makedirs(resistant_bnd_output_dir, exist_ok=True)

    # Read templates again if needed (same ones reused here)
    # You can skip this if using same content as above

    # Create personalized .cfg and .bnd for resistant patients
    for patient_id in resistant_patients_ids:
        patient_cfg = cfg_template_content.replace('{PATIENT_ID}', patient_id)
        patient_bnd = bnd_template_content.replace('{PATIENT_ID}', patient_id)

        cfg_output_path = os.path.join(resistant_cfg_output_dir, f'{patient_id}_{drug_interest}.cfg')
        bnd_output_path = os.path.join(resistant_bnd_output_dir, f'{patient_id}_{drug_interest}.bnd')

        with open(cfg_output_path, 'w') as file:
            file.write(patient_cfg)
        with open(bnd_output_path, 'w') as file:
            file.write(patient_bnd)

    print("All .cfg and .bnd files created for sensitive and resistant patients.")
