import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def create_table_patients_phenotypes(folder, dir_res_data, dir_sens_data):

    # Get list of file paths, excluding .DS_Store
    resistant_patients = [
        os.path.join(dir_res_data, f)
        for f in os.listdir(dir_res_data)
        if f != ".DS_Store" and os.path.isfile(os.path.join(dir_res_data, f))
    ]


    sensitive_patients = [
        os.path.join(dir_sens_data, f)
        for f in os.listdir(dir_sens_data)
        if f != ".DS_Store" and os.path.isfile(os.path.join(dir_sens_data, f))
    ]

    all_patients = resistant_patients + sensitive_patients



    # Dictionary to store dataframes for each patient
    patient_dataframes = {}

    # Load all patient dataframes
    for filepath in all_patients:
        patient_name = os.path.splitext(os.path.basename(filepath))[0]
        df = pd.read_csv(filepath, index_col=0)  
        patient_dataframes[patient_name] = df

        # Optional: print confirmation
        print(f"Loaded data for {patient_name}, shape: {df.shape}")

    # Debugging: Check the type of the patient_dataframes
    patient_names = list(patient_dataframes.keys())

  

    # Collect all unique conditions and phenotypes across all patients
    all_conditions = set()
    all_phenotypes = set()
    for df in patient_dataframes.values():
        all_conditions.update(df.index)
        all_phenotypes.update(df.columns)


    # Create an empty DataFrame for storing combined phenotype data, with patient names as index
    patients_phenot_table = pd.DataFrame(index=patient_names)

    # Populate the table with data
    for patient in patient_names:
        patient_df = patient_dataframes[patient]

        for condition in all_conditions:
            for phenotype in all_phenotypes:
                col_name = f"{condition}_{phenotype}"
                if condition in patient_df.index and phenotype in patient_df.columns:
                    value = patient_df.loc[condition, phenotype]
                    patients_phenot_table.at[patient, col_name] = value
                else:
                    print(f"Missing: {col_name} for patient '{patient}'")

                # Add the value in the patients_phenot_table with the patient as the index
                

    output_dir = os.path.join(folder, 'sensitive_resistant_results')
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, 'patients_phenot_table.csv')
    patients_phenot_table.to_csv(output_path, index=True)
    print(f"Saved final phenotype table to: {output_path}")

    return patients_phenot_table



