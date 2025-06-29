import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA



def create_combined_table_patients(res_tables_path, sens_tables_path, folder_results):
    # Combines individual patient phenotype result CSVs (from resistant and sensitive folders) into a single DataFrame.

    files_resistant = [
        f
        for f in os.listdir(res_tables_path)
        if os.path.isfile(os.path.join(res_tables_path, f))
    ]
    files_sensitive = [
        f
        for f in os.listdir(sens_tables_path)
        if os.path.isfile(os.path.join(sens_tables_path, f))
    ]
    # print(files_sensitive)

    patient_dataframes = {}

    for path in [res_tables_path, sens_tables_path]:
        if os.path.exists(path):
            for filename in os.listdir(path):
                file_path = os.path.join(path, filename)
                if (
                    os.path.isfile(file_path)
                    and filename.endswith(".csv")
                    and filename.startswith("SIDM")
                ):
                    patient_id = os.path.splitext(filename)[0].replace("_", "")
                    df = pd.read_csv(file_path, index_col=0)
                    patient_dataframes[patient_id] = df
    patient_names = list(patient_dataframes.keys())
    print(patient_names)

    # Collect all unique conditions and phenotypes across all patients
    all_conditions = set()
    all_phenotypes = set()
    for df in patient_dataframes.values():
        all_conditions.update(df.index)
        all_phenotypes.update(df.columns)

    # Create an empty DataFrame for storing combined phenotype data, with patient names as index
    patients_phenot_table = pd.DataFrame(index=patient_names)
    # patient_dataframes[patient_name] = df

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

    output_dir = os.path.join(folder_results, "sensitive_resistant_results")
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "patients_phenot_table.csv")
    patients_phenot_table.to_csv(output_path, index=True)
    print(f"Saved final phenotype table to: {output_path}")
    return patients_phenot_table
