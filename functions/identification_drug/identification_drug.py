import numpy as np
import pandas as pd
from identification_patients.get_patients_sens_res import get_patients

def identify_drug(drug_data, annotations_models, tissue_remove):
    results = {}
    drug_interests = drug_data["DRUG_NAME"].unique().tolist()

    for drug in drug_interests:
        top_resistant_ids, top_sensitive_ids, drug_data_filtered = get_patients(500, drug_data, annotations_models, drug, tissue_remove)
        patients_ids = top_resistant_ids + top_sensitive_ids
        patients_ids = top_resistant_ids + top_sensitive_ids

        drug_tissue_filtered = drug_data_filtered[
            drug_data_filtered["SANGER_MODEL_ID"].isin(patients_ids)
        ]
        results[drug] = {
            "name": drug,
            "<-1.5": float((drug_tissue_filtered["Z_SCORE"] < -1.5).sum()),
            ">1.5": float((drug_tissue_filtered["Z_SCORE"] > 1.5).sum()),
            "mean": float(drug_tissue_filtered["Z_SCORE"].mean()),
            "std": float(drug_tissue_filtered["Z_SCORE"].std()),
            "abs_zscore": float(drug_tissue_filtered["Z_SCORE"].abs().mean()),
        }

    pd_results = pd.DataFrame(results.values())
    pd_results.to_csv("drug_analysis.csv", index=False)
    return pd_results
