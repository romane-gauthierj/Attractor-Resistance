import pandas as pd


def get_cell_lines(patients_ids, patient_categ, annotations_models, folder_result):
    annotations_models_filtered = annotations_models[
        annotations_models["model_id"].isin(patients_ids)
    ][["model_id", "sample_site", "tissue"]]
    annotations_models_filtered.to_csv(
        f"{folder_result}/filtered_cell_lines_{patient_categ}.csv", index=False
    )
