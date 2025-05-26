import pandas as pd
import os
import re
import numpy as np





def get_patients(
    number_patients,
    drug_data,
    annotations_models,
    drug_interest,
    tissue_interest=None,
    tissue_remove=None,
):
    # extract the top 100 and bottom 100 patients resitant to the drug
    if tissue_interest is not None:
        tissue_interest = tissue_interest.upper()
        annotations_models["tissue"] = annotations_models["tissue"].str.upper()
        annotations_models = annotations_models[
            annotations_models["tissue"] == tissue_interest
        ]

    if tissue_remove is not None:
        tissue_remove = tissue_remove.upper()
        annotations_models["tissue"] = annotations_models["tissue"].str.upper()
        annotations_models = annotations_models[
            annotations_models["tissue"] != tissue_remove
        ]

    annotations_models_filtered = annotations_models[["tissue", "model_id"]]
    models_id = annotations_models["model_id"].tolist()

    # merge drug data with tissue
    annotations_models_filtered.rename(
        columns={"model_id": "SANGER_MODEL_ID"}, inplace=True
    )

    # print(annotations_models_filtered.columns.tolist())

    drug_tissue_data = pd.merge(
        drug_data, annotations_models_filtered, on="SANGER_MODEL_ID"
    )

    drug_data_filtered = drug_data[drug_data["SANGER_MODEL_ID"].isin(models_id)]

    drug_data_filtered = drug_data_filtered[
        ["SANGER_MODEL_ID", "DRUG_NAME", "PUTATIVE_TARGET", "AUC", "LN_IC50", "Z_SCORE"]
    ]
    drug_data_filtered = drug_data_filtered[
        drug_data_filtered["DRUG_NAME"] == drug_interest
    ]
    grouped_drug_data_filtered = (
        drug_data_filtered.groupby("SANGER_MODEL_ID")["Z_SCORE"].mean().reset_index()
    )

    # sort from lowest to highest
    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(
        by="Z_SCORE", ascending=True
    )

    top_sensitive = drug_data_filtered_ranked.iloc[0:number_patients, :]
    # print('sensitive')
    # print(top_sensitive.tail())

    top_sensitive_ids = top_sensitive["SANGER_MODEL_ID"].tolist()

    # sort from highest to lowest
    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(
        by="Z_SCORE", ascending=False
    )
    top_resistant = drug_data_filtered_ranked.iloc[0:number_patients, :]
    # print('resistant')
    # print(top_resistant.tail())

    top_resistant_ids = top_resistant["SANGER_MODEL_ID"].tolist()
    top_resistant_ids = list(set(top_resistant_ids))
    top_sensitive_ids = list(set(top_sensitive_ids))

    return top_resistant_ids, top_sensitive_ids, drug_data_filtered
