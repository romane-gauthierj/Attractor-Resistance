import pandas as pd
import os
import re
import numpy as np


def get_patients(
    number_patients,
    drug_data,
    annotations_models,
    drug_interest,
    tissue_interest,
    tissue_remove,
):
    if tissue_interest is not None:
        tissue_interest = tissue_interest.upper()

        # added - to check
        annotations_models["tissue"] = annotations_models["tissue"].astype(str)
        annotations_models["tissue"] = annotations_models["tissue"].str.upper()
        annotations_models = annotations_models[
            annotations_models["tissue"] == tissue_interest
        ]

    if tissue_remove is not None:
        tissue_remove = tissue_remove.upper()

        annotations_models["tissue"] = annotations_models["tissue"].astype(str)
        annotations_models["tissue"] = annotations_models["tissue"].str.upper()
        annotations_models = annotations_models[
            annotations_models["tissue"] != tissue_remove
        ]

    # control patients
    healthy_ids = list(
        set(
            annotations_models[annotations_models["tissue_status"] == "Normal"][
                "model_id"
            ]
        )
    )

    # remove models with unknown or healthy tissue status

    annotations_models = annotations_models[
        ~annotations_models["tissue_status"].isin(
            ["Unknown", "Normal", "Precancerous", "Papiloma"]
        )
    ]

    annotations_models_filtered = annotations_models[["tissue", "model_id"]]
    models_id = annotations_models["model_id"].tolist()

    # merge drug data with tissue
    annotations_models_filtered.rename(
        columns={"model_id": "SANGER_MODEL_ID"}, inplace=True
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

    # sort from lowest to highest (lowest values -> more sensitive)
    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(
        by="Z_SCORE", ascending=True
    )

    top_sensitive = drug_data_filtered_ranked.iloc[0:number_patients, :]
    # print('sensitive')
    # print(top_sensitive.tail())

    top_sensitive_ids = top_sensitive["SANGER_MODEL_ID"].tolist()

    # sort from highest to lowest (highest values -> more resistant)
    drug_data_filtered_ranked = grouped_drug_data_filtered.sort_values(
        by="Z_SCORE", ascending=False
    )
    top_resistant = drug_data_filtered_ranked.iloc[0:number_patients, :]
    # print('resistant')
    # print(top_resistant.tail())

    top_resistant_ids = top_resistant["SANGER_MODEL_ID"].tolist()
    top_resistant_ids = list(set(top_resistant_ids))
    top_sensitive_ids = list(set(top_sensitive_ids))

    return top_resistant_ids, top_sensitive_ids, healthy_ids
