import pandas as pd
import numpy as np


def get_patients_valid(phenotype_data, genes_data, cnv_data):
    phenotype_samples_id = list(phenotype_data["sample"])
    genes_samples_id = list(genes_data.columns)[1:]
    cnv_samples_id = list(cnv_data.columns)[1:]
    common_ids = list(
        set(phenotype_samples_id) & set(genes_samples_id) & set(cnv_samples_id)
    )
    phenotype_data_filtered = phenotype_data[phenotype_data["sample"].isin(common_ids)]
    phenotype_data_filtered = phenotype_data_filtered[
        ["sample", "diagnoses.tumor_stage"]
    ]

    # create groups based on the stages
    group_0 = ["stage 0"]
    group_1 = ["stage i", "stage ia", "stage ib"]
    group_2 = ["stage ii", "stage iia", "stage iib", "stage iic"]
    group_3 = ["stage iii", "stage iiia", "stage iiib", "stage iiic"]
    group_4 = ["stage iv", "stage iva", "stage ivb", "stage ivc"]
    conditions = [
        phenotype_data_filtered["diagnoses.tumor_stage"].isin(group_0),
        phenotype_data_filtered["diagnoses.tumor_stage"].isin(group_1),
        phenotype_data_filtered["diagnoses.tumor_stage"].isin(group_2),
        phenotype_data_filtered["diagnoses.tumor_stage"].isin(group_3),
        phenotype_data_filtered["diagnoses.tumor_stage"].isin(group_4),
    ]
    choices = [
        "Local Cancer",
        "Early Stage",
        "Larger Tumor",
        "Advanced Local Speed",
        "Metastatic",
    ]
    phenotype_data_filtered.loc[:, "Tumor Group"] = np.select(
        conditions, choices, default=""
    )
    sampled_df = phenotype_data_filtered.groupby("Tumor Group", group_keys=False).apply(
        lambda x: x.sample(n=min(len(x), 10), random_state=42)
    )
    sampled_df = sampled_df[
        sampled_df["Tumor Group"].notna() & (sampled_df["Tumor Group"] != "")
    ]
    samples_ids = list(sampled_df["sample"])
    return samples_ids
