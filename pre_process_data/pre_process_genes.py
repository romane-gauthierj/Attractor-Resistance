# Cell model passport - RNA seq expression

import numpy as np
import pandas as pd
import os


def classify_expression(z):
    if z > 1:
        return "high"
    elif z < -1:
        return "low"
    else:
        return "normal"


def process_genes(patients_ids, montagud_nodes, rna_seq_data):
    rna_seq_data["gene_symbol_upper"] = rna_seq_data["gene_symbol"].str.upper()
    rna_seq_data = rna_seq_data[rna_seq_data["gene_symbol_upper"].isin(montagud_nodes)]
    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]
    rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]
    return rna_seq_data


def create_table_rna_seq_patients(rna_seq_data):
    # filter only data with patient id and the nodes of the montagud_data

    rna_seq_data["z_score"] = rna_seq_data.groupby("gene_symbol")["rsem_tpm"].transform(
        lambda x: (x - x.mean()) / x.std()
    )

    rna_seq_data["gene_expression_level"] = rna_seq_data["z_score"].apply(
        classify_expression
    )

    rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "gene_expression_level"]]
    rna_seq_data.rename(columns={"gene_symbol": "gene_name"}, inplace=True)
    rna_seq_data = rna_seq_data[
        rna_seq_data["gene_expression_level"].isin(["low", "high"])
    ]

    table_rna_seq_patients = rna_seq_data.pivot_table(
        index="model_id",
        columns="gene_expression_level",
        values="gene_name",
        aggfunc=lambda x: ", ".join(sorted(set(x))),
    ).fillna("-")

    table_rna_seq_patients = table_rna_seq_patients.rename(
        columns={"low": "Low Gene Expression", "high": "High Gene Expression"}
    )

    return table_rna_seq_patients


def create_table_proteins_patients(proteins_data):
    # filter only data with patient id and the nodes of the montagud_data

    proteins_data["z_score"] = proteins_data.groupby("protein_symbol")[
        "rsem_tpm"
    ].transform(lambda x: (x - x.mean()) / x.std())

    proteins_data["protein_expression_level"] = proteins_data["z_score"].apply(
        classify_expression
    )

    proteins_data = proteins_data[
        ["model_id", "protein_symbol", "protein_expression_level"]
    ]
    proteins_data.rename(columns={"protein_symbol": "protein_name"}, inplace=True)
    proteins_data = proteins_data[
        proteins_data["protein_expression_level"].isin(["low", "high"])
    ]

    table_proteins_patients = proteins_data.pivot_table(
        index="model_id",
        columns="protein_expression_level",
        values="protein_name",
        aggfunc=lambda x: ", ".join(sorted(set(x))),
    ).fillna("-")

    table_proteins_patients = table_proteins_patients.rename(
        columns={"low": "Low Protein Abundance", "high": "High Protein Abundance"}
    )

    return table_proteins_patients
