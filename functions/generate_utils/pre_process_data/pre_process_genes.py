# Cell model passport - RNA seq expression

import numpy as np
import pandas as pd
import os


def identify_genes_synonyms(data_synonyms, uniprot_data, montagud_nodes):
    data_synonyms = data_synonyms[
        ~(data_synonyms["Gene name"].isna() & data_synonyms["Gene Synonym"].isna())
    ]

    # added - check ??
    data_synonyms["Gene name"] = data_synonyms["Gene name"].astype(str)
    data_synonyms["Gene Synonym"] = data_synonyms["Gene Synonym"].astype(str)


    data_synonyms["Gene name"] = data_synonyms["Gene name"].str.upper()
    data_synonyms["Gene Synonym"] = data_synonyms["Gene Synonym"].str.upper()

    data_synonyms = data_synonyms[
        data_synonyms["Gene name"].isin(montagud_nodes)
        | data_synonyms["Gene Synonym"].isin(montagud_nodes)
    ]

    def match_montagud_node(row, montagud_genes):
        if row["Gene name"] in montagud_genes:
            return row["Gene name"]
        elif row["Gene Synonym"] in montagud_genes:
            return row["Gene Synonym"]
        else:
            return None

    # Apply the function to each row
    data_synonyms["montagud_node"] = data_synonyms.apply(
        lambda row: match_montagud_node(row, montagud_nodes), axis=1
    )

    merged_df = pd.merge(data_synonyms, uniprot_data, on="Gene stable ID")
    merged_df["Gene Synonym"] = merged_df["Gene Synonym"].str.replace(
        r"[_-]", "", regex=True
    )
    merged_df = merged_df[["Gene name_x", "Gene Synonym", "montagud_node"]]

    return merged_df


def classify_expression(z):
    if z > 1:
        return "high"
    elif z < -1:
        return "low"
    else:
        return "normal"


def process_genes(patients_ids, montagud_nodes, rna_seq_data, synonyms_maps):
    rna_seq_data["gene_symbol_upper"] = rna_seq_data["gene_symbol"].str.upper()
    rna_seq_data = rna_seq_data[rna_seq_data["gene_symbol_upper"].isin(montagud_nodes)]
    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]
    rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]

    # remplace the gene symbol column names by its synonyms in the synonyms_maps dictionary
    rna_seq_data["gene_symbol"] = rna_seq_data["gene_symbol"].apply(
        lambda x: synonyms_maps.get(x, x)
    )

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
