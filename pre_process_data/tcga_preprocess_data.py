import numpy as np
import mygene


def pre_process_tcga_data(cnv_data, genes_data, patients_id, montagud_nodes):
    cnv_col = ["Unnamed: 0"] + patients_id
    genes_col = ["xena_sample"] + patients_id
    cnv_data_filtered = cnv_data[cnv_col]
    genes_data_filtered = genes_data[genes_col]

    # filter cnv data
    ensembl_ids_cnv = list(cnv_data_filtered["Unnamed: 0"])
    ensembl_ids_cnv = [eid.split(".")[0] for eid in ensembl_ids_cnv]
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        ensembl_ids_cnv, scopes="ensembl.gene", fields="symbol", species="human"
    )
    cnv_symbols = {item["query"]: item.get("symbol", None) for item in result}
    cnv_names = list(cnv_symbols.values())
    cnv_data_filtered["gene_name"] = cnv_names

    cnv_data_filtered = cnv_data_filtered[cnv_data_filtered["gene_name"].notna()]
    cnv_data_filtered = cnv_data_filtered.drop(columns=["Unnamed: 0"])

    cnv_data_filtered = cnv_data_filtered.set_index("gene_name")
    common_elements = list(set(cnv_data_filtered.index) & set(montagud_nodes))
    cnv_data_filtered = cnv_data_filtered.loc[common_elements]

    cnv_data_filtered = cnv_data_filtered.transpose()

    cnv_data_filtered.index.name = "samples_id"
    cnv_data_filtered = cnv_data_filtered.reset_index()

    df_melted_cnv = cnv_data_filtered.melt(
        id_vars=["samples_id"], var_name="gene_name", value_name="expression_value"
    )

    df_melted_cnv = df_melted_cnv.rename(
        columns={
            "samples_id": "model_id",
            "gene_name": "gene_symbol",
            "expression_value": "rsem_tpm",
        }
    )
    group_loss = [-1]
    group_normal = [0]
    group_gain = [1]

    conditions = [
        df_melted_cnv["rsem_tpm"].isin(group_loss),
        df_melted_cnv["rsem_tpm"].isin(group_normal),
        df_melted_cnv["rsem_tpm"].isin(group_gain),
    ]
    choices = ["Loss", "Normal", "Gain"]
    df_melted_cnv.loc[:, "effect"] = np.select(conditions, choices, default="")
    df_melted_cnv.to_csv("data/TCGA_data/filtered_data/cnv_samples_table.csv")

    # # Filter genes dataset
    ensembl_ids_gene = list(genes_data_filtered["xena_sample"])
    ensembl_ids_gene = [eid.split(".")[0] for eid in ensembl_ids_gene]
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        ensembl_ids_gene, scopes="ensembl.gene", fields="symbol", species="human"
    )
    gene_symbols = {item["query"]: item.get("symbol", None) for item in result}
    genes_names = list(gene_symbols.values())
    genes_data_filtered["gene_name"] = genes_names
    genes_data_filtered = genes_data_filtered[genes_data_filtered["gene_name"].notna()]
    genes_data_filtered = genes_data_filtered.drop(columns=["xena_sample"])

    genes_data_filtered = genes_data_filtered.set_index("gene_name")
    genes_data_filtered = genes_data_filtered.transpose()

    genes_data_filtered = genes_data_filtered.dropna(axis=1)
    genes_data_filtered = genes_data_filtered.dropna(axis=0)

    threshold = 0.5
    zero_proportion = (genes_data_filtered == 0).sum() / len(genes_data_filtered)
    genes_data_filtered_filtered = genes_data_filtered.loc[
        :, zero_proportion <= threshold
    ]

    # keep only the genes in montagud nodes
    genes_data_ids = list(genes_data_filtered_filtered.columns)

    common_elements = list(set(genes_data_ids) & set(montagud_nodes))
    genes_data_filtered = genes_data_filtered_filtered[common_elements]

    genes_data_filtered.index.name = "samples_id"
    genes_data_filtered = genes_data_filtered.reset_index()
    df_melted_genes = genes_data_filtered.melt(
        id_vars=["samples_id"], var_name="gene_name", value_name="expression_value"
    )

    df_melted_genes = df_melted_genes.rename(
        columns={
            "samples_id": "model_id",
            "gene_name": "gene_symbol",
            "expression_value": "rsem_tpm",
        }
    )
    df_melted_genes.to_csv("data/TCGA_data/filtered_data/genes_samples_table.csv")
    return df_melted_cnv, df_melted_genes
