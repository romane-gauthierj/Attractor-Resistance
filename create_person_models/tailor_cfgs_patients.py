# pipeline to tailor the generic lung adenocarcinoma boolean network to patient specific boolean

import pandas as pd
import re
import os


# function used for validation and general pipeline
def personalized_patients_genes_cfgs(
    rna_seq_data,
    montagud_nodes,
    original_data_dir,
    results_dir,
    patients_ids,
    rna_seq_data_filtered,
    drug_name,
):
    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]

    rna_seq_data = rna_seq_data[["model_id", "gene_symbol", "rsem_tpm"]]

    # Open each file in the directory
    os.makedirs(original_data_dir, exist_ok=True)

    # Directory where modified files will be saved
    os.makedirs(results_dir, exist_ok=True)

    rna_seq_data_max = rna_seq_data
    rna_seq_data_max["gene_max"] = rna_seq_data_max.groupby("gene_symbol")[
        "rsem_tpm"
    ].transform("max")

    # Loop through each file in the directory
    for filename in os.listdir(original_data_dir):
        file_path = os.path.join(original_data_dir, filename)
        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0].replace(
                f"_{drug_name}", ""
            )
            # Only proceed if model_id is in table_genes_patients
            if model_id_in_file in rna_seq_data_filtered.index:
                with open(file_path, "r") as file:
                    content = file.read()
                    # rna_seq_data_filtered_filtered = rna_seq_data_filtered[rna_seq_data_filtered.index == model_id_in_file]
                if model_id_in_file in rna_seq_data_filtered.index:
                    row = rna_seq_data_filtered.loc[model_id_in_file]
                    # Only modify if this row corresponds to the file being processed
                    high_genes = row["High Gene Expression"]
                    low_genes = row["Low Gene Expression"]
                    gene_high_list = [
                        gene.strip() for gene in high_genes.split(",") if gene.strip()
                    ]
                    gene_low_list = [
                        gene.strip() for gene in low_genes.split(",") if gene.strip()
                    ]
                    gene_list = list(set(gene_high_list + gene_low_list))
                    for gene in gene_list:  # add condition if gene is in the genes of the boolean generic network
                        if gene in montagud_nodes:
                            expr_max = rna_seq_data_max.loc[
                                rna_seq_data_max["gene_symbol"] == gene, "gene_max"
                            ].iloc[0]
                            gene_expr = rna_seq_data.loc[
                                (rna_seq_data["gene_symbol"] == gene)
                                & (rna_seq_data["model_id"] == model_id_in_file),
                                "rsem_tpm",
                            ].iloc[0]
                            prob_1 = min(max(gene_expr / expr_max, 0), 1)

                            u_pattern = re.compile(
                                rf"\$u_{gene}\s*=\s*(0|1);", re.DOTALL
                            )
                            d_pattern = re.compile(
                                rf"\$d_{gene}\s*=\s*(0|1);", re.DOTALL
                            )
                            u_line = f"$u_{gene} = {prob_1:.4f};"
                            d_line = f"$d_{gene} = {1 - prob_1:.4f};"
                            content = re.sub(u_pattern, u_line, content)
                            content = re.sub(d_pattern, d_line, content)
                # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
                modified_file_path = os.path.join(results_dir, filename)

                with open(modified_file_path, "w") as file:
                    file.write(content)
                print(f"Modified and saved: {modified_file_path}")


# function used for validation and general pipeline
def personalized_patients_proteins_cfgs(
    df_melted_proteins,
    montagud_nodes,
    original_data_dir,
    results_dir,
    patients_ids,
    table_proteins_patients,
    drug_name,
):
    df_melted_proteins = df_melted_proteins[
        df_melted_proteins["model_id"].isin(patients_ids)
    ]

    df_melted_proteins = df_melted_proteins[["model_id", "protein_symbol", "rsem_tpm"]]

    # Open each file in the directory
    os.makedirs(original_data_dir, exist_ok=True)

    # Directory where modified files will be saved
    os.makedirs(results_dir, exist_ok=True)

    protein_data_max = df_melted_proteins
    protein_data_max["protein_max"] = protein_data_max.groupby("protein_symbol")[
        "rsem_tpm"
    ].transform("max")

    # Loop through each file in the directory
    for filename in os.listdir(original_data_dir):
        if not filename.endswith(".cfg"):
            continue
        file_path = os.path.join(original_data_dir, filename)
        if os.path.isfile(file_path):  # Make sure it's a file, not a subdirectory
            # model_id_in_file = os.path.splitext(filename)[0]  # remove .cfg extension
            # model_id_in_file = os.path.splitext(filename)[0].replace('_AZD8931', '')
            model_id_in_file = os.path.splitext(filename)[0].replace(
                f"_{drug_name}", ""
            )
            # Only proceed if model_id is in table_proteins_patients
            if model_id_in_file in table_proteins_patients.index:
                with open(file_path, "r") as file:
                    content = file.read()
                    # rna_seq_data_filtered_filtered = rna_seq_data_filtered[rna_seq_data_filtered.index == model_id_in_file]
                if model_id_in_file in table_proteins_patients.index:
                    row = table_proteins_patients.loc[model_id_in_file]
                    # Only modify if this row corresponds to the file being processed
                    high_proteins = row["High Protein Abundance"]
                    low_proteins = row["Low Protein Abundance"]
                    protein_high_list = [
                        protein.strip()
                        for protein in high_proteins.split(",")
                        if protein.strip()
                    ]
                    protein_low_list = [
                        protein.strip()
                        for protein in low_proteins.split(",")
                        if protein.strip()
                    ]
                    protein_list = list(set(protein_high_list + protein_low_list))
                    for protein in protein_list:  # add condition if protein is in the proteins of the boolean generic network
                        if protein in montagud_nodes:
                            expr_max = protein_data_max.loc[
                                protein_data_max["protein_symbol"] == protein,
                                "protein_max",
                            ].iloc[0]

                            prot_expr = df_melted_proteins.loc[
                                (df_melted_proteins["protein_symbol"] == protein)
                                & (df_melted_proteins["model_id"] == model_id_in_file),
                                "rsem_tpm",
                            ].iloc[0]

                            prob_1 = min(max(prot_expr / expr_max, 0), 1)

                            u_pattern = re.compile(
                                rf"\$u_{protein}\s*=\s*(0|1);", re.DOTALL
                            )

                            d_pattern = re.compile(
                                rf"\$d_{protein}\s*=\s*(0|1);", re.DOTALL
                            )
                            u_line = f"$u_{protein} = {prob_1:.4f};"
                            d_line = f"$d_{protein} = {1 - prob_1:.4f};"
                            content = re.sub(u_pattern, u_line, content)
                            content = re.sub(d_pattern, d_line, content)

                            if protein == "CASPASE3":
                                print("Processing CASPASE3...")

                                u_pattern = re.compile(
                                    rf"\$u_{protein}\s*=\s*\d+\.?\d*;", re.DOTALL
                                )
                                d_pattern = re.compile(
                                    rf"\$d_{protein}\s*=\s*\d+\.?\d*;", re.DOTALL
                                )

                                u_line = f"$u_{protein} = {prob_1:.4f};"
                                d_line = f"$d_{protein} = {1 - prob_1:.4f};"

                                content, num_u = re.subn(u_pattern, u_line, content)
                                content, num_d = re.subn(d_pattern, d_line, content)

                                print(f"Replaced $u: {num_u} times, $d: {num_d} times")
                                print("Final content snippet:")
                                print(
                                    "\n".join(
                                        [
                                            line
                                            for line in content.splitlines()
                                            if "CASPASE3" in line
                                        ]
                                    )
                                )

                # modified_file_path = os.path.join(modified_output_dir, f'{filename}_{drug_name}')
                modified_file_path = os.path.join(results_dir, filename)

                with open(modified_file_path, "w") as file:
                    file.write(content)
                print(f"Modified and saved: {modified_file_path}")
