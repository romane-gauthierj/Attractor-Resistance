import pandas as pd
import os
import re


def tailor_bnd_cnv_validation(df_melted_cnv, bnd_dir, tissue):
    gain_group = df_melted_cnv[df_melted_cnv["effect"] == "Gain"]
    loss_group = df_melted_cnv[df_melted_cnv["effect"] == "Loss"]

    gain_ids = set(gain_group["model_id"])
    loss_ids = set(loss_group["model_id"])
    patients_all = list(gain_ids | loss_ids)

    for patient in patients_all:
        bnd_file = [
            file
            for file in os.listdir(bnd_dir)
            if file.endswith(".bnd") and patient in file
        ]
        if not bnd_file:
            print(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(bnd_dir, bnd_file[0])
        genes = df_melted_cnv[df_melted_cnv["model_id"] == patient]["gene_symbol"]
        gene_list = genes.tolist()

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0].replace(
            f"_{tissue}", ""
        )
        with open(bnd_file, "r") as file:
            content = file.read()

        modified_any = False
        for gene in gene_list:
            print(f"🔍 Processing patient {patient}, gene: {gene}")
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            if patient in gain_ids and patient in loss_ids:
                print(
                    f"Patient {patient} is in both gain and loss groups. Please review."
                )
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in gain_ids:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in loss_ids:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                print(f"Patient {patient} not found in gain or loss CNV groups.")
                continue

            gene_match = pattern.search(content)
            if gene_match:
                print(f"{gene} node found. Replacing...")
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    modified_any = True
            else:
                print(f"No {gene} node found in file for patient {patient_id}")

        if modified_any:
            print(f"{patient_id}: CNV — nodes modified")
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            print(f"{patient_id}: CNV — no nodes modified")


def tailor_bnd_cnv_cm(cnv_data_filtered, bnd_dir_res, bnd_dir_sens, drug_interest):
    gain_group = cnv_data_filtered[
        cnv_data_filtered["cn_category"].isin(["Amplification", "Gain"])
        & (cnv_data_filtered["total_copy_number"] > 2.0)
    ]

    loss_group = cnv_data_filtered[
        cnv_data_filtered["cn_category"].isin(["Deletion", "Loss"])
        & (cnv_data_filtered["total_copy_number"] < 2.0)
    ]

    gain_ids = set(gain_group["model_id"])
    loss_ids = set(loss_group["model_id"])
    patients_all = list(gain_ids | loss_ids)

    for patient in patients_all:
        bnd_file = None
        files_res = [
            file
            for file in os.listdir(bnd_dir_res)
            if file.endswith(".bnd") and patient in file
        ]
        if files_res:
            bnd_file = os.path.join(bnd_dir_res, files_res[0])
        else:
            files_sens = [
                file
                for file in os.listdir(bnd_dir_sens)
                if file.endswith(".bnd") and patient in file
            ]
            if files_sens:
                bnd_file = os.path.join(bnd_dir_sens, files_sens[0])
        if bnd_file is None:
            print(f"No .bnd file found for patient: {patient}")
            continue

        genes = cnv_data_filtered[cnv_data_filtered["model_id"] == patient][
            "gene_symbol"
        ]
        gene_list = genes.tolist()

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0].replace(
            f"_{drug_interest}", ""
        )
        with open(bnd_file, "r") as file:
            content = file.read()

        modified_any = False
        for gene in gene_list:
            print(f"🔍 Processing patient {patient}, gene: {gene}")
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            if patient in gain_ids and patient in loss_ids:
                print(
                    f"Patient {patient} is in both gain and loss groups. Please review."
                )
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in gain_ids:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in loss_ids:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                print(f"Patient {patient} not found in gain or loss CNV groups.")
                continue

            gene_match = pattern.search(content)
            if gene_match:
                print(f"{gene} node found. Replacing...")
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    modified_any = True
            else:
                print(f"No {gene} node found in file for patient {patient_id}")

        if modified_any:
            print(f"{patient_id}: CNV — nodes modified")
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            print(f"{patient_id}: CNV — no nodes modified")
