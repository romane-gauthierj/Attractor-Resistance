import pandas as pd
import os
import re
import logging
logger = logging.getLogger(__name__)


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
            logger.debug(f"No .bnd file found for patient: {patient}")
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
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            if patient in gain_ids and patient in loss_ids:
                logger.warning(
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
                logger.debug(f"Patient {patient} not found in gain or loss CNV groups.")
                continue

            gene_match = pattern.search(content)
            if gene_match:
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    modified_any = True
            else:
                logger.debug(f"No {gene} node found in file for patient {patient_id}")

        if modified_any:
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            logger.debug(f"{patient_id}: CNV — no nodes modified")









def tailor_bnd_cnv_cm(cnv_data_filtered, folder_models):
    gain_group = cnv_data_filtered[
        cnv_data_filtered["effect"]=='amplification'
    ]

    loss_group = cnv_data_filtered[
        cnv_data_filtered["effect"]=='deletion'
    ]

    gain_ids = set(gain_group["model_id"])
    loss_ids = set(loss_group["model_id"])
    patients_all = list(gain_ids | loss_ids)

    for patient in patients_all:
        files_categ = [
            file
            for file in os.listdir(folder_models)
            if file.endswith(".bnd") and patient in file
        ]
        if not files_categ:
            logger.debug(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(folder_models, files_categ[0])

        if bnd_file is None:
            logger.debug(f"No .bnd file found for patient: {patient}")
            continue

        genes = cnv_data_filtered[cnv_data_filtered["model_id"] == patient][
            "gene_symbol"
        ]
        gene_list = genes.tolist()

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0]
        with open(bnd_file, "r") as file:
            content = file.read()

        modified_any = False
        for gene in gene_list:
            
            # Check if this specific patient-gene combination is in both gain and loss
            patient_gene_in_gain = len(gain_group[
                (gain_group["model_id"] == patient) & 
                (gain_group["gene_symbol"] == gene)
            ]) > 0
            
            patient_gene_in_loss = len(loss_group[
                (loss_group["model_id"] == patient) & 
                (loss_group["gene_symbol"] == gene)
            ]) > 0
            
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            if patient_gene_in_gain and patient_gene_in_loss:
                logger.debug(
                    f"Patient {patient} with gene {gene} is in both gain and loss groups. Please review."
                )
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient_gene_in_gain:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient_gene_in_loss:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                logger.debug(f"Patient {patient} with gene {gene} not found in gain or loss CNV groups.")
                continue

            gene_match = pattern.search(content)
            if gene_match:
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    modified_any = True
            else:
                logger.debug(f"No {gene} node found in file for patient {patient_id}")

        if modified_any:
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            logger.debug(f"{patient_id}: CNV — no nodes modified")
