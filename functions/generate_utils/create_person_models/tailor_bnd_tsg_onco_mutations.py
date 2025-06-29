import pandas as pd
import os
import re


def tailor_bnd_mutat_validation(
    lof_mutations_tsg_filtered,
    gof_mutations_filtered,
    bnd_dir,
    tissue,
):
    patients_onco_mutations = list(gof_mutations_filtered["Sample_ID"])
    patients_tsg_mutations = list(lof_mutations_tsg_filtered["Sample_ID"])

    combined_df = pd.concat([gof_mutations_filtered, lof_mutations_tsg_filtered])
    combined_df = combined_df.set_index("Sample_ID")

    patients_all = list(
        set(patients_onco_mutations + patients_tsg_mutations)
    )  # unique patients

    for patient in patients_all:
        bnd_file = None
        files = [
            file
            for file in os.listdir(bnd_dir)
            if file.endswith(".bnd") and patient in file
        ]
        if files:
            bnd_file = os.path.join(bnd_dir, files[0])
        else:
            print(f"No .bnd file found for patient: {patient}")
            continue
        genes = combined_df.loc[patient, "gene"]

        if isinstance(genes, pd.Series):
            gene_list = genes.tolist()
        else:
            gene_list = [genes]
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
                re.DOTALL
            )
            if patient in patients_onco_mutations and patient in patients_tsg_mutations:
                print(
                    f"Patient {patient} is in both oncogene and TSG mutation groups. Please review."
                )
                # Prefer oncogene block if both
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in patients_onco_mutations:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in patients_tsg_mutations:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                print(f"Patient {patient} not found in oncogene or TSG mutation lists.")
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
            print(f"{patient_id}: mutations — nodes modified")
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            print(f"{patient_id}: mutations — no nodes modified")


def tailor_bnd_tsg_onco_mut(
    mutations_data_filtered_combined, bnd_dir_res, bnd_dir_sens, drug_interest
):
    patients_onco_mutations = []
    patients_tsg_mutations = []
    for patient, value in mutations_data_filtered_combined.iterrows():
        if value["OncogeneHighImpact"]:
            patients_onco_mutations.append(patient)
        if value["TumorSuppressorHighImpact"]:
            patients_tsg_mutations.append(patient)
    patients_all = list(
        set(patients_onco_mutations + patients_tsg_mutations)
    )  # unique patients
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
        if patient not in mutations_data_filtered_combined.index:
            print(f"Patient {patient} not found in mutation data.")
            continue
        genes = mutations_data_filtered_combined.loc[patient, "HugoSymbol"]
        if isinstance(genes, pd.Series):
            gene_list = genes.tolist()
        else:
            gene_list = [genes]
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
                re.DOTALL
            )
            if patient in patients_onco_mutations and patient in patients_tsg_mutations:
                print(
                    f"Patient {patient} is in both oncogene and TSG mutation groups. Please review."
                )
                # Prefer oncogene block if both
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in patients_onco_mutations:
                new_gene_block = f"""Node {gene} {{
        logic = 1;
        rate_up = 1;
        rate_down = 0;
    }}"""
            elif patient in patients_tsg_mutations:
                new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""
            else:
                print(f"Patient {patient} not found in oncogene or TSG mutation lists.")
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
            print(f"{patient_id}: mutations — nodes modified")
            with open(bnd_file, "w") as file:
                file.write(content)
        else:
            print(f"{patient_id}: mutations — no nodes modified")
