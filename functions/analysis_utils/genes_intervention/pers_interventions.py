import os
import re


def tailor_bnd_genes_intervention(
    gene_intervs, patients_group, models_folder, drug_interest
):
    """
    Modify the .bnd files for each patient in the patients_group to set the specified gene
    intervention. The gene node is set to logic = 0, rate_up = 0, and rate_down = 1.
    Args:
        gene_intervs (list): list of genes, The genes to be modified in the .bnd files.
        patients_group (list): List of patient identifiers to process.
        models_folder (str): Path to the folder containing the .bnd files.
        drug_interest (str): The drug of interest, used to identify the patient files.
    """
    if isinstance(gene_intervs, str):
        gene_intervs = [gene_intervs]

    modified_files = []
    for patient in patients_group:
        bnd_files = [
            file
            for file in os.listdir(models_folder)
            if file.endswith(".bnd") and patient in file
        ]
        if not bnd_files:
            print(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(models_folder, bnd_files[0])

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0].replace(
            f"_{drug_interest}", ""
        )

        with open(bnd_file, "r") as file:
            content = file.read()

        for gene_interv in gene_intervs:
            print(f"🔍 Processing patient {patient}, gene: {gene_interv}")
            pattern = re.compile(
                rf"Node\s+{re.escape(gene_interv)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            new_gene_block = f"""Node {gene_interv} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""

            gene_match = pattern.search(content)
            if gene_match:
                print(f"{gene_interv} node found. Replacing...")
                content, n_subs = re.subn(pattern, new_gene_block, content)
                if n_subs > 0:
                    with open(bnd_file, "w") as file:
                        file.write(content)
                    print(f"{patient_id}: CNV — nodes modified")
                    modified_files.append(bnd_file)

                else:
                    print(f"{patient_id}: CNV — no substitution made")
            else:
                print(f"No {gene_interv} node found in file for patient {patient_id}")

        with open(bnd_file, "w") as file:
            file.write(content)
        modified_files.append(bnd_file)

    return modified_files
