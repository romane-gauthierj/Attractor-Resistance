import os
import re
# import sys
# from pathlib import Path

# Add the parent directory of 'identification_patients' to sys.path
# parent_path = Path.cwd().parent  # Go up one level
# sys.path.append(str(parent_path))

# # Now you can import your function
# from identification_patients.get_patients_sens_res import get_patients


def tailor_bnd_genes_intervention(list_genes_interv, res_patients, bnd_dir, tissue):
    for patient in res_patients:
        bnd_file = [
            file
            for file in os.listdir(bnd_dir)
            if file.endswith(".bnd") and patient in file
        ]
        if not bnd_file:
            print(f"No .bnd file found for patient: {patient}")
            continue

        bnd_file = os.path.join(bnd_dir, bnd_file[0])

        patient_id = os.path.splitext(os.path.basename(bnd_file))[0].replace(
            f"_{tissue}", ""
        )
        with open(bnd_file, "r") as file:
            content = file.read()

        modified_any = False
        for gene in list_genes_interv:
            print(f"🔍 Processing patient {patient}, gene: {gene}")
            pattern = re.compile(
                rf"Node\s+{re.escape(gene)}\s*\{{[^{{}}]*?(?:\{{[^{{}}]*?\}}[^{{}}]*?)*\}}",
                re.DOTALL,
            )

            new_gene_block = f"""Node {gene} {{
        logic = 0;
        rate_up = 0;
        rate_down = 1;
    }}"""

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


# tissue_remove = 'Haematopoietic and Lymphoid'
# drug_interest = 'Refametinib'
# tissue_interest = 'PAN_CANCER'

# drug_data = pd.read_csv('../data/drug_sensitivity.csv')
# annotations_models = pd.read_csv('../data/model_list_20250407.csv')
# bnd_dir = '../models/personalized_boolean_Refametinib_PAN_CANCER_intervention/resistant_patient/personalized_boolean_modified/models_gene_expression'

# top_resistant_ids, top_sensitive_ids, drug_tissue_data= get_patients(20, drug_data, annotations_models, drug_interest, tissue_interest = None, tissue_remove = tissue_remove)


# list_genes_interv = ["FOXA1", "BCL2", "BIRC5"]

# tailor_bnd_genes_intervention(list_genes_interv, top_resistant_ids, bnd_dir, tissue_interest)
