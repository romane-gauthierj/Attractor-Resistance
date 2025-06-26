import os
import re

def tailor_bnd_genes_intervention(list_genes_interv, res_patients, bnd_dir, tissue):
    if isinstance(list_genes_interv, str):
        list_genes_interv = [list_genes_interv]
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