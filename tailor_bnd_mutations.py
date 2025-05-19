# Pipeline to tailor bnd files according to the KRAS and EGFR mutations
# ???? DID NOT WORK ???? !!!! cannot replace logical in bnd file for the KRAS and EGFR mutations


import os
import re
from identify_mutations_patients import identif_mutations_kras_egfr


def personalized_patients_mutations_bnds(mutations_data, patients_ids, bnd_dir, drug_interest):
    # New KRAS and EGFR blocks
    new_kras_block = """Node KRAS {
    logic = 1;
    rate_up = $u_KRAS;
    rate_down = 0;
    }"""

    new_egfr_block = """Node EGFR {
    logic = 1;
    rate_up = $u_EGFR;
    rate_down = 0;
    }"""

    # Improved regex (handles more complex logic and nested braces)
    kras_pattern = re.compile(r'Node\s+KRAS\s*\{[^{}]*?(?:\{[^{}]*?\}[^{}]*?)*\}', re.DOTALL)
    egfr_pattern = re.compile(r'Node\s+EGFR\s*\{[^{}]*?(?:\{[^{}]*?\}[^{}]*?)*\}', re.DOTALL)

    # Loop through files
    for filename in os.listdir(bnd_dir):
        if filename.endswith(f'{drug_interest}.bnd'):
            file_path = os.path.join(bnd_dir, filename)
            #patient_id = os.path.splitext(filename)[0].replace('_AZD8931', '')
            patient_id = os.path.splitext(filename)[0].replace(f'_{drug_interest}', '')
        

            with open(file_path, 'r') as file:
                content = file.read()
            
            original_content = content
            modified = False

            # Check if the file has a KRAS mutation
            mutations_kras_ids, mutations_egfr_ids = identif_mutations_kras_egfr(mutations_data, patients_ids)
        
           # why this is not working ??

            if (patient_id in mutations_kras_ids):
                print(patient_id)
                print('kras')
                
                kras_match = kras_pattern.search(content)
                print(kras_match)
                print()
                if kras_match:
                    print(f"KRAS node before replacement: {kras_match.group(0)}")
                    # Replace KRAS block if found
                    content, n_subs = re.subn(kras_pattern, new_kras_block, content)
                else:
                    print(f"No KRAS node found for patient {patient_id}")
                
                if n_subs > 0:
                    print(f'{patient_id}: KRAS mutation — node modified')
                    modified = True
                else:
                    print(f'{patient_id}: KRAS mutation — ⚠️ KRAS node not found!')

            # Check if the file has an EGFR mutation
            if patient_id in mutations_egfr_ids:
                egfr_match = egfr_pattern.search(content)
                if egfr_match:
                    print(f"EGFR node before replacement: {egfr_match.group(0)}")
                    # Replace EGFR block if found
                    content, n_subs = re.subn(egfr_pattern, new_egfr_block, content)
                else:
                    print(f"No EGFR node found for patient {patient_id}")

                if n_subs > 0:
                    print(f'{patient_id}: EGFR mutation — node modified')
                    modified = True
                else:
                    print(f'{patient_id}: EGFR mutation — ⚠️ EGFR node not found!')
                    
            # If modified, save the new content to the file
            if modified and content != original_content:
                print(f"Modified content for {filename}:")
                print(content)  # Debugging - print modified content of the file
                
                with open(file_path, 'w') as file:
                    file.write(content)