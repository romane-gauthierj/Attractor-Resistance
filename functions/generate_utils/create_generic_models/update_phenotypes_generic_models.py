import pandas as pd
import numpy as np
import os

import re


def generic_models_update_phenotypes(phenotype_interest, folder_models):
    for filename in os.listdir(folder_models):
        file_path = os.path.join(folder_models, filename)
        # Only process files ending with .cfg
        if os.path.isfile(file_path) and filename.endswith(".cfg"):
            with open(file_path, "r") as file:
                content = file.read()
            modified = False
            for node in phenotype_interest:
                # Regex pattern for is_internal assignment
                pattern = re.compile(rf"({re.escape(node)}\.is_internal\s*=\s*)[01]")
                if re.search(pattern, content):
                    # Update existing assignment to 0
                    content, n_subs = pattern.subn(rf"\g<1>0", content)
                    if n_subs > 0:
                        # print(f"Updated {node}.is_internal=0 in {filename}")
                        modified = True
                else:
                    # Add new assignment at the end if not present
                    content += f"\n{node}.is_internal=0"
                    # print(f"Added {node}.is_internal=0 to {filename}")
                    modified = True
            # Save only if modified
            if modified:
                modified_file_path = os.path.join(folder_models, filename)
                with open(modified_file_path, "w") as file:
                    file.write(content)
                    print(f"Modified and saved: {modified_file_path}")
