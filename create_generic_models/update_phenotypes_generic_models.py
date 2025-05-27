import pandas as pd
import numpy as np
import os

import re


def generic_models_update_phenotypes(
    phenotype_interest, original_data_dir, results_dir
):
    for filename in os.listdir(original_data_dir):
        file_path = os.path.join(original_data_dir, filename)
        # Only process files ending with .cfg
        if os.path.isfile(file_path) and filename.endswith(".cfg"):
            with open(file_path, "r") as file:
                content = file.read()
            # Find all nodes that have an is_internal assignment (0 or 1)
            node_names = set(re.findall(r"(\w+)\.is_internal\s*=\s*[01]", content))
            for node in node_names:
                # Set new value to 0 if node is in phenotype_interest, else 1
                new_internal_val = "0" if node in phenotype_interest else "1"
                # Create regex pattern to find this node's is_internal assignment
                pattern = re.compile(rf"({re.escape(node)}\.is_internal\s*=\s*)[01]")
                # Replacement string using group 1 and the new value
                replacement = r"\g<1>" + new_internal_val
                # Substitute all occurrences for this node
                content, n_subs = pattern.subn(replacement, content)
                if n_subs > 0:
                    print(
                        f"Updated {node}.is_internal={new_internal_val} in {filename}"
                    )
            # Save the modified content back to the results directory
            modified_file_path = os.path.join(results_dir, filename)
            with open(modified_file_path, "w") as file:
                file.write(content)
            print(f"Modified and saved: {modified_file_path}")


# node = "Proliferation"
# content = "Proliferation.is_internal=1;\nSomeOtherNode.is_internal=1"
# new_internal_val = "0"

# pattern = re.compile(rf"({re.escape(node)}\.is_internal\s*=\s*)[01]")
# replacement = r"\g<1>" + new_internal_val

# content, n_subs = pattern.subn(replacement, content)
# print(content)
# print("Subs made:", n_subs)


# # Modify phenotypes in the generic cfg files (then reuse them to create the personalized boolean networks)

# phenotype_interest = [
#     "Proliferation",
#     "Invasion",
#     "DNA_Repair",
#     "Migration",
#     "Apoptosis",
# ]

# original_data_dir = "../models/personalized_boolean_Refametinib_PAN_CANCER/resistant_patient/generic_models"
# results_dir = "../models/personalized_boolean_Refametinib_PAN_CANCER/resistant_patient/generic_models"

# for filename in os.listdir(original_data_dir):
#     file_path = os.path.join(original_data_dir, filename)

#     if os.path.isfile(file_path) and filename.endswith(".cfg"):
#         with open(file_path, "r") as file:
#             content = file.read()

#         # Find all nodes with is_internal statements
#         node_names = re.findall(r"(\w+)\.is_internal\s*=\s*[01]", content)

#         for node in node_names:
#             # Determine if node is a phenotype or not
#             new_internal_val = "1" if node in phenotype_interest else "0"

#             # Compile pattern to match the node's is_internal line
#             pattern = re.compile(
#                 rf"({re.escape(node)}\.is_internal\s*=\s*)[01]", re.DOTALL
#             )

#             # Replacement using backreference \1 and new value
#             replacement = r"\1" + new_internal_val

#             content, n_subs = pattern.subn(replacement, content)

#             if n_subs > 0:
#                 print(f"Updated {node}.is_internal={new_internal_val} in {filename}")

#         # Write the modified content back to file
#         modified_file_path = os.path.join(results_dir, filename)
#         with open(modified_file_path, "w") as file:
#             file.write(content)

#         print(f"Modified and saved: {modified_file_path}")
