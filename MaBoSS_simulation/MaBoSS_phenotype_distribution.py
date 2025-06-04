import matplotlib.pyplot as plt
import maboss
import os
import pandas as pd
import numpy as np


# change with the phenotypes columns of interest


# MODEL_BND = "../models/personalized_boolean_Refametinib_PAN_CANCER/resistant_patient/personalized_boolean_modified/models_gene_expression/SIDM01120_Refametinib.bnd"
# MODEL_CFG = "../models/personalized_boolean_Refametinib_PAN_CANCER/resistant_patient/personalized_boolean_modified/models_gene_expression/SIDM01120_Refametinib.cfg"

# inputs_list = ["EGF", "FGF", "TNFa"]


# model = maboss.load(MODEL_BND, MODEL_CFG)
# last_states = pd.DataFrame()
# name = "SIDM01120"
# for active_input in inputs_list:
#     model.network.set_istate(active_input, [0, 1])
#     for inactive_input in inputs_list:
#         if inactive_input != active_input:
#             model.network.set_istate(inactive_input, [1, 0])
#             res_model = model.run()
#             last_state = res_model.get_last_states_probtraj()
#             last_state.index = [name]
#             last_states = pd.concat([last_states, last_state])
# print(last_states)


# def compute_phenotypes_distribution(folder,dir_files, inputs_list, group_label, drug_interest):
#     output_dir = f'{folder}/{group_label}_results/only_gene_expression/single_input_on/phenotype_distribution_patients'
#     cfg_files = [f for f in os.listdir(dir_files) if f.endswith(f"{drug_interest}.cfg")]
#     patients_data_dict = {}
#     for cfg_file in cfg_files:
#         cfg_path = os.path.join(dir_files, cfg_file)
#         base_name = os.path.splitext(cfg_file)[0]
#         bnd_path = os.path.join(dir_files, base_name + ".bnd")
#         #model = maboss.load(bnd_path, cfg_path)
#         patients_results = pd.DataFrame()
#         print(f"\n--- Results for patient: {base_name} ---")

#         last_state = result.get_last_states_probtraj()
#         for active_input in inputs_list:
#             last_state.index = [f"{active_input}_ON"]
#         for inactive_input in inputs_list:
#             if inactive_input != active_input:
#                 model.network.set_istate(inactive_input, [1, 0])

#         patients_results = pd.concat([patients_results, last_state])
#         expected_cols = ['<nil>', 'Apoptosis', 'Proliferation', 'Metastasis']
#         available_cols = [col for col in expected_cols if col in patients_results.columns]

#         patients_results = patients_results[available_cols]
#         patients_results = patients_results.fillna(0)

#         output_path = os.path.join(output_dir, f"{base_name}_phenotypes.csv")
#         patients_results.to_csv(output_path)

#         patients_data_dict[base_name] = patients_results
#         print(patients_results)
#     return patients_data_dict


# why this function is not working?


# def compute_phenotypes_distribution(
#     phenotypes_list, folder, dir_files, inputs_list, drug_interest, group_label=None
# ):
#     if group_label is not None:
#         output_dir = f"{folder}/{group_label}_results/only_gene_expression/single_input_on/phenotype_distribution_patients"
#     else:
#         output_dir = f"{folder}/only_gene_expression/single_input_on/phenotype_distribution_patients"

#     cfg_files = [f for f in os.listdir(dir_files) if f.endswith(f"{drug_interest}.cfg")]
#     patients_data_dict = {}
#     for cfg_file in cfg_files:
#         cfg_path = os.path.join(dir_files, cfg_file)
#         base_name = os.path.splitext(cfg_file)[0]
#         bnd_path = os.path.join(dir_files, base_name + ".bnd")

#         if not os.path.exists(bnd_path):
#             print(f"Missing .bnd file: {bnd_path}")
#             continue
#         if not os.path.exists(cfg_path):
#             print(f"Missing .cfg file: {cfg_path}")
#             continue

#         patients_results = pd.DataFrame()

#         print(f"\n--- Results for patient: {base_name} ---")

#         for active_input in inputs_list:
#             try:
#                 model = maboss.load(bnd_path, cfg_path)  # Reload fresh
#                 model.network.set_istate(active_input, [0, 1])
#                 for inactive_input in inputs_list:
#                     if inactive_input != active_input:
#                         model.network.set_istate(inactive_input, [1, 0])
#                 result = model.run()

#                 last_state = result.get_last_states_probtraj()

#             except Exception as e:
#                 print(
#                     f"Error during simulation for {base_name} ({active_input} ON): {e}"
#                 )
#                 continue
#             last_state.index = [f"{active_input}_ON"]
#             patients_results = pd.concat([patients_results, last_state])
#             # available_cols = [
#             #     col for col in expected_cols if col in patients_results.columns
#             # ]
#             # patients_results = patients_results[available_cols]
#             # patients_results = patients_results.fillna(0)
#         print("patients columns are:", patients_results.columns)

#         for phenotype in phenotypes_list:
#             if phenotype not in patients_results.columns:
#                 patients_results[phenotype] = 0.0

#         patients_results = patients_results[phenotypes_list]

#         # Save to CSV
#         os.makedirs(output_dir, exist_ok=True)
#         output_path = os.path.join(output_dir, f"{base_name}_phenotypes.csv")
#         patients_results.to_csv(output_path)

#         # Add to dictionary
#         patients_data_dict[base_name] = patients_results
#         # print(patients_results)

#         print(patients_results)

#     return patients_data_dict



# def compute_mean_patients(phenotype_list, dic_patients_data):
#     stats_results_data = {}
#     results = {}
#     flags_end = [False, False, False]
#     df_results_mean = pd.DataFrame(
#         index=[
#             "EGF_ON",
#             "FGF_ON",
#             "TGFb_ON",
#             "Nutrients_ON",
#             "Hypoxia_ON",
#             "Acidosis_ON",
#             "Androgen_ON",
#             "TNFalpha_ON",
#             "Carcinogen_ON",
#         ],
#         columns=phenotype_list,
#     )
#     df_results_std = pd.DataFrame(
#         index=[
#             "EGF_ON",
#             "FGF_ON",
#             "TGFb_ON",
#             "Nutrients_ON",
#             "Hypoxia_ON",
#             "Acidosis_ON",
#             "Androgen_ON",
#             "TNFalpha_ON",
#             "Carcinogen_ON",
#         ],
#         columns=phenotype_list,
#     )

#     for idx_patient, (patient_id, patient_data) in enumerate(
#         dic_patients_data.items(), 1
#     ):
#         if len(dic_patients_data) == idx_patient:
#             flags_end[0] = True
#         for idx_condition, (condition_name, phenotypes) in enumerate(
#             patient_data.iterrows(), 1
#         ):
#             if patient_data.shape[0] == idx_condition:
#                 flags_end[1] = True
#             for idx_phenotype, (phenotype_name, phenotype_value) in enumerate(
#                 phenotypes.items(), 1
#             ):
#                 if len(phenotypes.values) == idx_phenotype:
#                     flags_end[2] = True
#                 if condition_name not in results.keys():
#                     results[condition_name] = {}
#                 if phenotype_name not in results[condition_name].keys():
#                     results[condition_name][phenotype_name] = []

#                 results[condition_name][phenotype_name].append(phenotype_value)
#                 temp_result = results[condition_name][phenotype_name]

#                 # Stats test table
#                 if condition_name not in stats_results_data.keys():
#                     stats_results_data[condition_name] = {}
#                 if phenotype_name not in stats_results_data[condition_name].keys():
#                     stats_results_data[condition_name][phenotype_name] = []

#                 stats_results_data[condition_name][phenotype_name].append(
#                     phenotype_value
#                 )

#                 if False not in flags_end:
#                     mean = np.mean(temp_result)
#                     std = np.std(temp_result, ddof=1)

#                     # results[condition_name][phenotype_name] = [mean, std]
#                     # Put mean and std in dataframe

#                     if pd.isna(df_results_mean.at[condition_name, phenotype_name]):
#                         df_results_mean.at[condition_name, phenotype_name] = mean
#                     else:
#                         df_results_mean.at[condition_name, phenotype_name] += mean

#                     # Same for std
#                     if pd.isna(df_results_std.at[condition_name, phenotype_name]):
#                         df_results_std.at[condition_name, phenotype_name] = std
#                     else:
#                         df_results_std.at[condition_name, phenotype_name] += std

#     stats_results_data_df = pd.DataFrame(stats_results_data)
#     print("Stats test table values")
#     print(stats_results_data_df)

#     return df_results_mean.round(2), df_results_std.round(2), stats_results_data_df
