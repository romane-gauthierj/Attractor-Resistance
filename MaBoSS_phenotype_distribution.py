import matplotlib.pyplot as plt
import maboss
import os
import pandas as pd
import numpy as np


def compute_phenotypes_distribution(folder,dir_files, inputs_list, group_label, drug_interest):
    output_dir = f'{folder}/{group_label}_results/only_gene_expression/single_input_on/phenotype_distribution_patients'
    cfg_files = [f for f in os.listdir(dir_files) if f.endswith(f"{drug_interest}.cfg")]
    patients_data_dict = {}
    for cfg_file in cfg_files:
        cfg_path = os.path.join(dir_files, cfg_file)
        base_name = os.path.splitext(cfg_file)[0]
        bnd_path = os.path.join(dir_files, base_name + ".bnd")
        #model = maboss.load(bnd_path, cfg_path)
        patients_results = pd.DataFrame()
        print(f"\n--- Results for patient: {base_name} ---")

        for active_input in inputs_list:
            model = maboss.load(bnd_path, cfg_path)
            model.network.set_istate(active_input, [0, 1])

            for inactive_input in inputs_list:
                if inactive_input != active_input:
                    model.network.set_istate(inactive_input, [1, 0])

            result = model.run()
            last_state = result.get_last_states_probtraj()
            last_state.index = [f"{active_input}_ON"]
            patients_results = pd.concat([patients_results, last_state])
            expected_cols = ['<nil>', 'Apoptosis', 'Proliferation', 'Metastasis']
            available_cols = [col for col in expected_cols if col in patients_results.columns]
            patients_results = patients_results[available_cols]
            patients_results = patients_results.fillna(0)

                # Save to CSV
            output_path = os.path.join(output_dir, f"{base_name}_phenotypes.csv")
            patients_results.to_csv(output_path)

            # Add to dictionary
            patients_data_dict[base_name] = patients_results
            # print(patients_results)

        patients_data_dict[base_name] = patients_results
        print(patients_results)
        

    return patients_data_dict




def compute_mean_patients(dic_patients_data):
    stats_results_data = {}
    results = {}
    flags_end = [False, False, False]
    df_results_mean = pd.DataFrame(index=['EGF_ON', 'FGF_ON', 'TGFb_ON', 'Nutrients_ON', 'Hypoxia_ON', 'Acidosis_ON', 'Androgen_ON', 'TNFalpha_ON', 'Carcinogen_ON'], columns=['<nil>', 'Apoptosis', 'Proliferation', 'Metastasis'])
    df_results_std = pd.DataFrame(index=['EGF_ON', 'FGF_ON', 'TGFb_ON', 'Nutrients_ON', 'Hypoxia_ON', 'Acidosis_ON', 'Androgen_ON', 'TNFalpha_ON', 'Carcinogen_ON'], columns=['<nil>', 'Apoptosis', 'Proliferation', 'Metastasis'])

    for idx_patient, (patient_id, patient_data) in enumerate(dic_patients_data.items(), 1):
        if len(dic_patients_data) == idx_patient:
            flags_end[0] = True
        for idx_condition, (condition_name, phenotypes) in enumerate(patient_data.iterrows(), 1):
            if patient_data.shape[0] == idx_condition:
                flags_end[1] = True
            for idx_phenotype, (phenotype_name, phenotype_value) in enumerate(phenotypes.items(), 1):
                if len(phenotypes.values) == idx_phenotype:
                    flags_end[2] = True
                if condition_name not in results.keys():
                    results[condition_name] = {}
                if phenotype_name not in results[condition_name].keys():
                    results[condition_name][phenotype_name] = []

                results[condition_name][phenotype_name].append(phenotype_value)
                temp_result = results[condition_name][phenotype_name]


                # Stats test table
                if condition_name not in stats_results_data.keys():
                    stats_results_data[condition_name] = {}
                if phenotype_name not in stats_results_data[condition_name].keys():
                    stats_results_data[condition_name][phenotype_name] = []
                
                stats_results_data[condition_name][phenotype_name].append(phenotype_value)

                if False not in flags_end:
                    mean = np.mean(temp_result)
                    std = np.std(temp_result, ddof=1)

                    #results[condition_name][phenotype_name] = [mean, std]
                    # Put mean and std in dataframe

                    if pd.isna(df_results_mean.at[condition_name, phenotype_name]):
                        df_results_mean.at[condition_name, phenotype_name] = mean
                    else:
                        df_results_mean.at[condition_name, phenotype_name] += mean
                

                    # Same for std
                    if pd.isna(df_results_std.at[condition_name, phenotype_name]):
                        df_results_std.at[condition_name, phenotype_name] = std
                    else:
                        df_results_std.at[condition_name, phenotype_name] += std

    stats_results_data_df = pd.DataFrame(stats_results_data)
    print('Stats test table values')
    print(stats_results_data_df)

    return df_results_mean.round(2), df_results_std.round(2), stats_results_data_df




