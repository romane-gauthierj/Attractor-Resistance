import maboss
import pandas as pd
import os
import glob
import ast
from collections import defaultdict
import json  # for explicit serialization


def compute_phenotype_table_generic(
    folder_model,
    input_nodes,
    phenotypes_interest,
    results_dir,
):
    model_bnd_list = glob.glob(os.path.join(folder_model, "*.bnd"))
    model_cfg_list = glob.glob(os.path.join(folder_model, "*.cfg"))

    if not model_bnd_list or not model_cfg_list:
        raise FileNotFoundError(
            "BND or CFG file not found in the folder_model directory."
        )

    model_bnd = model_bnd_list[0]
    model_cfg = model_cfg_list[0]
    model = maboss.load(model_bnd, model_cfg)

    results = pd.DataFrame(index=input_nodes, columns=phenotypes_interest)

    threshold_diff = 0.01
    initial_max_time = 5
    max_time_limit = 30

    for active_node in input_nodes:
        model.network.set_istate(active_node, [0, 1])
        for inactive_node in input_nodes:
            if inactive_node != active_node:
                model.network.set_istate(inactive_node, [1, 0])
        converged = False
        current_max_time = initial_max_time
        while not converged and current_max_time < max_time_limit:
            model.update_parameters(
                time_tick=0.1, max_time=current_max_time, sample_count=500
            )
            res_lung_pers = model.run()
            last_state = res_lung_pers.get_nodes_probtraj()
            if len(last_state) < 10:
                print(
                    f"[{active_node}] Warning: Not enough data for convergence check."
                )
                break
            diff = last_state.diff().abs()
            max_change_recent = diff.tail(10).max()
            # print(f"Max change over last 10 time steps for {active_node}:\n{max_change_recent}\n")
            if (max_change_recent > threshold_diff).any():
                current_max_time += 1
            else:
                converged = True
        if not converged:
            print(
                f"[{active_node}] Did not fully converge by max_time = {max_time_limit}."
            )
        final_probs = last_state.iloc[-1]
        # print(f"[{active_node}] Final probabilities:\n{final_probs}\n")
        for phenotype in phenotypes_interest:
            if phenotype in final_probs.index:
                results.loc[active_node, phenotype] = final_probs[phenotype]

        results.loc["Overal Mean"] = results.mean()
    os.makedirs(results_dir, exist_ok=True)
    results.to_csv(f"{results_dir}/phenotype_distribution_generic.csv")
    return results


# loop over max time to identify the one for which the algo converges
def compute_phenotype_table(
    folder_save_results,
    folder_models,
    patient_id,
    input_nodes,
    phenotypes_interest,
    context_label,
):
    # context_label: can be either the tissue name or the drug name

    # model_pers_bnd = f"{folder_models}/pers_models/{patient_id}_{context_label}.bnd"
    # model_pers_cfg = f"{folder_models}/pers_models/{patient_id}_{context_label}.cfg"

    model_pers_bnd = f"{folder_models}/{patient_id}_{context_label}.bnd"
    model_pers_cfg = f"{folder_models}/{patient_id}_{context_label}.cfg"

    if not os.path.exists(model_pers_bnd) or not os.path.exists(model_pers_cfg):
        print(f"Missing model files for {patient_id}")
        return None

    model_pers_lung = maboss.load(model_pers_bnd, model_pers_cfg)

    results = pd.DataFrame(index=input_nodes, columns=phenotypes_interest)

    threshold_diff = 0.01
    initial_max_time = 5
    max_time_limit = 30

    for active_node in input_nodes:
        model_pers_lung.network.set_istate(active_node, [0, 1])
        for inactive_node in input_nodes:
            if inactive_node != active_node:
                model_pers_lung.network.set_istate(inactive_node, [1, 0])

        converged = False
        current_max_time = initial_max_time
        while not converged and current_max_time < max_time_limit:
            model_pers_lung.update_parameters(
                time_tick=0.1, max_time=current_max_time, sample_count=500
            )
            res_lung_pers = model_pers_lung.run()

            last_state = res_lung_pers.get_nodes_probtraj()
            if len(last_state) < 10:
                print(
                    f"[{active_node}] Warning: Not enough data for convergence check."
                )
                break
            diff = last_state.diff().abs()
            max_change_recent = diff.tail(10).max()
            # print(f"Max change over last 10 time steps for {active_node}:\n{max_change_recent}\n")
            if (max_change_recent > threshold_diff).any():
                current_max_time += 1
            else:
                converged = True
        if not converged:
            print(
                f"[{active_node}] Did not fully converge by max_time = {max_time_limit}."
            )
        final_probs = last_state.iloc[-1]
        # print(f"[{active_node}] Final probabilities:\n{final_probs}\n")
        for phenotype in phenotypes_interest:
            if phenotype in final_probs.index:
                results.loc[active_node, phenotype] = final_probs[phenotype]
    results.to_csv(f"{folder_save_results}/{patient_id}.csv")
    return results


def compute_phenotype_table_custom_inputs(
    folder_save_results,
    folder_models,
    patient_id,
    input_nodes,
    phenotypes_interest,
    context_label,
    active_inputs,
):
    # function to run MaBoSS with more than one Input ON
    model_pers_bnd = f"{folder_models}/pers_models/{patient_id}_{context_label}.bnd"
    model_pers_cfg = f"{folder_models}/pers_models/{patient_id}_{context_label}.cfg"

    if not os.path.exists(model_pers_bnd) or not os.path.exists(model_pers_cfg):
        print(f"Missing model files for {patient_id}")
        return None

    model_pers_lung = maboss.load(model_pers_bnd, model_pers_cfg)
    results = pd.DataFrame(columns=phenotypes_interest)

    threshold_diff = 0.01
    initial_max_time = 5
    max_time_limit = 30

    # Ensure active_inputs is a list of lists
    if isinstance(active_inputs[0], str):
        active_inputs = [active_inputs]

    for inputs_on in active_inputs:
        row_key = ",".join(sorted(inputs_on))

        # Set ON/OFF states
        for node in input_nodes:
            if node in inputs_on:
                model_pers_lung.network.set_istate(node, [0, 1])
            else:
                model_pers_lung.network.set_istate(node, [1, 0])

        converged = False
        current_max_time = initial_max_time
        while not converged and current_max_time < max_time_limit:
            model_pers_lung.update_parameters(
                time_tick=0.1, max_time=current_max_time, sample_count=500
            )
            res_lung_pers = model_pers_lung.run()
            last_state = res_lung_pers.get_nodes_probtraj()
            if len(last_state) < 10:
                print(
                    f"[{active_inputs}] Warning: Not enough data for convergence check."
                )
                break
            diff = last_state.diff().abs()
            max_change_recent = diff.tail(10).max()
            # print(f"Max change over last 10 time steps for {active_node}:\n{max_change_recent}\n")
            if (max_change_recent > threshold_diff).any():
                current_max_time += 1
            else:
                converged = True
        if not converged:
            print(
                f"[{active_inputs}] Did not fully converge by max_time = {max_time_limit}."
            )
        final_probs = last_state.iloc[-1]
        results.loc[row_key] = [final_probs.get(ph, None) for ph in phenotypes_interest]

    results.to_csv(f"{folder_save_results}/{patient_id}.csv")
    return results


def compute_phenotype_mean_group_validation(
    stages_groups,
    folder_groups_means,
):
    for group in stages_groups:
        folder_path = f"{folder_groups_means}/{group}"
        os.makedirs(folder_path, exist_ok=True)
        files = [
            f
            for f in os.listdir(folder_path)
            if f.startswith("TCGA") and f.endswith(".csv")
        ]

        dfs = []

        for file in files:
            file_path = os.path.join(folder_path, file)
            df = pd.read_csv(
                file_path, index_col=0
            )  # Assuming first column is input IDs
            dfs.append(df)

        if dfs:
            # Now concatenate all dataframes along a new axis to create a 3D structure (like a panel)
            # Then compute mean across that new axis (axis=0 means across files)
            mean_df = pd.concat(dfs, axis=0).groupby(level=0).mean()
            mean_row = mean_df.mean(axis=0)
            mean_row.name = "Overall_Mean"
            mean_df = pd.concat([mean_df, mean_row.to_frame().T])

            mean_df.to_csv(
                f"{folder_groups_means}/{group}/{group}_mean_phenotype_values.csv"
            )

        else:
            print(f"No CSV files found in {folder_path}")


def compute_mean_phenotype_values(df_patients_combined):
    # compute mean for each patient in the resistant group
    conditions = df_patients_combined.index
    phenotypes = df_patients_combined.columns
    patient_mean = pd.DataFrame(index=conditions, columns=phenotypes)

    for condition, values in df_patients_combined.iterrows():
        for phenotype in df_patients_combined.columns:
            df_patients_values = df_patients_combined.loc[condition][phenotype]
            df_patients_values = ast.literal_eval(df_patients_values)
            mean = sum(df_patients_values) / len(df_patients_values)
            patient_mean.loc[condition][phenotype] = mean
    return patient_mean


def collect_group_data(group_folder_path, patients_id):
    # Aggregates all patient CSVs in a group folder (e.g., resistant or sensitive).
    # For each input/phenotype pair, collects all values across patients into a list.
    # this is run separately for each group (resistant and sensitive)

    combined_data = defaultdict(lambda: defaultdict(list))

    # List all files and print for debugging
    all_files = os.listdir(group_folder_path)

    prefix = patients_id[0][:3]

    # Adjust file filter if needed
    patient_files = [
        f for f in all_files if f.endswith(".csv") and f.startswith(prefix)
    ]

    for file in patient_files:
        file_path = os.path.join(group_folder_path, file)

        df = pd.read_csv(file_path, index_col=0)

        for input_name in df.index:
            for phenotype in df.columns:
                value = df.at[input_name, phenotype]
                try:
                    combined_data[input_name][phenotype].append(float(value))
                except Exception as e:
                    print(f"Error converting value '{value}' in {file_path}: {e}")

    combined_results = pd.DataFrame.from_dict(combined_data, orient="index")
    # Serialize lists as JSON strings for CSV output
    combined_results_serialized = combined_results.apply(
        lambda x: json.dumps(x) if isinstance(x, list) else x
    )
    output_path = os.path.join(group_folder_path, "combined_results.csv")
    combined_results_serialized.to_csv(output_path)
