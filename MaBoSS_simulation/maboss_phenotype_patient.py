import maboss
import pandas as pd
import os


# loop over max time to identify the one for which the algo converges
def compute_phenotype_table(
    folder_save_results,
    folder_models,
    patient_id,
    input_nodes,
    phenotypes_interest,
    tissue,
):
    model_pers_bnd = f"{folder_models}/{patient_id}_{tissue}.bnd"
    model_pers_cfg = f"{folder_models}/{patient_id}_{tissue}.cfg"
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
    results.to_csv(f"{folder_save_results}/_{patient_id}.csv")
    return results


def compute_phenotype_mean_group_validation(
    stages_groups,
    folder_groups_means,
):
    # compute mean

    for group in stages_groups:
        folder_path = f"{folder_groups_means}/{group}"
        os.makedirs(folder_path, exist_ok=True)
        files = [
            f
            for f in os.listdir(folder_path)
            if f.startswith("_TCGA") and f.endswith(".csv")
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
            print(mean_df)

            mean_df.to_csv(
                f"{folder_groups_means}/{group}/{group}_mean_phenotype_values.csv"
            )

        else:
            print(f"No CSV files found in {folder_path}")
    return mean_df


# def combine_groups_values(folder_to_group, base_path):
#     # Folder and grouping setup
#     all_folders = list(folder_to_group.keys())
#     groups = list(set(folder_to_group.values()))

#     # Storage
#     data_combined = {}

#     # Load data
#     for folder in all_folders:
#         group_label = folder_to_group[folder]
#         folder_path = os.path.join(base_path, folder)

#         for file in os.listdir(folder_path):
#             if file.endswith(".csv"):
#                 df = pd.read_csv(os.path.join(folder_path, file), index_col=0)
#                 for input_name in df.index:
#                     for phenotype in df.columns:
#                         key = (input_name, phenotype)
#                         if key not in data_combined:
#                             data_combined[key] = {g: [] for g in groups}
#                         value = df.at[input_name, phenotype]

#                         data_combined[key][group_label].append(float(value))

#     return data_combined


# def compute_phenotype_mean(group, folder_groups_means, results_mean_folder):
#     # compute mean

#     folder_path = f"{folder_groups_means}"
#     os.makedirs(folder_path, exist_ok=True)
#     files = [f for f in os.listdir(folder_path) if f.endswith(".csv")]

#     dfs = []

#     for file in files:
#         file_path = os.path.join(folder_path, file)
#         df = pd.read_csv(file_path, index_col=0)  # Assuming first column is input IDs
#         dfs.append(df)

#     if dfs:
#         # Now concatenate all dataframes along a new axis to create a 3D structure (like a panel)
#         # Then compute mean across that new axis (axis=0 means across files)
#         mean_df = pd.concat(dfs, axis=0).groupby(level=0).mean()
#         mean_row = mean_df.mean(axis=0)
#         mean_row.name = "Overall_Mean"
#         mean_df = pd.concat([mean_df, mean_row.to_frame().T])
#         print(mean_df)

#         mean_df.to_csv(f"{results_mean_folder}/{group}_mean_phenotype_values.csv")

#     else:
#         print(f"No CSV files found in {folder_path}")
#     return mean_df
