import pandas as pd
import os
import shutil
import logging
import ast

from dotenv import load_dotenv

from functions.generate_models import generate_models_re, pre_process_re

from functions.analysis import downstream_analysis

from functions.validation_utils.validation_Breast import (
    correlate_boolean_predictions_with_gene_signatures
)
from functions.analysis_utils.stats.stats_proba import compute_stats_test_after_ko



logger = logging.getLogger(__name__)

# Load environment variables from .dotenv file
load_dotenv('config.env')




def load_config():

    mutations_data = pd.read_csv(os.getenv('mutations_data_path', 'data/cellmodel_data/mutations_all_20250318.csv'))
    annotations_models = pd.read_csv(os.getenv('annotations_models_path', 'data/model_list_20250407.csv'))
    
    drug_data = pd.read_csv(os.getenv('drug_data_path', 'data/drug_sensitivity.csv'))

    montagud_original_path = os.getenv('montagud_original_path', 'data/montagud_models/Montagud_inter_nodes_data.csv')
    montagud_original_data_df = pd.read_csv(montagud_original_path, header=1)
    montagud_original_data_df = montagud_original_data_df.loc[:, ['Target node', 'Interaction type', 'Source']]




    nodes_montagud_synonyms = pd.read_csv(os.getenv('nodes_montagud_synonyms_path', 'data/montagud_models/nodes_processed.csv'))
    rna_seq_data = pd.read_csv(os.getenv('rna_seq_data_path', 'data/cellmodel_data/rnaseq_merged_20250117/rnaseq_merged_20250117.csv'))
    cnv_data = pd.read_csv(os.getenv('cnv_data_path', 'data/cellmodel_data/cnv_summary_20250207.csv'))
    proteins_data = pd.read_csv(os.getenv('proteins_data_path', 'data/cellmodel_data/proteomics_all_20250211.csv'))
    onco_tsg_data = pd.read_csv(os.getenv('onco_tsg_data_path', 'data/oncogenes_tsg.tsv'), sep='\t')
    
    generic_model_path = os.getenv('generic_model_path', 'data/montagud_models/Montagud2022_Prostate_Cancer')

    number_patients = int(os.getenv('number_patients', 35))
    genetic_intervention = os.getenv('genetic_intervention', None)
    drug_name = os.getenv('drug_name', 'Pictilisib')


    continuous_variable = os.getenv("continuous_variable", 'genes')
    discrete_variable = os.getenv("discrete_variable", 'mutations')
    normalization_techniques = os.getenv("normalization_techniques", 'sigmoid')
    intervention_gene = os.getenv("intervention_gene", None)

    return (mutations_data, annotations_models, drug_data, montagud_original_data_df, nodes_montagud_synonyms, rna_seq_data, cnv_data, proteins_data, onco_tsg_data, generic_model_path, number_patients, genetic_intervention, drug_name,
            continuous_variable, discrete_variable, normalization_techniques, intervention_gene)



def main():
    mutations_data, annotations_models, drug_data, montagud_original_data_df, nodes_montagud_synonyms, rna_seq_data, cnv_data, proteins_data, onco_tsg_data, generic_model_path, number_patients, genetic_intervention, drug_name, continuous_variable, discrete_variable, normalization_technique, interventions = load_config()
    inputs_list = ['EGF', 'FGF', 'TGFb', 'Nutrients', 'Hypoxia', 'Acidosis', 'Androgen', 'TNFalpha', 'SPOP', 'Carcinogen']
    phenotype_interest = ['Invasion', 'Apoptosis']
    folder_name = f'{discrete_variable}_{continuous_variable}'

    type_model=f'{continuous_variable}_models'

    drugs_dict = {'Pictilisib': ['PI3K'], 'AZD7762': ['CHK1_2']}



    if interventions is not None:
        intervention_list = [None] + [[i.strip()] for i in interventions.split(',')]
    else:
        intervention_list = [None]



    drug_targets = drugs_dict[drug_name]


    for intervention_gene in intervention_list:
        print(f"Processing intervention_gene: {intervention_gene}")
    

        # for drug in drugs_dict:
        #     drug_targets = drugs_dict[drug]

        # for norm_technique in normalization_techniques:
        if intervention_gene is not None:
            subdir = f"{'_'.join(drug_targets)}_target_{normalization_technique}/intervention_{'_'.join(intervention_gene)}"
        else: 
            subdir = f"{'_'.join(drug_targets)}_target_{normalization_technique}"


        print(f"Using subdir: {subdir}")


        folder_generic_models = f"analysis/{drug_name}/{folder_name}/{subdir}/models/generic/"
        folder_models = f"analysis/{drug_name}/{folder_name}/{subdir}/models"
        folder_results = f"analysis/{drug_name}/{folder_name}/{subdir}"


        patients_categ = ['resistant', 'sensitive', 'healthy']
        

        for patient_categ in patients_categ:
            os.makedirs(f"analysis/{drug_name}/{folder_name}/{subdir}/results/{patient_categ}", exist_ok=True)
            os.makedirs(f"analysis/{drug_name}/{folder_name}/{subdir}/models/{patient_categ}", exist_ok=True)

        dest_dir = f"analysis/{drug_name}/{folder_name}/{subdir}/models/generic"

        os.makedirs(dest_dir, exist_ok=True)

        # Copy the files
        shutil.copy(f'{generic_model_path}.bnd', dest_dir)
        shutil.copy(f'{generic_model_path}.cfg', dest_dir)


        top_resistant_ids, top_sensitive_ids, top_healthy_ids, montagud_node_model, all_montagud_nodes, rna_seq_data_models_filtered, cnv_data_filtered, df_melted_protein, mutations_data_filtered = pre_process_re(
        montagud_original_data_df.copy(),
        nodes_montagud_synonyms,
        rna_seq_data.copy(),
        cnv_data.copy(),
        mutations_data.copy(),
        number_patients,
        drug_data.copy(),
        annotations_models.copy(),
        drug_name,
        proteins_data.copy(),
        type_model,
        onco_tsg_data,
        tissue_interest=None,
        tissue_remove=None,
    )
        
        patients_ids = top_resistant_ids + top_sensitive_ids + top_healthy_ids


        os.makedirs(f"analysis/{drug_name}/{folder_name}/{subdir}/data_filtered", exist_ok=True)
        rna_seq_data_models_filtered.to_csv(f"analysis/{drug_name}/{folder_name}/{subdir}/data_filtered/rna_seq_data_filtered.csv")
        cnv_data_filtered.to_csv(f"analysis/{drug_name}/{folder_name}/{subdir}/data_filtered/cnv_data_filtered.csv")

        # Convert list to DataFrame then save
        all_montagud_nodes_df = pd.DataFrame(all_montagud_nodes, columns=['gene_symbol'])
        all_montagud_nodes_df.to_csv(f"analysis/{drug_name}/{folder_name}/{subdir}/data_filtered/all_montagud_nodes.csv", index=False)

        generate_models_re(
        normalization_technique,
        discrete_variable,
        continuous_variable,
        folder_generic_models,
        folder_models,
        top_resistant_ids,
        top_sensitive_ids,
        top_healthy_ids,
        drug_name,
        drug_targets,
        phenotype_interest,
        rna_seq_data_models_filtered,
        montagud_node_model,
        cnv_data_filtered,
        mutations_data_filtered,
        type_model,
        df_melted_protein,
        amplif_factor = 100, 
        intervention_gene = intervention_gene,  
        genetic_intervention = genetic_intervention,
    )

        downstream_analysis(
        folder_name, 
        subdir,
        folder_results,
        folder_models,
        drug_name,
        top_resistant_ids,
        top_sensitive_ids,
        top_healthy_ids,
        patients_categ,
        inputs_list,
        phenotype_interest,
        annotations_models,
        intervention_gene=intervention_gene, # instead of None
        list_active_inputs=None,
    )

        proba_phenotype = pd.read_csv(f'analysis/{drug_name}/{folder_name}/{subdir}/results/sensitive_resistant_results/patients_phenot_table.csv', index_col=0)

        results_corr_invasion_df, correlation_invasion_data = correlate_boolean_predictions_with_gene_signatures(False, proba_phenotype, 'Epithelial Mesenchymal Transition', 'Invasion', rna_seq_data, patients_ids)
        results_corr_prolif_df, correlation_prolif_data = correlate_boolean_predictions_with_gene_signatures(False, proba_phenotype, 'G2-M Checkpoint', 'Proliferation', rna_seq_data, patients_ids)
        results_corr_apoptosis_df, correlation_prolif_data = correlate_boolean_predictions_with_gene_signatures(False, proba_phenotype, 'Apoptosis', 'Apoptosis', rna_seq_data, patients_ids)


        results_corr_invasion_df.to_csv(f'analysis/{drug_name}/{folder_name}/{subdir}/results/output/results_corr_invasion_df.csv')
        results_corr_prolif_df.to_csv(f'analysis/{drug_name}/{folder_name}/{subdir}/results/output/results_corr_prolif_df.csv')
        results_corr_apoptosis_df.to_csv(f'analysis/{drug_name}/{folder_name}/{subdir}/results/output/results_corr_apoptosis_df.csv')

        
        if intervention_gene is not None:
            res_values = pd.read_csv(
            f'analysis/{drug_name}/{discrete_variable}_{continuous_variable}/'
            f"{'_'.join(drug_targets)}_target_{normalization_technique}/results/resistant/combined_results.csv",
            index_col=0
        )   
            sens_values = pd.read_csv(
            f'analysis/{drug_name}/{discrete_variable}_{continuous_variable}/'
            f"{'_'.join(drug_targets)}_target_{normalization_technique}/results/sensitive/combined_results.csv",
            index_col=0
        )
            


            # after KO
            res_values_knockout_ctnnb1 = pd.read_csv(
            f'analysis/{drug_name}/{discrete_variable}_{continuous_variable}/{subdir}/results/resistant/combined_results.csv',
            index_col=0
        )
            sens_values_knockout_ctnnb1 = pd.read_csv(
            f'analysis/{drug_name}/{discrete_variable}_{continuous_variable}/{subdir}/results/sensitive/combined_results.csv',
            index_col=0
            )

            output_folder = f'analysis/{drug_name}/{discrete_variable}_{continuous_variable}/{subdir}/results/output'



            compute_stats_test_after_ko(res_values,  sens_values, res_values_knockout_ctnnb1, sens_values_knockout_ctnnb1 , output_folder, 'EGF','Invasion')
            
            

            logger.debug(f"Completed {drug_name} with {normalization_technique} normalization successfully!")


if __name__ == "__main__":
    main()
