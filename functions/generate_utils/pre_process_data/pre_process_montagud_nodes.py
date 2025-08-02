import pandas as pd


# def process_montagud_nodes_synonyms(montagud_synonyms_data):

#     montagud_node_synonyms = montagud_synonyms_data.dropna(subset=['Node'])
#      # Create multiple rows at once
#     new_rows = pd.DataFrame({
#         'Node': ['cFLAR', 'eEF2', 'eEF2K', 'Rheb'],  
#         'HGNC symbols': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'], 
#         'unique': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'] 
#     })

#     montagud_node_synonyms = pd.concat([montagud_node_synonyms, new_rows], ignore_index=True)


#     montagud_node_synonyms = montagud_node_synonyms.rename(columns={'HGNC symbols': 'Node_synonyms'})
    
#     montagud_node_synonyms['Node_synonyms'] = montagud_node_synonyms['Node_synonyms'].str.split(',')
#     montagud_node_synonyms = montagud_node_synonyms.explode('Node_synonyms')

#     synonyms_to_nodes_dict = montagud_node_synonyms.set_index('Node_synonyms')['Node'].to_dict()

#     return montagud_node_synonyms, synonyms_to_nodes_dict


def process_montagud_nodes_synonyms(montagud_synonyms_data):

    montagud_node_synonyms = montagud_synonyms_data.dropna(subset=['Node'])

     # Create multiple rows at once
    new_rows = pd.DataFrame({
        'Node': ['cFLAR', 'eEF2', 'eEF2K', 'Rheb'],  
        'HGNC symbols': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'], 
        'unique': ['CFLAR', 'EEF2', 'EEF2K', 'RHEB'] 
    })

    montagud_node_synonyms = pd.concat([montagud_node_synonyms, new_rows], ignore_index=True)


    montagud_node_synonyms = montagud_node_synonyms.rename(columns={'HGNC symbols': 'Node_synonyms'})
    
    montagud_node_synonyms['Node_synonyms'] = montagud_node_synonyms['Node_synonyms'].str.split(',')
    montagud_node_synonyms = montagud_node_synonyms.explode('Node_synonyms')

     # Strip whitespace from Node column
    montagud_node_synonyms['Node_synonyms'] = montagud_node_synonyms['Node_synonyms'].str.strip()
    montagud_node_synonyms['Node'] = montagud_node_synonyms['Node'].str.strip()


    montagud_node_synonyms = montagud_node_synonyms[montagud_node_synonyms['Node_synonyms'] != '']
    montagud_node_synonyms = montagud_node_synonyms[montagud_node_synonyms['Node'] != '']

    

    synonyms_to_nodes_dict = montagud_node_synonyms.set_index('Node_synonyms')['Node'].to_dict()

    return montagud_node_synonyms, synonyms_to_nodes_dict





def process_montagud_nodes(
    montagud_original_data_df, montagud_node_synonyms,
):

    # Create list of genes of interest (in Montagud data)
    # montagud_node_model are the nodes of the model 
    montagud_node_model = list(
        set(montagud_original_data_df["Target node"].tolist() + montagud_original_data_df["Source"].tolist())
    )
    montagud_node_model = [node for node in montagud_node_model if node != "0/1"]
    all_montagud_nodes = list(set(montagud_node_synonyms['Node_synonyms'])) + list(set(montagud_node_synonyms['Node'])) + list(set(montagud_node_model))



    return montagud_node_model, all_montagud_nodes


