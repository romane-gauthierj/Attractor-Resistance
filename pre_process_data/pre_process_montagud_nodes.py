

def process_montagud_nodes(montagud_data, node_to_remove=None):
    montagud_nodes = list(set(montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()))
    montagud_nodes = [node for node in montagud_nodes if node != '0/1']
    montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]
    if node_to_remove != None:
        montagud_nodes = [node for node in montagud_nodes if node not in node_to_remove]
        montagud_nodes = list(set(montagud_nodes))
    return montagud_nodes








# to_remove = ['FUSED_EVENT', 'AR_ERG']