# def process_montagud_nodes(montagud_data, node_to_remove=None):
#     montagud_nodes = list(set(montagud_data['Target node'].tolist() + montagud_data['Source'].tolist()))
#     montagud_nodes = [node for node in montagud_nodes if node != '0/1']
#     montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]
#     if node_to_remove != None:
#         montagud_nodes = [node for node in montagud_nodes if node not in node_to_remove]
#         montagud_nodes = list(set(montagud_nodes))
#     return montagud_nodes


def process_montagud_nodes(montagud_data, name_montagud_maps, nodes_to_add):
    # transform nodes_to_add to list if single element
    if isinstance(nodes_to_add, str):
        nodes_to_add = [nodes_to_add]

    # Create list of genes of interest (in Montagud data)
    montagud_nodes = list(
        set(montagud_data["Target node"].tolist() + montagud_data["Source"].tolist())
    )
    montagud_nodes = [node for node in montagud_nodes if node != "0/1"]
    montagud_nodes = [node.upper() for node in montagud_nodes if isinstance(node, str)]

    # Apply mapping
    montagud_nodes = [name_montagud_maps.get(x, x) for x in montagud_nodes]
    # Add any additional nodes
    montagud_nodes = list(set(montagud_nodes + nodes_to_add))

    return montagud_nodes


# to_remove = ['FUSED_EVENT', 'AR_ERG']
