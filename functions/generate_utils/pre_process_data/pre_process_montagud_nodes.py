def process_montagud_nodes(
    montagud_data, name_montagud_maps, nodes_to_add, nodes_to_remove
):
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
    for node in nodes_to_remove:
        if node in montagud_nodes:
            montagud_nodes.remove(node)

    return montagud_nodes

