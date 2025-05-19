def simulate_single_mutants(model, list_nodes=[], sign="BOTH", cmaboss=False):
    """
        Simulates a batch of single mutants and return an array of results

        :param model: the model on which to perform the simulations
        :param list_nodes: the node(s) which are to be mutated (default: all nodes)
        :param sign: which mutations to perform. "ON", "OFF", or "BOTH" (default)
    """ 
    
    if len(list_nodes) == 0:
        list_nodes = list(model.network.keys())
    
    list_single_mutants = []
    if (sign == "BOTH" or sign == "ON"):
        list_single_mutants += [(node, "ON") for node in list_nodes]
    if (sign == "BOTH" or sign == "OFF"):
        list_single_mutants += [(node, "OFF") for node in list_nodes]
    
    res = {}
    for single_mutant in list_single_mutants:
        t_model = model.copy()
        t_model.mutate(*single_mutant)
        res.update({single_mutant: t_model.run(cmaboss=cmaboss)})
        
    return res
    
    
def simulate_double_mutants(model, list_nodes, sign="BOTH", cmaboss=False):
    """
        Simulates a batch of double mutants and return an array of results

        :param model: the model on which to perform the simulations
        :param list_nodes: the node(s) which are to be mutated (default: all nodes)
        :param sign: which mutations to perform. "ON", "OFF", or "BOTH" (default)
    """ 
    
    if len(list_nodes) == 0:
        list_nodes = list(model.network.keys())
    
    list_single_mutants = []
    if (sign == "BOTH" or sign == "ON"):
        list_single_mutants += [(node, "ON") for node in list_nodes]
    if (sign == "BOTH" or sign == "OFF"):
        list_single_mutants += [(node, "OFF") for node in list_nodes]
    
    list_double_mutants = [(a, b) for idx, a in enumerate(list_single_mutants) for b in list_single_mutants[idx + 1:] if a[0] != b[0]]

    res = {}
    for double_mutant in list_double_mutants:
        t_model = model.copy()
        t_model.mutate(*(double_mutant[0]))
        t_model.mutate(*(double_mutant[1]))
        res.update({double_mutant: t_model.run(cmaboss=cmaboss)})
        
    return res
    
def filter_sensitivity(results, state=None, node=None, minimum=None, maximum=None):
    """
        Filter a list of results by state of nodes value

        :param results: the list of results to filter
        :param state: the state on which to apply the filter (default None)
        :param node: the state on which to apply the filter (default None)
        :param minumum: the minimal value of the node (default None)
        :param maximum: the maximal value of the node (default None)

        Example : 

        Filtering results showing more than 50% for Proliferation node
        >>> res_ensemble = filter_sensitivity(results, node='Proliferation', maximum=0.5)

        Filtering results showing more than 10% for Apoptosis -- NonACD state
        >>> res_ensemble = filter_sensitivity(results, state='Apoptosis -- NonACD', minimum=0.1)
    """ 
    ret_res = {}
    for (mutant, res) in results.items():
        
        if state is not None:
            t_res = res.get_last_states_probtraj()
            if state in t_res.columns:
                if minimum is not None and maximum is not None:
                    if t_res[state].values[0] >= minimum and t_res[state].values[0] <= maximum:
                        ret_res.update({mutant: res})
                elif minimum is not None and t_res[state].values[0] >= minimum:
                    ret_res.update({mutant: res})
                elif maximum is not None and t_res[state].values[0] <= maximum:
                    ret_res.update({mutant: res})
            elif maximum is not None and minimum is None:
                ret_res.update({mutant:res})

        elif node is not None:
            
            t_res = res.get_last_nodes_probtraj()

            if node in t_res.columns:
                if minimum is not None and maximum is not None:
                    if t_res[node].values[0] >= minimum and t_res[node].values[0] <= maximum:
                        ret_res.update({mutant: res})
                elif minimum is not None and t_res[node].values[0] >= minimum:
                    ret_res.update({mutant: res})
                elif maximum is not None and t_res[node].values[0] <= maximum:
                    ret_res.update({mutant: res})
            elif maximum is not None and minimum is None:
                ret_res.update({mutant: res})

    return ret_res
