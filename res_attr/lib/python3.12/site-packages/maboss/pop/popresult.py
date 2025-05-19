"""
Class that contains the results of a PopMaBoSS simulation.
"""

import pandas, numpy, os

class PopMaBoSSResult:
            
    def __init__(self, sim):
        self._sim = sim
        
    def get_last_states_probtraj(self):
        
        raw_res = self.get_raw_last_states_probtraj()
        df = pandas.Series(
            raw_res[0][0], 
            index=raw_res[2], 
            name=raw_res[1][0]
        )
        df.sort_index(inplace=True)

        return df

    def get_states_probtraj(self, prob_cutoff=None):

        raw_res = self.get_raw_states_probtraj()
        df = pandas.DataFrame(*(raw_res[0:3]))
        # df.sort_index(axis=1, inplace=True)

        if prob_cutoff is not None:
            maxs = df.max(axis=0)
            return df[maxs[maxs>prob_cutoff].index]

        return df

    def get_states_probtraj_by_index(self, index):
        
        raw_res = self.get_raw_states_probtraj_by_index(index)
        df = pandas.Series(
            raw_res[0][0], 
            index=raw_res[2], 
            name=raw_res[1][0]
        )
        df.sort_index(inplace=True)

        return df

    ########### Simple Last Probtraj

    def get_last_simple_states_probtraj(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.Series(raw_res[0][1:], index=raw_res[1], name=raw_res[2][1:])
        df.sort_index(inplace=True)
        return df
    
    def get_last_simple_states_probtraj_errors(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.Series(raw_res[3][1:], index=raw_res[1], name=raw_res[2][1:])
        df.sort_index(inplace=True)
        return df

    ########### Simple Probtraj
    
    def get_simple_states_probtraj(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.DataFrame(raw_res[0][:,1:], index=raw_res[1], columns=raw_res[2][1:])
        df.sort_index(axis=1, inplace=True)
        return df

    def get_simple_states_probtraj_errors(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.DataFrame(raw_res[3][:,1:], index=raw_res[1], columns=raw_res[2][1:])
        df.sort_index(axis=1, inplace=True)
        return df

    def get_custom_states_probtraj(self):
        raw_res = self.get_raw_custom_probtraj()
        df = pandas.DataFrame(raw_res[0][:,], index=raw_res[1], columns=[int(x) for x in raw_res[2]])
        df.sort_index(axis=1, inplace=True)
        df.sort_index(axis=0, inplace=True)
        return df
        
    def get_custom_states_probtraj_errors(self):
        raw_res = self.get_raw_custom_probtraj()
        df = pandas.DataFrame(raw_res[3][:,], index=raw_res[1], columns=[int(x) for x in raw_res[2]])
        df.sort_index(axis=1, inplace=True)
        df.sort_index(axis=0, inplace=True)
        
        return df
        
    def get_last_custom_states_probtraj(self, prob_cutoff=None, rescale=True):
        raw_res = self.get_raw_custom_last_probtraj()
        df = pandas.Series(raw_res[0][0], index=[int(x) for x in raw_res[2]], name=raw_res[1][0])
        
        if rescale:
            df /= df.sum()
        
        if prob_cutoff is not None:
            df = df[df > prob_cutoff]
        
        df.sort_index(inplace=True)
        
        return df
    
    
    ########### Popsize
    
    def get_simple_popsize(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.Series(raw_res[0][:,0], index=raw_res[1], name=raw_res[2][0]) 
        return df
    
    def get_simple_popsize_errors(self):
        raw_res = self.get_raw_simple_probtraj()
        df = pandas.Series(raw_res[3][:,0], index=raw_res[1], name=raw_res[2][0]) 
        return df
    
    def get_simple_states_popsize(self):
        return self.get_simple_states_probtraj().multiply(self.get_simple_popsize(), axis=0)
    
    # def get_simple_states_popsize_errors(self):
    #     return self.get_simple_states_probtraj_errors().multiply(self.get_simple_popsize_errors(), axis=0)
    
    ########### Simple nodes probtraj
    
    def get_simple_nodes_probtraj(self):
        
        raw_probas, indexes, states, _ = self.get_raw_simple_probtraj()
        nodes = self._sim.get_nodes()
        nodes_indexes = {node:index for index, node in enumerate(nodes)}
        
        nodes_probas = numpy.zeros((len(indexes), len(nodes)))
        for i in range(len(indexes)):
            for j, state in enumerate(states[1:]):
                if state != "<nil>":
                    for node in state.split(" -- "):
                        nodes_probas[i, nodes_indexes[node]] += raw_probas[i][j+1]
    
        df = pandas.DataFrame(nodes_probas, index=indexes, columns=nodes)
        return df
    
    def get_simple_nodes_errors(self):
        
        _, indexes, states, raw_errors = self.get_raw_simple_probtraj()
        nodes = self._sim.get_nodes()
        nodes_indexes = {node:index for index, node in enumerate(nodes)}
        
        nodes_errors = numpy.zeros((len(indexes), len(nodes)))
        for i in range(len(indexes)):
            for j, state in enumerate(states[1:]):
                if state != "<nil>":
                    for node in state.split(" -- "):
                        nodes_errors[i, nodes_indexes[node]] += raw_errors[i][j+1]
    
        df = pandas.DataFrame(nodes_errors, index=indexes, columns=nodes)
        return df
    
    
    ########### Nodes Probtraj
 
    def get_nodes_probtraj(self, nodes=[], prob_cutoff=None):
        
        # if len(nodes) == 0:
        #     _nodes = self.sim.get_nodes()
            
        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
                
        self._popnodes = set()  #{"{%s:%s}" % (node,0) for node in _nodes}
        for t_state in states:
            for subpopstates in t_state[1:-1].split(","):
                state, pop = subpopstates[1:-1].split(":")
                if state != "<nil>":
                    for node in state.split(" -- "):
                        self._popnodes.add("{%s:%s}" % (node, pop))
                        
                        
        def sort_popnodes(popnode):
            node, pop = popnode[1:-1].split(":")
            return (node, int(pop))
        
        self._popnodes = sorted(list(self._popnodes), key=lambda key: sort_popnodes(key))
        self._popnodes_index = {popnode:index for index, popnode in enumerate(self._popnodes)}                
        
        self._popnodes_by_nodes = {node:[] for node in self._sim.get_nodes()}
        for popnode in self._popnodes:
            node, pop = popnode[1:-1].split(":")
            self._popnodes_by_nodes[node].append(popnode)
                
        new_probs = numpy.zeros((len(indexes), len(self._popnodes)))
        for i in range(len(indexes)):
            for j, state in enumerate(states):
                for subpopstates in state[1:-1].split(","):
                    state, pop = subpopstates[1:-1].split(":")
                    if state != "<nil>":
                        for node in state.split(" -- "):
                            new_probs[i, self._popnodes_index["{%s:%s}" % (node, pop)]] += raw_probas[i][j]
        
        self.nd_probtraj = pandas.DataFrame(new_probs, columns=self._popnodes, index=indexes)
        
        _list_popnodes = []
        if len(nodes) > 0:
            for node in nodes:
                _list_popnodes += self._popnodes_by_nodes[node]
        
        if prob_cutoff is not None:
            if len(_list_popnodes) > 0:
                maxs = self.nd_probtraj[_list_popnodes].max(axis=0)
            else:
                maxs = self.nd_probtraj.max(axis=0)
             
            return self.nd_probtraj[maxs[maxs>prob_cutoff].index]

        if len(_list_popnodes) > 0:
            return self.nd_probtraj[_list_popnodes]
        
        return self.nd_probtraj

    # def get_last_nodes_probtraj(self):
    #     raw_res = self.cmaboss_result.get_last_nodes_probtraj()
        
    #     df = pandas.Series(raw_res[0][0], index=raw_res[2], name=raw_res[1][0])
    #     df.sort_index(inplace=True)
    #     return df


    # def get_fptable(self):
    #     raw_res = self.cmaboss_result.get_fp_table()

    #     df = pandas.DataFrame(["#%d" % fp for fp in sorted(raw_res.keys())], columns=["FP"])

    #     df["Proba"] = [raw_res[fp][0] for fp in sorted(raw_res.keys())]
    #     df["State"] = [raw_res[fp][1] for fp in sorted(raw_res.keys())]

    #     for node in self.simul.network.keys():
    #         df[node] = [1 if node in raw_res[fp][1].split(" -- ") else 0 for fp in sorted(raw_res.keys())]

    #     return df

    def get_state_dist_by_index(self, index, fun):
        dist_state = {}
        for raw_popstates, proba in self.get_states_probtraj_by_index(index).items():
            value = fun(raw_popstates)
            if value is not None:
                if value in dist_state.keys():
                    dist_state[value] += proba
                else:
                    dist_state[value] = proba
            
        serie = pandas.Series(dist_state).sort_index()
        serie *= 1/serie.sum()
        return serie

    def get_last_state_dist(self, networkstate=None, rescale=True):
        if networkstate is None:
            return self.get_last_states_probtraj()
            
        dist_state = {}
        for raw_popstates, value in self.get_last_states_probtraj().items():
            popstates = parse_pop_state(raw_popstates)
            if networkstate in popstates.keys():
                if popstates[networkstate] in dist_state.keys():
                    dist_state[popstates[networkstate]] += value
                else:
                    dist_state[popstates[networkstate]] = value
                    
        df = pandas.Series(dist_state).sort_index()
        if rescale:
            df /= df.sum()
        return df


    def plot_last_state_dist(self, networkstate=None, **args):
        if networkstate is None:
            self.get_last_state_dist().plot.pie(autopct='%1.1f%%', legend=False)
        else:
            self.get_last_state_dist(networkstate).plot.bar(**args)
        
    # def get_node_dist(self, node):
    #     dist_node = {}
    #     for raw_popnode, value in serie.items():
    #         state, pop = raw_popnode[1:-1].split(":")
    #         dist_node.update({int(pop): value})
    #     return pandas.Series(dist_node).sort_index()
    def get_nb_dists(self):

        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
        
        sizes = set()
        for state in states:
            popstates = parse_pop_state(state)
            nb = sum(popstates.values())
            sizes.add(nb)
            
        sizes = sorted(list(sizes))
        sizes_indexes = {size:index for index, size in enumerate(sizes)}
        
        raw_nb_dists = numpy.zeros((len(indexes), len(sizes)))
        for i in range(len(indexes)):
            for j, state in enumerate(states):
                popstates = parse_pop_state(state)
                nb = sum(popstates.values())
                raw_nb_dists[i, sizes_indexes[nb]] += raw_probas[i][j]
                
        return pandas.DataFrame(raw_nb_dists, columns=sizes, index=indexes)
    
    def get_last_nb_dist(self):
        dist_state = {}
        for raw_popstates, value in self.get_last_states_probtraj().items():
            if value > 0:
                popstates = parse_pop_state(raw_popstates)
                nb = sum(popstates.values())
                if nb in dist_state.keys():
                    dist_state[nb] += value
                else:
                    dist_state[nb] = value
        return pandas.Series(dist_state).sort_index()

    def plot_last_nb_dist(self, **xargs):
        self.get_last_nb_dist().plot.bar(**xargs)
    
    def get_last_log_nb_dist(self, num=20, vmin=None, vmax=None, lspace=None):
        nbdist = self.get_last_nb_dist()
        ldist = numpy.log(numpy.array(nbdist.index))
        values = {}
        if lspace is None:
            if vmin is None:
                lmin = numpy.floor(ldist[0])
            else:
                lmin = numpy.log(vmin)
                
            if vmax is None:
                lmax = numpy.ceil(ldist[-1])
            else:
                lmax = numpy.log(vmax)
                
            lspace = numpy.linspace(lmin, lmax, num)
        for i, lognb in enumerate(lspace):
            if i < (len(lspace)-1):
                where = numpy.where((ldist >= lognb) & (ldist < lspace[i+1]))
                proba = nbdist.iloc[where].sum()
                values[numpy.exp(lognb)] = proba
        return pandas.Series(values)

    def plot_last_log_nb_dist(self, num=20, vmin=None, vmax=None, lspace=None, **xargs):
        self.get_last_log_nb_dist(num, vmin, vmax, lspace).plot.bar(**xargs)

    ########### Activity ratios

    def get_activity_ratio_expected(self, node):
        
        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
        activity_ratio = []
        for i in range(len(indexes)):
            activity = 0
            for j, state in enumerate(states):
                popstates = parse_pop_state(state)
                nb = sum(popstates.values())
                
                for state, pop in popstates.items():
                    if state != "<nil>" and node in state.split(" -- "):                 
                        activity += raw_probas[i][j]*pop/nb
                        
            activity_ratio.append(activity)
        return pandas.Series(activity_ratio, index=indexes, name=node)
        
    def get_activity_ratio_stdev(self, node):
        
        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
        activity_ratio = self.get_activity_ratio_expected(node)
        
        stdevs = []
        for i in range(len(indexes)):
            stdev = 0
            for j, state in enumerate(states):
                    popstates = parse_pop_state(state)
                    nb = sum(popstates.values())
                    t_proba = 0
                    t_stdev = 0
                    for state, pop in popstates.items():
                        if state != "<nil>" and node in state.split(" -- "):                 
                            t_proba += raw_probas[i][j]
                            t_stdev += pop/nb
                    stdev += ((t_stdev-activity_ratio.iloc[i])**2)*t_proba            
            stdevs.append((stdev)**0.5)
        return pandas.Series(stdevs, index=indexes, name=node)
    
    def get_activity_ratio_nodes(self, node0, node1, time=None, nbins=20):
        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
        
        if time is None:
            raw_probas = raw_probas[-1]
            indexes = indexes[-1]
        else:
            raw_probas = raw_probas[indexes.index(time)]
            indexes = time
            
        space0 = numpy.linspace(0, 1, nbins+1)
        space1 = numpy.linspace(0, 1, nbins+1)
        probas = numpy.zeros((len(space0), len(space1)))
        for j, state in enumerate(states):
            popstates = parse_pop_state(state)
            nb = sum(popstates.values())
            pop0 = 0
            pop1 = 0
            for state, pop in popstates.items():
                if state != "<nil>":
                    popnb = pop/nb
                    if node0 in state.split(" -- "):                 
                        pop0 += popnb
                    if node1 in state.split(" -- "):
                        pop1 += popnb
            pos0 =  numpy.floor(pop0*nbins).astype(int)
            pos1 =  numpy.floor(pop1*nbins).astype(int)
            probas[pos0, pos1] += raw_probas[j]

        return pandas.DataFrame(probas, index=["%.2g" % val for val in space1], columns=["%.2g" % val for val in space0])
    
    def get_activity_ratio_states(self, state0, state1, time=None, nbins=20):
        raw_probas, indexes, states, _ = self.get_raw_states_probtraj()
        
        if time is None:
            raw_probas = raw_probas[-1]
            indexes = indexes[-1]
        else:
            raw_probas = raw_probas[indexes.index(time)]
            indexes = time
            
        space0 = numpy.linspace(0, 1, nbins+1)
        space1 = numpy.linspace(0, 1, nbins+1)
        probas = numpy.zeros((len(space0), len(space1)))
        for j, state in enumerate(states):
            popstates = parse_pop_state(state)
            nb = sum(popstates.values())
            pop0 = 0
            pop1 = 0
            for state, pop in popstates.items():
                if state == state0:
                    pop0 += pop/nb
                if state == state1:
                    pop1 += pop/nb
            pos0 =  numpy.floor(pop0*nbins).astype(int)
            pos1 =  numpy.floor(pop1*nbins).astype(int)
            probas[pos0, pos1] += raw_probas[j]

        return pandas.DataFrame(probas, index=["%.2g" % val for val in space1], columns=["%.2g" % val for val in space0])
    
def parse_pop_state(cols):    
    pops = {}
    if len(cols) > 2:
        for pop in cols[1:-1].split(","):
            name, value = pop[1:-1].split(":")
            pops[name] = int(value)
    return pops

__all__ = ["PopMaBoSSResult"]
