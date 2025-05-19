"""
Class that contains the results of a MaBoSS simulation.
"""
import numpy as np
import os
from .popresult import PopMaBoSSResult
import itertools
class StoredPopResult(PopMaBoSSResult):
    
    def __init__(self, sim, workdir=None, prefix="res", hexfloat=True):

        PopMaBoSSResult.__init__(self, sim)

        self._sim = sim

        self.raw_states_probtraj = None
        self.raw_last_states_probtraj = None
        
        self.raw_simple_probtraj = None
        self.raw_simple_last_probtraj = None

        self.state_probtraj = None
        self.state_probtraj_errors = None


        self._workdir = workdir
        self._prefix = prefix
        
        
        self._raw_data = None
        self._first_state_index = None

        self._raw_states = None
        self._raw_probas = None
        self._raw_errors = None
        self._raw_entropy = None
        self._raw_last_data = None
        self._raw_custom_last_data = None

        self.indexes = None
        self.states = None
        self.nodes = None
        self.states_indexes = None
        self.nodes_indexes = None
        self._hexfloat = hexfloat
        
        
        
        self._simple_raw_data = None
        self._simple_first_state_index = None
        self._simple_raw_states = None
        self._simple_raw_probas = None
        self._simple_raw_errors = None
        self.simple_state_probtraj = None
        self.simple_state_probtraj_errors = None
        self.simple_nd_probtraj = None
        self.simple_nd_probtraj_error = None
        self.simple_indexes = None
        self.simple_states = None
        self.simple_states_indexes = None
        self.simple_nodes = None
        self.simple_nodes_indexes = None
        self._simple_raw_popsize = None
        self._simple_raw_popsize_errors = None
        self.simple_popsize = None
        self.simple_popsize_errors = None
        
        
    def get_raw_last_states_probtraj(self):
                
        data, first_col, hexfloat = self._get_raw_last_data()

        time = data[0]
        states = [s for s in data[first_col::3]]
        if hexfloat:
            probs = np.array([float.fromhex(v) for v in data[first_col+1::3]])
        else:
            probs = np.array([float(v) for v in data[first_col+1::3]])
        
        return ([probs], [time], states)
        
    def get_raw_states_probtraj(self):
    
        raw_states = self._get_raw_states()
        raw_probas = self._get_raw_probas()
        raw_errors = self._get_raw_errors()
        indexes, states = self._get_indexes()
        states_indexes = self._get_states_indexes()
        
        new_data = np.zeros((len(raw_probas), len(states)))
        for i, t_probas in enumerate(raw_probas):
            for j, proba in enumerate(t_probas):
                new_data[i, states_indexes[raw_states[i][j]]] = proba    
        
        new_errors = np.zeros((len(raw_errors), len(states)))
        for i, t_probas in enumerate(raw_errors):
            for j, proba in enumerate(t_probas):
                new_errors[i, states_indexes[raw_states[i][j]]] = proba    
            
        return new_data, indexes, states, new_errors

    def get_raw_states_probtraj_by_index(self, index):
        
        raw_line, first_index, hexfloat = self._get_raw_data_by_index(index)
        time = None
        states = []
        for i, token in enumerate(split_tab_gen(raw_line)):
            if i == 0:
                time = token
            elif i >= self._first_state_index and (i-first_index) % 3 == 0:
                states.append(token)
                
        states_indexes = {state:i for i, state in enumerate(states)}
                
        new_data = np.zeros((len(states)))
        raw_data = split_tab_gen(raw_line)
        state = None
        value = None
        for i, token in enumerate(raw_data):
            if i >= first_index and (i-first_index) % 3 == 0:
                state = token
            elif i >= first_index and (i-first_index) % 3 == 1:
                value = token
                if hexfloat:
                    new_data[states_indexes[state]] = float.fromhex(value)
                else:
                    new_data[states_indexes[state]] = float(value)
                    
        return [new_data], [time], states
    ########### Simple Last Probtraj

    def get_raw_simple_last_probtraj(self):
        if self.raw_simple_last_probtraj is None:
            self.raw_simple_last_probtraj = self.cmaboss_result.get_simple_last_probtraj()
        return self.raw_simple_last_probtraj

    ########### Simple Probtraj
    
    def get_raw_simple_probtraj(self):
        if self.raw_simple_probtraj is None:
            self.raw_simple_probtraj = self.cmaboss_result.get_simple_probtraj()
        return self.raw_simple_probtraj

    ########### Custom Probtraj
    
    def get_raw_custom_last_probtraj(self):
        data, first_col, hexfloat = self._get_raw_custom_last_data()

        time = data[0]
        states = [s for s in data[first_col::3]]
        if hexfloat:
            probs = np.array([float.fromhex(v) for v in data[first_col+1::3]])
        else:
            probs = np.array([float(v) for v in data[first_col+1::3]])
        
        return ([probs], [time], states)

    def _get_raw_custom_last_data(self):
        if self._raw_custom_last_data is None:
            with self._get_custom_probtraj_fd() as probtraj:

                if self._first_state_index is None:
                    first_line = probtraj.readline()
                    self._first_state_index = next(i for i, col in enumerate(first_line.strip("\n").split("\t")) if col == "State")

                last_line = probtraj.readlines()[-1]
                self._raw_custom_last_data = last_line.strip("\n").split("\t")
                self._hexfloat = self._raw_custom_last_data[self._first_state_index+1].startswith("0x")
        
        return self._raw_custom_last_data, self._first_state_index, self._hexfloat

    

    def get_probtraj_file(self):
        return os.path.join(self._workdir, "%s_pop_probtraj.csv" % self._prefix)

    def get_simple_probtraj_file(self):
        return os.path.join(self._workdir, "%s_simple_pop_probtraj.csv" % self._prefix)
    
    def get_custom_probtraj_file(self):
        return os.path.join(self._workdir, "%s_custom_pop_probtraj.csv" % self._prefix)
    
    def _get_probtraj_fd(self):
        return open(self.get_probtraj_file(), 'r')

    def _get_custom_probtraj_fd(self):
        return open(self.get_custom_probtraj_file(), 'r')

    def _get_raw_data(self):
    
        if self._raw_data is None:
            with self._get_probtraj_fd() as probtraj:

                raw_lines = probtraj.readlines()

                if self._first_state_index is None:                    
                    self._first_state_index = next(i for i, col in enumerate(raw_lines[0].strip("\n").split("\t")) if col == "State")

                self._raw_data = [line.strip("\n").split("\t") for line in raw_lines[1:]]
                self._hexfloat = self._raw_data[0][self._first_state_index+1].startswith("0x")

        return self._raw_data, self._first_state_index, self._hexfloat
        
    def _get_raw_last_data(self):
    
        if self._raw_last_data is None:
            with self._get_probtraj_fd() as probtraj:

                if self._first_state_index is None:
                    first_line = probtraj.readline()
                    self._first_state_index = next(i for i, col in enumerate(first_line.strip("\n").split("\t")) if col == "State")

                last_line = probtraj.readlines()[-1]
                self._raw_last_data = last_line.strip("\n").split("\t")
                self._hexfloat = self._raw_last_data[self._first_state_index+1].startswith("0x")
        
        return self._raw_last_data, self._first_state_index, self._hexfloat

    def _get_raw_data_by_index(self, index):
    
        with self._get_probtraj_fd() as probtraj:

            for i, line in enumerate(probtraj):
                if i == 0 and self._first_state_index is None:
                    self._get_first_state_index(line)
                if i == index+1:
                    raw_data = split_tab_gen(line)
                    
                    first_value = next(itertools.islice(raw_data, self._first_state_index+1, None))
                    self._hexfloat = first_value.startswith("0x")
        
                    return line, self._first_state_index, self._hexfloat


    def _get_raw_states(self):
        if self._raw_states is None:
            data, first_state_index, _ = self._get_raw_data()
            self._raw_states = [[s for s in t_data[first_state_index::3]] for t_data in data]
        return self._raw_states

    def _get_raw_probas(self):
        if self._raw_probas is None:
            data, first_state_index, hexfloat = self._get_raw_data()
            if hexfloat:
                self._raw_probas = [[float.fromhex(p) for p in t_data[first_state_index+1::3]] for t_data in data]
            else:
                self._raw_probas = [[np.float64(p) for p in t_data[first_state_index+1::3]] for t_data in data]
        return self._raw_probas

    def _get_raw_errors(self):
        if self._raw_errors is None:
            data, first_state_index, hexfloat = self._get_raw_data()
            if hexfloat:
                self._raw_errors = [[float.fromhex(p) for p in t_data[first_state_index+2::3]] for t_data in data]
            else:
                self._raw_errors = [[np.float64(p) for p in t_data[first_state_index+2::3]] for t_data in data]
        return self._raw_errors
        
    def _get_raw_entropy(self):
        if self._raw_entropy is None:
            data, _, hexfloat = self._get_raw_data()
            if hexfloat:
                self._raw_entropy = [[float.fromhex(p) for p in t_data[1:4]] for t_data in data]
            else:
                self._raw_entropy = [[np.float64(p) for p in t_data[1:4]] for t_data in data]
        return self._raw_entropy

    def _get_indexes(self):

        if self.indexes is None:
            data, _, hexfloat = self._get_raw_data()
            self.indexes = [float(t_data[0]) for t_data in data] 
            
        if self.states is None:
            self.states = set()
            for t_states in self._get_raw_states():
                self.states.update(t_states)
            self.states = list(self.states)

        return self.indexes, self.states
    
    def _get_first_state_index(self, line):
        if self._first_state_index is None:
            self._first_state_index = next(i for i, col in enumerate(split_tab_gen(line)) if col == "State")
        return self._first_state_index
    
    def _get_indexes_nocache(self):

        states = set()
        indexes = []
        with self._get_probtraj_fd() as probtraj:

            for line in probtraj.readlines()[1:]:
                raw_data = line.strip("\n").split("\t")
                
                indexes.append(float(raw_data[0]))
                states.update([s for s in raw_data[first_state_index::3]])
                
        return indexes, list(states)
    
    def _get_nodes(self):
        
        if self.nodes is None:
            self.nodes = set()
            for t_states in self._get_raw_states():
                for state in t_states:
                    if state != "<nil>":
                        self.nodes.update([node for node in state.split(" -- ")])
            self.nodes = list(self.nodes)
        
        return self.nodes

    def _get_states_indexes(self):
        if self.states_indexes is None:
            _, states = self._get_indexes()
            self.states_indexes = {state:index for index, state in enumerate(states)}
        return self.states_indexes

    def _get_nodes_indexes(self):
        if self.nodes_indexes is None:
            nodes = self.output_nodes if self.output_nodes is not None else self._get_nodes()
            self.nodes_indexes = {node:index for index, node in enumerate(nodes)}
        return self.nodes_indexes

        
        
        
        
        
        
        
        
    def _get_simple_probtraj_fd(self):
        return open(self.get_simple_probtraj_file(), 'r')

    def _get_simple_raw_data(self):
    
        if self._simple_raw_data is None:
            with self._get_simple_probtraj_fd() as probtraj:

                raw_lines = probtraj.readlines()

                if self._simple_first_state_index is None:                    
                    self._simple_first_state_index = next(i for i, col in enumerate(raw_lines[0].strip("\n").split("\t")) if col == "State")
                
                self._simple_raw_data = [line.strip("\n").split("\t") for line in raw_lines[1:]]
                self._hexfloat = self._simple_raw_data[0][self._simple_first_state_index+1].startswith("0x")

        return self._simple_raw_data, self._simple_first_state_index, self._hexfloat
        
    def _get_simple_raw_popsize(self):
        if self._simple_raw_popsize is None:
            data, _, hexfloat = self._get_simple_raw_data()
            if hexfloat:
                self._simple_raw_popsize = [float.fromhex(t_data[5]) for t_data in data]
            else:
                self._simple_raw_popsize = [np.float64(t_data[5]) for t_data in data]
    
        return self._simple_raw_popsize
    
    def _get_simple_raw_popsize_errors(self):
        if self._simple_raw_popsize_errors is None:
            data, _, hexfloat = self._get_simple_raw_data()
            if hexfloat:
                self._simple_raw_popsize_errors = [float.fromhex(t_data[6]) for t_data in data]
            else:
                self._simple_raw_popsize_errors = [np.float64(t_data[6]) for t_data in data]
    
        return self._simple_raw_popsize_errors
        
    def _get_simple_raw_states(self):
        if self._simple_raw_states is None:
            data, first_state_index, _ = self._get_simple_raw_data()
            self._simple_raw_states = [[s for s in t_data[first_state_index::3]] for t_data in data]
        return self._simple_raw_states

    def _get_simple_raw_probas(self):
        if self._simple_raw_probas is None:
            data, first_state_index, hexfloat = self._get_simple_raw_data()
            if hexfloat:
                self._simple_raw_probas = [[float.fromhex(p) for p in t_data[first_state_index+1::3]] for t_data in data]
            else:
                self._simple_raw_probas = [[np.float64(p) for p in t_data[first_state_index+1::3]] for t_data in data]
        return self._simple_raw_probas

    def _get_simple_raw_errors(self):
        if self._simple_raw_errors is None:
            data, first_state_index, hexfloat = self._get_simple_raw_data()
            if hexfloat:
                self._simple_raw_errors = [[float.fromhex(p) for p in t_data[first_state_index+2::3]] for t_data in data]
            else:
                self._simple_raw_errors = [[np.float64(p) for p in t_data[first_state_index+2::3]] for t_data in data]
        return self._simple_raw_errors

    def _get_simple_indexes(self):

        if self.simple_indexes is None:
            data, _, _ = self._get_simple_raw_data()
            self.simple_indexes = [float(t_data[0]) for t_data in data] 
            
        if self.simple_states is None:
            self.simple_states = set()
            for t_states in self._get_simple_raw_states():
                self.simple_states.update(t_states)
            self.simple_states = list(self.simple_states)

        return self.simple_indexes, self.simple_states
    
    def _get_simple_states_indexes(self):
        if self.simple_states_indexes is None:
            _, states = self._get_simple_indexes()
            self.simple_states_indexes = {state:index for index, state in enumerate(states)}
        return self.simple_states_indexes
    
    def _get_simple_nodes(self):
        
        if self.simple_nodes is None:
            self.simple_nodes = set()
            for t_states in self._get_simple_raw_states():
                for state in t_states:
                    if state != "<nil>":
                        self.simple_nodes.update([node for node in state.split(" -- ")])
            self.simple_nodes = list(self.simple_nodes)
        
        return self.simple_nodes

    def _get_simple_nodes_indexes(self):
        if self.simple_nodes_indexes is None:
            nodes = self.output_nodes if self.output_nodes is not None else self._get_simple_nodes()
            self.simple_nodes_indexes = {node:index for index, node in enumerate(nodes)}
        return self.simple_nodes_indexes
        
    def get_raw_simple_probtraj(self, prob_cutoff=None):
        """
            Returns the simplified population state probability vs time, as a pandas dataframe.

            :param float prob_cutoff: returns only the states with proba > cutoff
        """

        raw_states = self._get_simple_raw_states()
        raw_probas = self._get_simple_raw_probas()
        raw_errors = self._get_simple_raw_errors()
        raw_popsizes = self._get_simple_raw_popsize()
        raw_popsizes_errors = self._get_simple_raw_popsize_errors()
        
        indexes, states = self._get_simple_indexes()
        states_indexes = self._get_simple_states_indexes()
        new_data = np.zeros((len(raw_probas), len(states)+1))
        for i, t_probas in enumerate(raw_probas):
            new_data[i, 0] = raw_popsizes[i]
            for j, proba in enumerate(t_probas):
                new_data[i, states_indexes[raw_states[i][j]]+1] = proba    
        
        new_error = np.zeros((len(raw_errors), len(states)+1))
        for i, t_probas in enumerate(raw_errors):
            new_error[i, 0] = raw_popsizes_errors[i]
            for j, proba in enumerate(t_probas):
                new_error[i, states_indexes[raw_states[i][j]]+1] = proba    
                
        return new_data, indexes, ["Population"] + states, new_error
    
    # def get_simple_states_probtraj_errors(self, prob_cutoff=None):
    #     """
    #         Returns the simplified population state probability vs time, as a pandas dataframe.

    #         :param float prob_cutoff: returns only the states with proba > cutoff
    #     """
    #     if self.simple_state_probtraj_errors is None:

    #         raw_states = self._get_simple_raw_states()
    #         raw_errors = self._get_simple_raw_errors()
            
    #         indexes, states = self._get_simple_indexes()
    #         states_indexes = self._get_simple_states_indexes()
    #         new_data = np.zeros((len(raw_errors), len(states)))
    #         for i, t_probas in enumerate(raw_errors):
    #             for j, proba in enumerate(t_probas):
    #                 new_data[i, states_indexes[raw_states[i][j]]] = proba    
            
    #         self.simple_state_probtraj_errors = pd.DataFrame(
    #             data=new_data,
    #             columns=states,
    #             index=indexes

    #         )
    #         self.simple_state_probtraj_errors.sort_index(axis=1, inplace=True)

    #     if prob_cutoff is not None:
    #         maxs = self.simple_state_probtraj_errors.max(axis=0)
    #         return self.simple_state_probtraj_errors[maxs[maxs>prob_cutoff].index]

    #     return self.simple_state_probtraj_errors


    # def get_simple_nodes_probtraj(self, nodes=None, prob_cutoff=None):
    #     """
    #         Returns the node probability vs time, as a pandas dataframe.

    #         :param float prob_cutoff: returns only the nodes with proba > cutoff
    #     """
    #     if self.simple_nd_probtraj is None:

    #         raw_states = self._get_simple_raw_states()
    #         raw_probas = self._get_simple_raw_probas()
    #         indexes, states = self._get_simple_indexes()
    #         nodes = self.output_nodes if self.output_nodes is not None else self._get_simple_nodes()
    #         nodes_indexes = self._get_simple_nodes_indexes()

    #         new_probs = np.zeros((len(indexes), len(nodes)))
    #         for i, t_probas in enumerate(raw_probas):
    #             for j, proba in enumerate(t_probas):
    #                 if raw_states[i][j] != "<nil>":
    #                     for node in raw_states[i][j].split(" -- "):
    #                         new_probs[i, nodes_indexes[node]] += proba

    #         self.simple_nd_probtraj = pd.DataFrame(new_probs, columns=nodes, index=indexes)
    #         self.simple_nd_probtraj.sort_index(axis=1, inplace=True)

    #     if prob_cutoff is not None:
    #         maxs = self.simple_nd_probtraj.max(axis=0)
    #         return self.simple_nd_probtraj[maxs[maxs>prob_cutoff].index]

    #     if nodes is not None:
    #         return self.simple_nd_probtraj[nodes]

    #     return self.simple_nd_probtraj

    # def get_simple_nodes_probtraj_error(self):
    #     """
    #         Returns the node probability error vs time, as a pandas dataframe.
    #     """
    #     if self.simple_nd_probtraj_error is None:

    #         raw_states = self._get_simple_raw_states()
    #         raw_errors = self._get_simple_raw_errors()
    #         indexes, states = self._get_simple_indexes()
    #         nodes = self.output_nodes if self.output_nodes is not None else self._get_simple_nodes()
    #         nodes_indexes = self._get_simple_nodes_indexes()
            
    #         new_errors = np.zeros((len(indexes), len(nodes)))
    #         for i, t_raw_errors in enumerate(raw_errors):
    #             for j, error in enumerate(t_raw_errors):
    #                 if raw_states[i][j] != "<nil>":
    #                     for node in raw_states[i][j].split(" -- "):
    #                         new_errors[i, nodes_indexes[node]] += error

    #         self.simple_nd_probtraj_error = pd.DataFrame(new_errors, columns=nodes, index=indexes)
    #         self.simple_nd_probtraj_error.sort_index(axis=1, inplace=True)

    #     return self.simple_nd_probtraj_error


    # def get_simple_popsize(self):
    #     if self.simple_popsize is None:
    #         indexes, _ = self._get_simple_indexes()
    #         popsizes = self._get_simple_raw_popsize()
            
    #         self.simple_popsize = pd.Series(popsizes, index=indexes, name="Popsize")
    #     return self.simple_popsize
    
    # def get_simple_popsize_errors(self):
    #     if self.simple_popsize_errors is None:
    #         indexes, _ = self._get_simple_indexes()
    #         popsizes_errors = self._get_simple_raw_popsize_errors()
            
    #         self.simple_popsize_errors = pd.Series(popsizes_errors, index=indexes, name="Popsize")
    #     return self.simple_popsize_errors
           
    # def get_simple_states_popsize(self):
    #     return self.get_simple_states_probtraj().multiply(self.get_simple_popsize(), axis=0)
       
       
    # def get_simple_nodes_popsize(self, nodes=None):
    #     return self.get_simple_nodes_probtraj(nodes).multiply(self.get_simple_popsize(), axis=0)
       
def split_tab_gen(string):
    start = 0
    for i, c in enumerate(string):
        if c == "\t":
            yield string[start:i]
            start = i+1
    yield string[start:]