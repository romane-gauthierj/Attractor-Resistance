"""
Class that contains the results of a MaBoSS simulation.
"""

import pandas
import tempfile
import os
from .results.baseresult import BaseResult
from contextlib import ExitStack
from time import time
class CMaBoSSResult2(BaseResult):

    def __init__(self, sim, workdir=None, only_final_state=False):

        BaseResult.__init__(self, "some_path")
        self.cmaboss_simulation = sim.cmaboss.MaBoSSSim(net=sim.cmaboss_net, cfg=sim.cmaboss_cfg)
        self.cmaboss_result = self.cmaboss_simulation.run(only_last_state=only_final_state)

        self.workdir = workdir
        self.only_final_state = only_final_state

        if workdir is not None:
            if not os.path.exists(workdir):
                os.makedirs(workdir)

            self.cmaboss_result.display_run(os.path.join(workdir, "%s_run.txt" % prefix))

            if self.only_final_state:
                self.cmaboss_result.display_final_states(os.path.join(workdir, "%s_finalprob.csv" % prefix))
            else:
                self.cmaboss_result.display_probtraj(os.path.join(workdir, "%s_probtraj.csv" % prefix))
                self.cmaboss_result.display_fp(os.path.join(workdir, "%s_fp.csv" % prefix))
                self.cmaboss_result.display_statdist(os.path.join(workdir, "%s_statdist.csv" % prefix))


    def get_last_states_probtraj(self, as_series=False):
        
        raw_res = self.cmaboss_result.get_last_probtraj()
        if not as_series:
            df = pandas.DataFrame(*raw_res)
            df.sort_index(axis=1, inplace=True)
        else:
            df = pandas.Series(raw_res[0][0], index=raw_res[2], name=raw_res[1][0])
            df.sort_index(inplace=True)

        return df

    def get_states_probtraj(self, prob_cutoff=None):
        if not self.only_final_state:

            raw_res = self.cmaboss_result.get_probtraj()
            df = pandas.DataFrame(*raw_res)
            df.sort_index(axis=1, inplace=True)

            if prob_cutoff is not None:
                maxs = df.max(axis=0)
                return df[maxs[maxs>prob_cutoff].index]

            return df

    def get_nodes_probtraj(self, prob_cutoff=None):
        if not self.only_final_state:
            raw_res = self.cmaboss_result.get_nodes_probtraj()
            df = pandas.DataFrame(*raw_res)
            df.sort_index(axis=1, inplace=True)

            if prob_cutoff is not None:
                maxs = df.max(axis=0)
                return df[maxs[maxs>prob_cutoff].index]

            return df

    def get_last_nodes_probtraj(self, as_series=False):
        raw_res = self.cmaboss_result.get_last_nodes_probtraj()
        if not as_series:
            df = pandas.DataFrame(*raw_res)
            df.sort_index(axis=0, inplace=True)
        else:
            df = pandas.Series(raw_res[0][0], index=raw_res[2], name=raw_res[1][0])
            df.sort_index(inplace=True)
        return df

    def get_observed_graph(self, prob_cutoff=None):
        raw_res = self.cmaboss_result.observed_graph()
        if raw_res is not None:
            graph = pandas.DataFrame(raw_res[0], columns=raw_res[1], index=raw_res[1])
            for state, values in graph.iterrows():
                if values.sum() > 0:
                    graph.loc[state, :] = values/values.sum()
            
            if prob_cutoff is not None:
                graph[graph < prob_cutoff] = 0
        
            return graph
        
    def get_observed_durations(self, prob_cutoff=None):
        raw_res = self.cmaboss_result.observed_durations()
        if raw_res is not None:
            graph = pandas.DataFrame(raw_res[0], columns=raw_res[1], index=raw_res[1])
            # for state, values in graph.iterrows():
                # if values.sum() > 0:
                #     graph.loc[state, :] = values/values.sum()
            
            if prob_cutoff is not None:
                graph[graph < prob_cutoff] = 0
        
            return graph
        
    def get_fptable(self):
        if not self.only_final_state:
            raw_res = self.cmaboss_result.get_fp_table()

            df = pandas.DataFrame(["#%d" % fp for fp in sorted(raw_res.keys())], columns=["FP"])

            df["Proba"] = [raw_res[fp][0] for fp in sorted(raw_res.keys())]
            df["State"] = [raw_res[fp][1] for fp in sorted(raw_res.keys())]

            for node in self.simul.network.keys():
                df[node] = [1 if node in raw_res[fp][1].split(" -- ") else 0 for fp in sorted(raw_res.keys())]

            return df


__all__ = ["CMaBoSSResult2"]
