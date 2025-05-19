"""
Class that contains the results of a PopMaBoSS simulation.
"""

import os

from .popresult import PopMaBoSSResult

class cMaBoSSPopMaBoSSResult(PopMaBoSSResult):

    def __init__(self, sim, workdir=None, prefix="res", overwrite=False, hexfloat=True):

        PopMaBoSSResult.__init__(self, sim)
        
        self._sim = sim        
        self.cmaboss_result = sim.cmaboss_sim.run()

        self.raw_states_probtraj = None
        self.raw_last_states_probtraj = None
        
        self.raw_simple_probtraj = None
        self.raw_simple_last_probtraj = None

        self.raw_custom_probtraj = None
        self.raw_custom_last_probtraj = None
        
        self.workdir = workdir

        if workdir is not None and (not os.path.exists(os.path.join(workdir, "%s_run.txt" % prefix)) or overwrite):

            self.cmaboss_result.display_run(os.path.join(workdir, "%s_run.txt" % prefix), hexfloat)
            self.cmaboss_result.display_probtraj(
                os.path.join(workdir, "%s_pop_probtraj.csv" % prefix), 
                os.path.join(workdir, "%s_simple_pop_probtraj.csv" % prefix),
                hexfloat
            )
            self.cmaboss_result.display_fp(os.path.join(workdir, "%s_fp.csv" % prefix), hexfloat)
            # self.cmaboss_result.display_statdist(os.path.join(workdir, "%s_statdist.csv" % prefix))

    def get_raw_last_states_probtraj(self):
        
        if self.raw_last_states_probtraj is None:
            self.raw_last_states_probtraj = self.cmaboss_result.get_last_probtraj()
        return self.raw_last_states_probtraj
    
    def get_raw_states_probtraj(self):
        if self.raw_states_probtraj is None:
            self.raw_states_probtraj = self.cmaboss_result.get_probtraj()
        return self.raw_states_probtraj

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
        
    def get_raw_custom_probtraj(self):
        if self.raw_custom_probtraj is None:
            self.raw_custom_probtraj = self.cmaboss_result.get_custom_probtraj()
        return self.raw_custom_probtraj
        
    def get_raw_custom_last_probtraj(self):
        if self.raw_custom_last_probtraj is None:
            self.raw_custom_last_probtraj = self.cmaboss_result.get_custom_last_probtraj()
        return self.raw_custom_last_probtraj
