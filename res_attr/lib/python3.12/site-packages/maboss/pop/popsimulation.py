"""
Class that contains the cMaBoSS Population simulation.
"""
# from .sbmlcmabossresult import SBMLCMaBoSSResult
# from sys import stdout
import os, sys
from .popresult import PopMaBoSSResult
from .cmabosspopresult import cMaBoSSPopMaBoSSResult
from .storedpopresult import StoredPopResult

class PopSimulation(object):

    def __init__(self, network=None, cfgs=None, maboss_model=None, cmaboss=None, cmaboss_sim=None):

        if maboss_model is not None:
            self.cmaboss = maboss_model.cmaboss
            self.cmaboss_sim = self.cmaboss.PopMaBoSSSim(network_str=maboss_model.str_bnd(), config_str=maboss_model.str_cfg())
        
        elif cmaboss is not None and cmaboss_sim is not None:
            self.cmaboss = cmaboss
            self.cmaboss_sim = cmaboss_sim
        else:
            nb_nodes = self.count_nodes(network)
            self.cmaboss = self.get_cmaboss(nb_nodes)

            if cfgs is None:
                self.cmaboss_sim = self.cmaboss.PopMaBoSSSim(network)
            else:
                if isinstance(cfgs, str):
                    self.cmaboss_sim = self.cmaboss.PopMaBoSSSim(network, config=cfgs)
                elif isinstance(cfgs, list):
                    self.cmaboss_sim = self.cmaboss.PopMaBoSSSim(network, configs=cfgs)

        self.param = self.cmaboss_sim.param
        self.network = self.cmaboss_sim.network
        
    def count_nodes(self, network):
        
        res = 0
        with open(network, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.lower().startswith("node"):
                    res += 1
        return res
    
    def get_cmaboss(self, nb_nodes):

        assert nb_nodes <= 1024, "Models with more than 1024 nodes are not compatible with this version of MaBoSS"

        if nb_nodes <= 64:
            return __import__("cmaboss")
        elif nb_nodes <= 128:
            return __import__("cmaboss_128n")
        elif nb_nodes <= 256:
            return __import__("cmaboss_256n")
        elif nb_nodes <= 512:
            return __import__("cmaboss_512n")
        else:
            return __import__("cmaboss_1024n")

    def update_parameters(self, **kwargs):
        self.cmaboss_sim.update_parameters(**kwargs)

    def set_custom_pop_output(self, formula):
        self.cmaboss_sim.set_custom_pop_output(formula)
    
    def copy(self):
        return PopSimulation(cmaboss=self.cmaboss, cmaboss_sim=self.cmaboss_sim.copy())
        
    # def get_nodes(self):
    #     return self.cmaboss_sim.get_nodes()

    def run(self, workdir=None, prefix="res", overwrite=False, hexfloat=True, cmaboss=True):
        if workdir is None or not os.path.exists(os.path.join(workdir, "%s_run.txt" % prefix)) or overwrite:
            return cMaBoSSPopMaBoSSResult(self, workdir, prefix, overwrite, hexfloat)
        else:
            return StoredPopResult(self, workdir, prefix, hexfloat)

    def str_bnd(self):
        return self.cmaboss_sim.str_bnd()
        
    def str_cfg(self):
        return self.cmaboss_sim.str_cfg()

    def print_bnd(self, out=sys.stdout):
        """Produce the content of the bnd file associated to the simulation."""
        print(self.str_bnd(), file=out)

    def print_cfg(self, out=sys.stdout):
        """Produce the content of the cfg file associated to the simulation."""
        print(self.str_cfg(), file=out)
__all__ = ["BNetSimulation"]
