"""
Class that contains the cMaBoSS simulation.
"""
from .cmabossresult2 import CMaBoSSResult2
from sys import stdout

class CMaBoSSSimulation(object):

    def __init__(self, bnd=None, cfgs=None, cmaboss=None, cmaboss_sim=None):
        
        if cmaboss is not None and cmaboss_sim is not None:
            self.cmaboss = cmaboss
            self.cmaboss_sim = cmaboss_sim
            self.network = self.cmaboss_sim.network
        else:
            self.bnd = bnd
            self.cfgs = cfgs

            self.nb_nodes = self.count_nodes()
            self.cmaboss = self.get_cmaboss()

            self.cmaboss_net = self.cmaboss.MaBoSSNet(self.bnd)
            self.cmaboss_cfg = self.cmaboss.MaBoSSCfg(self.cmaboss_net, self.cfgs)
            
            self.cmaboss_sim = self.cmaboss.MaBoSSSim(net=self.cmaboss_net, cfg=self.cmaboss_cfg)
            self.network = self.cmaboss_net    
            
        self.param = self.cmaboss_sim.param
        

    def count_nodes(self):
        
        res = 0
        with open(self.bnd, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.lower().startswith("node"):
                    res += 1
        return res


    def get_cmaboss(self):

        assert self.nb_nodes <= 1024, "Models with more than 1024 nodes are not compatible with this version of MaBoSS"

        if self.nb_nodes <= 64:
            return __import__("cmaboss")
        elif self.nb_nodes <= 128:
            return __import__("cmaboss_128n")
        elif self.nb_nodes <= 256:
            return __import__("cmaboss_256n")
        elif self.nb_nodes <= 512:
            return __import__("cmaboss_512n")
        else:
            return __import__("cmaboss_1024n")

    def get_logical_rules(self):

        out = self.cmaboss_sim.get_logical_rules()    
        rules = {}
        for line in out.split("\n"):
            if ":" in line:
                node, rule = line.split(" : ", 1)
                rules.update({node.strip(): rule.strip()})

        return rules
        
    def update_parameters(self, **kwargs):
        self.cmaboss_sim.update_parameters(**kwargs)

    def run(self, workdir=None, only_final_state=False):
        return CMaBoSSResult2(self, workdir, only_final_state)

    def check(self):
        return self.cmaboss.MaBoSSSim(net=self.cmaboss_net, cfg=self.cmaboss_cfg)

    def str_bnd(self):
        return self.cmaboss_sim.str_bnd()
        
    def str_cfg(self):
        return self.cmaboss_sim.str_cfg()

    def copy(self):
        return CMaBoSSSimulation(cmaboss=self.cmaboss, cmaboss_sim=self.cmaboss_sim.copy())
        
    def print_bnd(self, out=stdout):
        """Produce the content of the bnd file associated to the simulation."""
        print(self.str_bnd(), file=out)

    def print_cfg(self, out=stdout):
        """Produce the content of the cfg file associated to the simulation."""
        print(self.str_cfg(), file=out)

__all__ = ["CMaBoSSSimulation"]
