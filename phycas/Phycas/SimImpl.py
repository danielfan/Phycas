import os,sys,math,random
from phycas import *
from MCMCManager import LikelihoodCore
from phycas.Utilities.PhycasCommand import *

class SimImpl(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Simulates DNA sequence data.
    
    """
    def __init__(self, phycas, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the SimulateImpl object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas               = phycas
        self.opts                 = opts
        self.starting_tree_source = None
        self.starting_tree        = None
        self.ntax                 = None
        self.sim_model_tree       = None
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simulates a DNA dataset and stores it in NEXUS format in the supplied
        filename.
        
        """
        self.starting_tree_source = self.opts.tree_source
        self.starting_tree        = self.opts.tree_topology
        self.ntax                 = len(self.opts.taxon_labels)
        self.phycas.phycassert(self.ntax > 3, 'Must specify labels for at least four taxa')
        core = LikelihoodCore(self)
        core.setupCore()
        if not core.tree.hasEdgeLens():
            tm = TreeManip(core.tree)
            tm.setRandomEdgeLengths(core.starting_edgelen_dist)
        self.sim_model_tree = core.tree
        core.prepareForSimulation()
        sim_data = core.simulate()
        sim_data.saveToNexusFile(self.opts.file_name, self.opts.taxon_labels, 'dna', ('a','c','g','t'))
        
