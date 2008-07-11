import os,sys,math,random
from phycas import *
from MCMCManager import LikelihoodCore
from phycas.Phycas.PhycasCommand import *
from phycas.ReadNexus import NexusReader

class LikeImpl(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    To be written.
    
    """
    def __init__(self, phycas, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the LikeImpl object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas                = phycas
        self.opts                  = opts
        self.starting_tree         = None
        self.taxon_labels          = None
        self.data_matrix           = None
        self.ntax                  = None
        self.nchar                 = None
        self.reader                = NexusReader()
        
    def readDataFromFile(self):
        self.reader.readFile(self.opts.data_file_name)
        self.taxon_labels = self.reader.getTaxLabels()
        self.data_matrix = self.reader.getLastDiscreteMatrix(True)
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current tree and current
        model.
        
        """
        if self.opts.tree_topology.__class__.__name__ == 'str':
            self.starting_tree = Newick(self.opts.tree_topology)
        else:
            self.phycas.phycassert(self.opts.tree_topology.__class__.__name__ == 'Newick', 'expecting tree_topology to be either a string or a Newick object')
            self.starting_tree = self.opts.tree_topology
        self.phycas.phycassert(self.opts.data_file_name is not None, "specify data_file_name before calling like()")
        self.readDataFromFile()
        core = LikelihoodCore(self)
        core.setupCore()
        core.prepareForLikelihood()
        return core.calcLnLikelihood()
        
        