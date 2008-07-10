import os,sys,math,random
from phycas import *

def cloneDistribution(d):
    if d == None:
        return None
    else:
        return d.clone()

class UnimapManager(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Orchestrates a uniformized mapping MCMC analysis.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes UnimapManager object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas = phycas
        self.model                              = None
        self.likelihood                         = None
        self.tree                               = Phylogeny.Tree()
        self.r                                  = ProbDist.Lot()
        self.starting_edgelen_dist              = cloneDistribution(self.phycas.starting_edgelen_dist)
        
    def setup(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Copied largely from MCMCManager::setupCore.
        
        """
        # Set seed if user has supplied one
        if self.phycas.random_seed != 0:
            self.r.setSeed(int(self.phycas.random_seed))

        self.starting_edgelen_dist.setLot(self.r)

        # Create a substitution model
        if self.phycas.default_model in ['gtr','hky']:
            if self.phycas.default_model == 'gtr':
                self.model = Likelihood.GTRModel()
                self.model.setRelRates(self.phycas.relrates)
                if self.phycas.fix_relrates:
                    self.model.fixRelRates()
            else:
                self.model = Likelihood.HKYModel()
                self.model.setKappa(self.phycas.kappa)
                if self.phycas.fix_kappa:
                    self.model.fixKappa()
            self.phycas.phycassert(self.phycas.base_freqs, 'base_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.phycas.phycassert(len(self.phycas.base_freqs) == 4, 'base_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(self.phycas.base_freqs))
            self.phycas.phycassert(self.phycas.base_freqs[0] >= 0.0, 'base_freqs[0] cannot be negative (%f was specified)' % self.phycas.base_freqs[0])
            self.phycas.phycassert(self.phycas.base_freqs[1] >= 0.0, 'base_freqs[1] cannot be negative (%f was specified)' % self.phycas.base_freqs[1])
            self.phycas.phycassert(self.phycas.base_freqs[2] >= 0.0, 'base_freqs[2] cannot be negative (%f was specified)' % self.phycas.base_freqs[2])
            self.phycas.phycassert(self.phycas.base_freqs[3] >= 0.0, 'base_freqs[3] cannot be negative (%f was specified)' % self.phycas.base_freqs[3])
            self.model.setNucleotideFreqs(self.phycas.base_freqs[0], self.phycas.base_freqs[1], self.phycas.base_freqs[2], self.phycas.base_freqs[3])  #POL should be named setStateFreqs?
            if self.phycas.fix_freqs:
                self.model.fixStateFreqs()
        else:
            self.model = Likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
        self.phycas.phycassert(not self.phycas.use_flex_model, 'Flex model not allowed in uniformized mapping MCMC')
        self.phycas.phycassert(self.phycas.num_rates == 1, 'Among-site rate heterogeneity submodels are not allowed in uniformized mapping MCMC at the present time')
        self.model.setNGammaRates(1)

        self.phycas.phycassert(not self.phycas.pinvar_model, 'Proportion of invariable sites submodel is not allowed in uniformized mapping MCMC at the present time')
        self.model.setNotPinvarModel()

        if self.phycas.fix_edgelens:
            self.model.fixEdgeLengths()
            
        # Create the TreeLikelihood object
        self.likelihood = Likelihood.TreeLikelihood(self.model)
        self.likelihood.useUnimap(True) # ! new function
        self.likelihood.setUFNumEdges(self.phycas.uf_num_edges)
        self.phycas.phycassert(self.phycas.data_source == 'file', "At this time, only data_source 'file' allowed when using uniformized mapping MCMC")
        self.phycas.readDataFromFile()
        # ! modified copyDataFromDiscreteMatrix to avoid compressing data patterns in case of unimap MCMC
        self.likelihood.copyDataFromDiscreteMatrix(self.phycas.data_matrix)
        self.phycas.npatterns = self.likelihood.getNPatterns()

        # Build the starting tree
        if self.phycas.starting_tree_source == 'random':
            # Build a random tree
            Phylogeny.TreeManip(self.tree).randomTree(
                self.phycas.ntax,           # number of tips
                self.r,                     # pseudorandom number generator
                self.starting_edgelen_dist, # distribution from which to draw starting edge lengths
                False)                      # Yule tree if True, edge lengths independent if False
            self.phycas.warn_tip_numbers = False
            self.phycas.starting_tree = self.tree.makeNewick()
        elif self.phycas.starting_tree_source == 'usertree':
            # Build user-specified tree
            self.phycas.starting_tree = self.phycas.tree_topology
            self.phycas.starting_tree.buildTree(self.tree)
            if not self.tree.tipNumbersSetUsingNames():
                self.phycas.warn_tip_numbers = True
        else:
            self.phycas.phycassert(False, "starting_tree_source must be either 'random' or 'usertree'")
                
        # prepare tree for unimap MCMC
        self.likelihood.prepareForLikelihood(self.tree) 

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initiates a uniformized mapping MCMC analysis.
        
        """
        self.setup()
        self.likelihood.nielsenMapping(self.tree, self.r) 
        self.phycas.output('Unimap MCMC analysis finished.')
       

        
        
    
