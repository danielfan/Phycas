#
# NOTE: This file no longer being used as of Lawrence meeting May 2009
#

import os,sys,math,random
from phycas import *
from MCMCManager import LikelihoodCore
from phycas.Utilities.PhycasCommand import *

def cloneDistribution(d):
    if d == None:
        return None
    else:
        return d.clone()

class UnimapImpl(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Orchestrates a uniformized mapping MCMC analysis.
    
    """
    def __init__(self, phycas, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes UnimapManager object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas                 = phycas
        self.opts                   = opts
        self.model                  = None
        self.likelihood             = None
        self.tree                   = Phylogeny.Tree()
        self.r                      = ProbDist.Lot()
        self.starting_edgelen_dist  = cloneDistribution(self.phycas.starting_edgelen_dist)
        self.ntax                   = 0         # Will hold the actual number of taxa after data file read
        self.reader                 = NexusReader()
        
    def readDataFromFile(self):
        self.reader.readFile(self.opts.data_file_name)
        self.taxon_labels = self.reader.getTaxLabels()
        self.data_matrix = self.reader.getLastDiscreteMatrix(True)
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only

    def setup(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is for parts of the setup that should occur right before
        run() is called. Setup is deferred until this point to give the
        user a chance to change the default settings before the call to
        run(). This function does these things: 
        1) reads the data and ensures that the taxon_labels list is filled
        with the correct number of taxon labels; 
        2) creates a starting tree description; 
        3) creates an appropriate heat_vector (for either Metropolis-coupling
        or discrete path sampling); 
        4) calls MCMCManager's createChains function to handle setup for each
        individual chain; and 
        5) opens the parameter and tree files.
        
        """
        # Section below copied largely from MCMCImpl.setup
        
        # REVISIT LATER
        # Hack until PhycasCommand implements properties
        if self.opts.internal_edgelen_prior is None:
            self.opts.internal_edgelen_prior = self.opts.edgelen_prior
        if self.opts.external_edgelen_prior is None:
            self.opts.external_edgelen_prior = self.opts.edgelen_prior
        
        # Read the data
        if self.opts.data_source == 'file':
            self.readDataFromFile()
        elif (len(self.taxon_labels) != self.opts.ntax):
            # Create and store a list of default taxon labels
            self.ntax = self.opts.ntax
            for i in range(self.ntax):
                s = 'taxon_%d' % (i + 1)
                self.taxon_labels.append(s)
        self.phycas.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")

        # Create a tree description to be used for building starting trees
        if self.opts.starting_tree_source == 'file':
            self.phycas.phycas.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")
            tree_defs = self.reader.getTrees()
            self.phycas.phycassert(len(tree_defs) > 0, 'a trees block defining at least one tree must be stored in the nexus data file')
            # Grab first tree description in the data file
            # TODO allow some other tree than the first
            self.starting_tree = tree_defs[0]
        elif self.opts.starting_tree_source == 'usertree':
            self.starting_tree = self.tree_topology
        elif self.opts.starting_tree_source == 'random':
            self.phycas.phycassert(self.ntax > 0, 'expecting ntax to be greater than 0')
            self.starting_tree = None
        else:
            self.phycas.phycassert(False, "starting_tree_source should equal 'random', 'file', or 'usertree', but instead it was this: %s" % self.starting_tree_source)
        
        # Determine heating levels if multiple chains
        if self.heat_vector == None:
        	self.phycassert(self.opts.nchains)
            if self.opts.nchains == 1:
                self.heat_vector = [1.0]
            else:
                # Determine vector of powers for each chain
                self.heat_vector = []
                for i in range(self.opts.nchains):
                    # Standard heating 
                    # 0 1.000 = 1/1.0 cold chain explores posterior
                    # 1 0.833 = 1/1.2
                    # 2 0.714 = 1/1.4
                    # 3 0.625 = 1/1.6
                    temp = 1.0/(1.0 + float(i)*self.heating_lambda)
                    self.heat_vector.append(temp)
        else:
        	self.phycassert(len(self.heat_vector) == 1)
            # User supplied his/her own heat_vector; perform sanity checks
            self.opts.nchains = len(self.heat_vector)
            self.phycas.phycassert(self.heat_vector.index(1.0) < self.opts.nchains, 'user-supplied heat_vector does not allow for a cold chain (one power must be 1.0)')

        self.mcmc_manager.createChains()
        self.openParameterAndTreeFiles()
        
        # Section below copied largely from MCMCManager.setupCore
        
        # Set seed if user has supplied one
        if self.opts.random_seed != 0:
            self.r.setSeed(int(self.opts.random_seed))

        self.starting_edgelen_dist.setLot(self.r)

        # Create a substitution model
        if self.opts.default_model in ['gtr','hky']:
            if self.opts.default_model == 'gtr':
                self.model = Likelihood.GTRModel()
                self.model.setRelRates(self.opts.relrates)
                if self.opts.fix_relrates:
                    self.model.fixRelRates()
            else:
                self.model = Likelihood.HKYModel()
                self.model.setKappa(self.opts.kappa)
                if self.opts.fix_kappa:
                    self.model.fixKappa()
            self.phycas.phycassert(self.opts.state_freqs, 'state_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.phycas.phycassert(len(self.opts.state_freqs) == 4, 'state_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(self.opts.state_freqs))
            self.phycas.phycassert(self.opts.state_freqs[0] >= 0.0, 'state_freqs[0] cannot be negative (%f was specified)' % self.opts.state_freqs[0])
            self.phycas.phycassert(self.opts.state_freqs[1] >= 0.0, 'state_freqs[1] cannot be negative (%f was specified)' % self.opts.state_freqs[1])
            self.phycas.phycassert(self.opts.state_freqs[2] >= 0.0, 'state_freqs[2] cannot be negative (%f was specified)' % self.opts.state_freqs[2])
            self.phycas.phycassert(self.opts.state_freqs[3] >= 0.0, 'state_freqs[3] cannot be negative (%f was specified)' % self.opts.state_freqs[3])
            self.model.setNucleotideFreqs(self.opts.state_freqs[0], self.opts.state_freqs[1], self.opts.state_freqs[2], self.opts.state_freqs[3])  #POL should be named setStateFreqs?
            if self.opts.fix_freqs:
                self.model.fixStateFreqs()
        else:
            self.model = Likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
		
        self.phycas.phycassert(self.opts.num_rates == 1, 'Among-site rate heterogeneity submodels are not allowed in uniformized mapping MCMC at the present time')
        self.model.setNGammaRates(1)

        self.phycas.phycassert(not self.opts.pinvar_model, 'Proportion of invariable sites submodel is not allowed in uniformized mapping MCMC at the present time')
        self.model.setNotPinvarModel()

        if self.opts.fix_edgelens:
            self.model.fixEdgeLengths()
            
        # Create the TreeLikelihood object
        self.likelihood = Likelihood.TreeLikelihood(self.model)
        self.likelihood.useUnimap(True) # ! new function
        self.likelihood.setUFNumEdges(self.opts.uf_num_edges)
        self.phycas.phycassert(self.phycas.data_source == 'file', "At this time, only data_source 'file' allowed when using uniformized mapping MCMC")
        self.phycas.readDataFromFile()
        # ! modified copyDataFromDiscreteMatrix to avoid compressing data patterns in case of unimap MCMC
        self.likelihood.copyDataFromDiscreteMatrix(self.opts.data_matrix)
        self.opts.npatterns = self.likelihood.getNPatterns()

        # Build the starting tree
        if self.opts.starting_tree_source == 'random':
            # Build a random tree
            Phylogeny.TreeManip(self.tree).randomTree(
                self.ntax,           # number of tips
                self.r,                     # pseudorandom number generator
                self.starting_edgelen_dist, # distribution from which to draw starting edge lengths
                False)                      # Yule tree if True, edge lengths independent if False
            self.warn_tip_numbers = False
            self.starting_tree = self.tree.makeNewick()
        elif self.phycas.starting_tree_source == 'usertree':
            # Build user-specified tree
            self.starting_tree = self.opts.tree_topology
            self.starting_tree.buildTree(self.tree)
            if not self.tree.tipNumbersSetUsingNames():
                self.warn_tip_numbers = True
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
        self.likelihood.fullRemapping(self.tree, self.r) 
        self.phycas.output('Unimap MCMC analysis finished.')
       
