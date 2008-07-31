import os,sys,math
from phycas import *
import phycas.Phylogeny as Phylogeny
import phycas.ProbDist as ProbDist
import phycas.Likelihood as Likelihood

def cloneDistribution(d):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This function clones (deep copies) a ProbabilityDistribution object.
    Cloning is needed because while separate chains in an MCMC analysis 
    need to begin with the same random number seed, it is not good for
    them to share a pseudorandom number generator if (for example) 
    separate chains are run in different threads or on different
    processors.
    
    """
    if d == None:
        return None
    else:
        return d.clone()

class LikelihoodCore:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The LikelihoodCore class comprises only those features of a 
    MarkovChain that are needed to compute likelihoods. Thus, one could 
    construct a stand-alone LikelihoodCore object if all you wanted to do
    was calculate the likelihood of a tree under a particular model. The
    LikelihoodCore class serves as the base class for MarkovChain.
    
    """
    def __init__(self, parent):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The constructor for the LikelihoodCore class stores the supplied 
        parent object in the data member self.parent. It also creates a Tree
        object and stores it in self.tree and a pseudorandom number generator
        (Lot object), storing it in self.r. Finally, it clones the starting
        edge length distribution, storing it in self.starting_edgelen_dist.
        
        """
        self.parent                 = parent
        self.model                  = None
        self.likelihood             = None
        self._tree                  = None
        self.r                      = self.parent._getLot()
    def getTree(self):
        if self._tree is None:
            self.setupCore()
        return self._tree
    def setTree(self, t):
        self._tree = t
    def delTree(self, t):
        if self._tree:
            del self._tree
            self._tree = None
    tree = property(getTree, setTree, delTree)

    def setupCore(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The setupCore function does the following based on information stored
        in self.parent: 1) sets the random number seed of the local pseudo-
        random number generator (self.r); 2) creates a substitution model,
        storing it in self.model; 3) creates a TreeLikelihood object, storing
        it in self.likelihood; and 4) builds the starting tree, storing it in
        self.tree.
        
        """
        # Set seed if user has supplied one
        self.r = self.parent._getLot()
        #self.starting_edgelen_dist.setLot(self.r)

        # Create a substitution model
        if self.parent.opts.model.type in ['gtr','hky']:
            if self.parent.opts.model.type == 'gtr':
                self.model = Likelihood.GTRModel()
                self.model.setRelRates(self.parent.opts.model.relrates)
                if self.parent.opts.model.fix_relrates:
                    self.model.fixRelRates()
            else:
                self.model = Likelihood.HKYModel()
                self.model.setKappa(self.parent.opts.model.kappa)
                if self.parent.opts.model.fix_kappa:
                    self.model.fixKappa()
            self.parent.phycassert(self.parent.opts.model.base_freqs, 'base_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.parent.phycassert(len(self.parent.opts.model.base_freqs) == 4, 'base_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(self.parent.opts.model.base_freqs))
            self.parent.phycassert(self.parent.opts.model.base_freqs[0] >= 0.0, 'base_freqs[0] cannot be negative (%f was specified)' % self.parent.opts.model.base_freqs[0])
            self.parent.phycassert(self.parent.opts.model.base_freqs[1] >= 0.0, 'base_freqs[1] cannot be negative (%f was specified)' % self.parent.opts.model.base_freqs[1])
            self.parent.phycassert(self.parent.opts.model.base_freqs[2] >= 0.0, 'base_freqs[2] cannot be negative (%f was specified)' % self.parent.opts.model.base_freqs[2])
            self.parent.phycassert(self.parent.opts.model.base_freqs[3] >= 0.0, 'base_freqs[3] cannot be negative (%f was specified)' % self.parent.opts.model.base_freqs[3])
            self.model.setNucleotideFreqs(self.parent.opts.model.base_freqs[0], self.parent.opts.model.base_freqs[1], self.parent.opts.model.base_freqs[2], self.parent.opts.model.base_freqs[3])  #POL should be named setStateFreqs?
            if self.parent.opts.model.fix_freqs:
                self.model.fixStateFreqs()
        else:
            self.model = Likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
        # Note must defer setting up pattern specific rates model until we know number of patterns
        if self.parent.opts.model.use_flex_model:
            self.model.setNGammaRates(self.parent.opts.model.num_rates)
            self.model.setFlexModel()
            self.model.setFlexRateUpperBound(self.parent.opts.model.flex_L)
        elif self.parent.opts.model.num_rates > 1:
            self.model.setNGammaRates(self.parent.opts.model.num_rates)
            self.model.setPriorOnShapeInverse(self.parent.opts.model.use_inverse_shape)    #POL should be named useInverseShape rather than setPriorOnShapeInverse
            self.model.setShape(self.parent.opts.model.gamma_shape)
            if self.parent.opts.model.fix_shape:
                self.model.fixShape()
        else:
            self.model.setNGammaRates(1)
            
        if self.parent.opts.model.pinvar_model:
            assert not self.parent.opts.model.use_flex_model, 'Cannot currently use flex model with pinvar'
            self.model.setPinvarModel()
            self.model.setPinvar(self.parent.opts.model.pinvar)
            if self.parent.opts.model.fix_pinvar:
                self.model.fixPinvar()
        else:
            self.model.setNotPinvarModel()

        if self.parent.opts.model.fix_edgelens:
            self.model.fixEdgeLengths()
            
        # Create the likelihood object
        self.likelihood = Likelihood.TreeLikelihood(self.model)
        self.likelihood.setUFNumEdges(self.parent.opts.uf_num_edges)
        self.likelihood.useUnimap(self.parent.opts.use_unimap)
        if self.parent.data_matrix:
            self.likelihood.copyDataFromDiscreteMatrix(self.parent.data_matrix)
            #POL changed self.npatterns to self.parent.npatterns below
            self.parent.npatterns = self.likelihood.getNPatterns()

        # Build the starting tree
        self.tree = self.parent.getStartingTree()

    def prepareForSimulation(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds data structures (conditional likelihood arrays and transition
        probability matrices) to nodes in tree to facilitate simulating data.
        Unlike prepareForLikelihood function, does not add data to the tips.
        
        """
        self.likelihood.prepareForSimulation(self.tree)
        
    def simulate(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a simulation, returning a new SimData object containing the
        simulated data.
        
        """
        sim_data = Likelihood.SimData()
        self.parent.phycassert(self.parent.opts.nchar > 0, 'nchar must be greater than zero in order to perform simulations')
        self.likelihood.simulateFirst(sim_data, self.tree, self.r, self.parent.opts.nchar)
        return sim_data
        
    def prepareForLikelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds data structures (conditional likelihood arrays and transition
        probability matrices) to nodes in tree to facilitate likelihood
        calculations. Also adds the data to the tip nodes.
        
        """
        self.likelihood.prepareForLikelihood(self.tree)
        
    def calcLnLikelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes and returns the log-likelihood.
        
        """
        return self.likelihood.calcLnL(self.tree)

    def debugCheckForUncachedCLAs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Scans tree for cached CLAs. If any are found, returns True. If no CLAs
        are currently cached, returns False. False is the expected result,
        because there should only be cached CLAs in the tree if we are in the
        middle of performing a Metropolis-Hastings move. If a move has been
        proposed, but not yet accepted or rejected, cached CLAs provide a way
        to return to the pre-move state after rejection without having to
        recalculate the likelihood.
        
        """
        return self.likelihood.debugCheckForUncachedCLAs(self.tree)
        
class MarkovChain(LikelihoodCore):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The MarkovChain class encapsulates the notion of a Markov chain used
    in Bayesian analyses. The MarkovChain class has the ability to
    orchestrate a Markov chain Monte Carlo (MCMC) analysis. In 
    Metropolis-coupled MCMC, the cold chain and each of the heated chains
    are MarkovChain objects. The primary addition MarkovChain adds to the
    base class LikelihoodCore are prior distributions and the 
    self.heating_power data member.
    
    """
    def __init__(self, parent, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The MarkovChain constructor makes a copy of the supplied parent
        object, clones all of the prior ProbabilityDistribution objects in the
        supplied parent object, and sets self.heating_power to the supplied
        power.
        
        """
        LikelihoodCore.__init__(self, parent)
        
        self.parent                  = parent
        self.heating_power           = power
        self.relrate_prior           = cloneDistribution(self.parent.opts.model.relrate_prior)
        self.base_freq_param_prior   = cloneDistribution(self.parent.opts.model.base_freq_param_prior)
        self.gamma_shape_prior       = cloneDistribution(self.parent.opts.model.gamma_shape_prior)
        self.external_edgelen_dist   = cloneDistribution(self.parent.opts.model.external_edgelen_dist)
        self.internal_edgelen_dist   = cloneDistribution(self.parent.opts.model.internal_edgelen_dist)
        self.kappa_prior             = cloneDistribution(self.parent.opts.model.kappa_prior)
        self.pinvar_prior            = cloneDistribution(self.parent.opts.model.pinvar_prior)
        self.flex_prob_param_prior   = cloneDistribution(self.parent.opts.model.flex_prob_param_prior)
        self.edgelen_hyperprior      = None
        if self.parent.opts.model.edgelen_hyperprior is not None:
            self.edgelen_hyperprior  = cloneDistribution(self.parent.opts.model.edgelen_hyperprior)
        self.chain_manager           = None

        self.setupChain()

    def resetNEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the resetNEvals function of self.likelihood. This resets the 
        number of likelihood evaluations performed to 0.

        """
        return self.likelihood.resetNEvals()
    
    def getNEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the getNEvals function of self.likelihood. This returns the 
        number of likelihood evaluations performed since the last call of
        resetNEvals.

        """
        return self.likelihood.getNEvals()
        
    def paramFileHeader(self, paramf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the parameter file. The
        parameter paramf should be an (already open for writing) file object.
        The parameter file corresponds to the *.p file produced by MrBayes. It
        begins with a line containing the initial random number seed, followed
        by a line of column titles for the (tab-delimited) sampled parameter
        information on the subsequent lines.

        """
        paramf.write('[ID: %d]\n' % self.r.getInitSeed())
        if self.parent.opts.doing_path_sampling:
            paramf.write('Gen\tbeta\tLnL\tTL')
        else:
            paramf.write('Gen\tLnL\tTL')
        paramf.write(self.model.paramHeader())
        if self.parent.opts.model.edgelen_hyperprior is not None:
            if self.parent.opts.model.internal_edgelen_dist is self.parent.opts.model.external_edgelen_dist:
                paramf.write('\thyper(all)')
            else:
                paramf.write('\thyper(external)')
                paramf.write('\thyper(internal)')
        if self.parent.opts.model.use_flex_model:
            paramf.write('\trates_probs')

    def treeFileHeader(self, treef):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the tree file. The
        parameter treef should be an (already open for writing) file object.
        The tree file corresponds to the *.t file produced by MrBayes, and is
        in NEXUS format. This function writes the opening "#NEXUS" line, a 
        comment containing the initial random number seed, and the beginning
        of a NEXUS TREES block, including the TRANSLATE command.

        """
        treef.write('#NEXUS\n')
        treef.write('[ID: %d]\n' % self.r.getInitSeed())
        treef.write('begin trees;\n')
        treef.write('   translate\n')
        for i in range(self.parent.ntax):
            if self.parent.taxon_labels[i].find(' ') < 0:
                # no spaces found in name
                treef.write('       %d %s%s\n' % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))
            else:
                # at least one space in taxon name, so enclose name in quotes
                treef.write("       %d '%s'%s\n" % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))

    def setPower(self, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the heating_power data member and calls the setPower method for
        every updater so that all updaters are immediately informed that the
        power for this chain has changed.
        
        """
        self.heating_power = power
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(power)

    def setupChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The setupChain method prepares a MarkovChain object before an MCMC
        analysis is started. This method is called at the end of the 
        MarkovChain constructor. The duties performed include:
        1) calling the base class setupCore method; 
        2) preparing the starting tree by adding data structures needed for
        computing likelihood (transition matrices and conditional likelihood
        arrays); 
        3) setting up prior distributions for model parameters; 
        4) creating an MCMCManager object and adding all relevant updaters to
        it; 
        and 5) makes sure each updater knows the type of heating and the 
        power.
        
        """
        LikelihoodCore.setupCore(self)
        LikelihoodCore.prepareForLikelihood(self)
        
        # Make sure that each prior is using the same pseudorandom number generator object
        self.relrate_prior.setLot(self.r)
        self.base_freq_param_prior.setLot(self.r)
        self.flex_prob_param_prior.setLot(self.r)
        self.gamma_shape_prior.setLot(self.r)
        self.pinvar_prior.setLot(self.r)
        if self.edgelen_hyperprior is not None:
            self.edgelen_hyperprior.setLot(self.r)
        self.external_edgelen_dist.setLot(self.r)
        if self.parent.opts.model.internal_edgelen_dist:
            self.internal_edgelen_dist.setLot(self.r)
        
        # Define priors for the model parameters
        if self.parent.opts.model.type == 'gtr':
            self.model.setRelRatePrior(self.relrate_prior)
            self.model.setStateFreqParamPrior(self.base_freq_param_prior)   #POL should be named state_freq_param_prior
        elif self.parent.opts.model.type == 'hky':
            self.model.setKappaPrior(self.kappa_prior)
            self.model.setStateFreqParamPrior(self.base_freq_param_prior)   #POL should be named state_freq_param_prior

        # If rate heterogeneity is to be assumed, add priors for these model parameters here
        if self.parent.opts.model.use_flex_model:
            self.model.setNumFlexSpacers(self.parent.opts.model.flex_num_spacers)
            self.model.setFLEXProbParamPrior(self.flex_prob_param_prior)
        elif self.parent.opts.model.num_rates > 1:
            self.model.setDiscreteGammaShapePrior(self.gamma_shape_prior)
        if self.parent.opts.model.pinvar_model:
            self.model.setPinvarPrior(self.pinvar_prior)
        
        # Define edge length prior distributions
        separate_edge_len_dists = self.parent.opts.model.internal_edgelen_dist is not self.parent.opts.model.external_edgelen_dist
        self.model.separateInternalExternalEdgeLenPriors(separate_edge_len_dists)
        self.model.setExternalEdgeLenPrior(self.external_edgelen_dist)
        self.model.setInternalEdgeLenPrior(self.internal_edgelen_dist)

        if self.parent.opts.model.edgelen_hyperprior is not None:
            #self.edgelen_hyperprior.setMeanAndVariance(1.0, 10.0)
            self.model.setEdgeLenHyperPrior(self.edgelen_hyperprior)
            #todo self.model.starting_edgelen_hyperparam
            if self.parent.opts.model.fix_edgelen_hyperparam:
                self.model.fixEdgeLenHyperprior()   #POL should be named fixEdgeLenHyperparam
        else:
            self.model.setEdgeLenHyperPrior(None)

        self.likelihood.replaceModel(self.model)            

        if self.parent.opts.data_source is None:
            self.likelihood.setNoData() # user apparently wants to run MCMC with no data
        
        # Create an MCMCChainManager object and add all necessary updaters
        self.chain_manager = Likelihood.MCMCChainManager()
        self.chain_manager.addMCMCUpdaters(
            self.model,                     # substitution model
            self.tree,                      # tree
            self.likelihood,                # likelihood calculation machinery
            self.r,                         # pseudorandom number generator
            #POLPY_NEWWAY False,            # separate_edgelen_params (deprecated: always False)
            self.parent.opts.slice_max_units,    # maximum number of slice units allowed
            self.parent.opts.slice_weight)       # weight for each parameter added

        # Create a TreeScalerMove object to handle scaling the entire tree to allow faster
        # convergence in edge lengths. This move is unusual in using slice sampling rather
        # than Metropolis-Hastings updates: most "moves" in parent are Metropolis-Hastings.
        if self.parent.opts.tree_scaler_weight > 0:
            self.tree_scaler_move = Likelihood.TreeScalerMove()
            self.tree_scaler_move.setName("Tree scaler move")
            self.tree_scaler_move.setWeight(self.parent.opts.tree_scaler_weight)
            self.tree_scaler_move.setTree(self.tree)
            self.tree_scaler_move.setModel(self.model)
            self.tree_scaler_move.setTreeLikelihood(self.likelihood)
            self.tree_scaler_move.setLot(self.r)
            if self.model.edgeLengthsFixed():
                self.tree_scaler_move.fixParameter()
            self.chain_manager.addMove(self.tree_scaler_move)

        if self.parent.opts.use_unimap:
            # Create a MappingMove (to refresh the mapping for all sites)
            self.unimapping_move = Likelihood.MappingMove()
            self.unimapping_move.setName("Univent mapping move")
            self.unimapping_move.setWeight(self.parent.opts.mapping_move_weight)
            self.unimapping_move.setTree(self.tree)
            self.unimapping_move.setModel(self.model)
            self.unimapping_move.setTreeLikelihood(self.likelihood)
            self.unimapping_move.setLot(self.r)
            self.chain_manager.addMove(self.unimapping_move)

            # Create a UnimapNNIMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_nni_move = Likelihood.UnimapNNIMove()
            self.unimap_nni_move.setName("Unimap NNI move")
            self.unimap_nni_move.setWeight(self.parent.opts.unimap_nni_move_weight)
            self.unimap_nni_move.setTree(self.tree)
            self.unimap_nni_move.setModel(self.model)
            self.unimap_nni_move.setTreeLikelihood(self.likelihood)
            self.unimap_nni_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_nni_move)

        elif self.parent.opts.fix_topology:
            # Create an EdgeMove object to handle Metropolis-Hastings
            # updates to the edge lengths only (does not change the topology)
            self.edge_move = Likelihood.EdgeMove()
            self.edge_move.setName("Edge length move")
            self.edge_move.setWeight(self.parent.opts.edge_move_weight)
            self.edge_move.setTree(self.tree)
            self.edge_move.setModel(self.model)
            self.edge_move.setTreeLikelihood(self.likelihood)
            self.edge_move.setLot(self.r)
            self.edge_move.setLambda(self.parent.opts.edge_move_lambda)
            if self.model.edgeLengthsFixed():
                self.edge_move.fixParameter()
            self.chain_manager.addMove(self.edge_move)
        else:
            # Create a LargetSimonMove object to handle Metropolis-Hastings
            # updates to the tree topology and edge lengths
            self.larget_simon_move = Likelihood.LargetSimonMove()
            self.larget_simon_move.setName("Larget-Simon move")
            self.larget_simon_move.setWeight(self.parent.opts.ls_move_weight)
            self.larget_simon_move.setTree(self.tree)
            self.larget_simon_move.setModel(self.model)
            self.larget_simon_move.setTreeLikelihood(self.likelihood)
            self.larget_simon_move.setLot(self.r)
            self.larget_simon_move.setLambda(self.parent.opts.ls_move_lambda)
            if self.model.edgeLengthsFixed():
                self.larget_simon_move.fixParameter()
            self.chain_manager.addMove(self.larget_simon_move)

        # If requested, create an NCatMove object to allow the number of rate categories to change
        if self.parent.opts.model.use_flex_model:
            # Create an NCatMove object
            self.ncat_move = Likelihood.NCatMove()
            
            # Set up features specific to NCatMove
            self.ncat_move.setCatProbPrior(self.flex_prob_param_prior)
            self.ncat_move.setL(self.parent.opts.model.flex_L)
            self.ncat_move.setS(self.parent.opts.model.flex_num_spacers)
            self.ncat_move.setLambda(self.parent.opts.model.flex_lambda)
            self.ncat_move.setPhi(self.parent.opts.model.flex_phi)

            # Continue setting up NCatMove object
            self.ncat_move.setName("NCat move")
            self.ncat_move.setWeight(self.parent.opts.model.flex_ncat_move_weight)
            self.ncat_move.setTree(self.tree)
            self.ncat_move.setModel(self.model)
            self.ncat_move.setTreeLikelihood(self.likelihood)
            self.ncat_move.setLot(self.r)
            
            self.chain_manager.addMove(self.ncat_move)
            
        # If requested, create a BushMove object to allow polytomous trees
        if self.parent.opts.allow_polytomies:
            # Create a BushMove object
            self.bush_move = Likelihood.BushMove()

            # Set up the topology prior
            self.topo_prior_calculator = self.bush_move.getTopoPriorCalculator()
            self.topo_prior_calculator.chooseUnrooted()
            self.topo_prior_calculator.setC(self.parent.opts.topo_prior_C)
            if self.parent.opts.polytomy_prior:
                self.topo_prior_calculator.choosePolytomyPrior()
            else:
                self.topo_prior_calculator.chooseResolutionClassPrior()
                
            # Continue setting up BushMove object
            self.bush_move.setName("Bush move")
            self.bush_move.setWeight(self.parent.opts.bush_move_weight)
            self.bush_move.setTree(self.tree)
            self.bush_move.setModel(self.model)
            self.bush_move.setTreeLikelihood(self.likelihood)
            self.bush_move.setLot(self.r)
            self.bush_move.setEdgeLenDistMean(self.parent.opts.bush_move_edgelen_mean)
            #self.bush_move.viewProposedMove(self.parent.bush_move_debug)
            if self.model.edgeLengthsFixed():
                self.bush_move.fixParameter()
            self.bush_move.finalize()
            
            self.chain_manager.addMove(self.bush_move)

        # REVISIT LATER
        #if self.parent.doing_samc:
        #    self.samc_move = Likelihood.SamcMove(self.starting_edgelen_dist)
        #
        #    # Continue setting up SAMC move object
        #    self.samc_move.setName("SAMC move")
        #    self.samc_move.setWeight(self.parent.samc_move_weight)
        #    self.samc_move.setTree(self.tree)
        #    self.samc_move.setModel(self.model)
        #    self.samc_move.setTreeLikelihood(self.likelihood)
        #    self.samc_move.setLot(self.r)
        #    if self.model.edgeLengthsFixed():
        #        self.samc_move.fixParameter()
        #    #self.samc_move.finalize()
        #    self.chain_manager.addMove(self.samc_move)

        self.chain_manager.finalize()

        # Calculate relative rates if among-site rate heterogeneity is part of the model
        # Currently, the unnormalized rates are stored in the model; the function call
        # below normalizes the rate means and probabilities so that they will be accurately
        # recorded in the first line of the parameter file
        self.likelihood.recalcRelativeRates()

        # Make sure each updater knows the heating power and heating type
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(self.heating_power)
            #if self.parent.opts.is_standard_heating:
            updater.setStandardHeating()
            #else:
            #    updater.setLikelihoodHeating()

class MCMCManager:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This class manages one to many Markov chains participating in a
    Bayesian MCMC analysis. It creates the MarkovChain objects and manages
    swapping and sampling the chains. If likelihood heating is employed,
    it also has the ability to estimate (albeit crudely) the marginal
    likelihood of the model. A single instance of MCMCManager is created
    by parent in the parent contructor.
    
    """
    def __init__(self, parent):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Stores the parent object passed into the constructor and creates an
        empty self.chains list.
        
        """
        self.parent = parent
        self.chains = []

    def paramFileHeader(self, paramf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simply passes the file object paramf on to the paramFileHeader method
        of the MarkovChain object representing the cold chain.
        
        """
        m = self.getColdChain()
        m.paramFileHeader(paramf)

    def treeFileHeader(self, treef):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simply passes the file object treef on to the treeFileHeader method
        of the MarkovChain object representing the cold chain.
                
        """
        m = self.getColdChain()
        m.treeFileHeader(treef)

    def createChains(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a separate MarkovChain for every element in parent.heat_vector
        and adds it to self.chains. This function is also responsible for 
        performing several sanity checks (last chance to abort before creating
        the MCMC machinery).
        
        """
        # Sanity checks
        unimap_and_flex            = (self.parent.opts.use_unimap and self.parent.opts.model.use_flex_model)
        unimap_and_ratehet         = (self.parent.opts.use_unimap and self.parent.opts.model.num_rates > 1)
        unimap_and_polytomies      = (self.parent.opts.use_unimap and self.parent.opts.allow_polytomies)
        unimap_and_multiple_chains = (self.parent.opts.use_unimap and self.parent.opts.nchains > 1)
        #unimap_and_samc            = (self.parent.opts.use_unimap and self.parent.opts.doing_samc)
        #self.parent.phycassert(not unimap_and_samc, 'SAMC cannot (yet) be used in conjunction with use_unimap')
        self.parent.phycassert(not unimap_and_polytomies, 'Allowing polytomies cannot (yet) be used in conjunction with use_unimap')
        self.parent.phycassert(not unimap_and_flex, 'Flex model cannot (yet) be used in conjunction with use_unimap')
        self.parent.phycassert(not unimap_and_ratehet, 'Rate heterogeneity cannot (yet) be used in conjunction with use_unimap')
        self.parent.phycassert(not unimap_and_multiple_chains, 'Multiple chains cannot (yet) be used in conjunction with use_unimap')
        
        # Create the chains
        for heating_power in self.parent.heat_vector:
            markov_chain = MarkovChain(self.parent, heating_power)
            self.chains.append(markov_chain)

    def getNumChains(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of chains being managed by this object
        (i.e. returns the length of the self.chains list).
        
        """
        return len(self.chains)

    def getColdChainIndex(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the index (0..nchains-1) of the chain corresponding to the
        cold chain (i.e. the first chain whose power equals 1.0).
        
        """
        return self.parent.heat_vector.index(1.0)

    def getColdChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the cold chain MarkovChain object. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        i = self.parent.heat_vector.index(1.0)
        return self.chains[i]

    def getColdChainManager(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the chain manager for the cold chain. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        i = self.parent.heat_vector.index(1.0)
        return self.chains[i].chain_manager

    def setChainPower(self, chain_index, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Loops through all updaters in the chain having the specified 
        chain_index, setting the power for each updater to the supplied power.
        
        """
        parent.phycassert(len(chains) > chain_index, 'chain index specified (%d) too large for number of chains (%d)' % (chain_index, len(chains)))
        parent.phycassert(len(self.parent.heat_vector) == len(self.chains), 'length of heat vector (%d) not equal to number of chains (%d)' % (len(self.heat_vector, len(self.nchains))))
        self.parent.heat_vector[chain_index] = power
        self.chains[chain_index].setPower(power)

    def resetNEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls resetNEvals method for every MarkovChain in the self.chains
        list.
        
        """
        for c in self.chains:
            c.resetNEvals()

    def setRandomSeedAllChains(self, rnseed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the master random number seed to rnseed, then passes rnseed to
        the setSeed method of every MarkovChain in the self.chains list.
        
        """
        self.parent.opts.random_seed = rnseed
        for c in self.chains:
            c.r.setSeed(int(rnseed))

    def getTotalEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of likelihood evaluations over all 
        MarkovChain objects in the self.chains list by calling the getNEvals
        method for each chain.
        
        """
        total = 0
        for c in self.chains:
            total += c.getNEvals()
        return total

    def recordSample(self, cycle):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Records the current tree topology and edge lengths by adding a line to
        the tree file, and records tree length and substitution parameters
        by adding a line to the parameter file.
        
        """
        # Gather log-likelihoods, and if path sampling save in path_sample list for later
        lnLikes = []
        #multichain_path_sampling = self.parent.opts.nchains > 1 and not self.parent.opts.is_standard_heating
        for i,c in enumerate(self.chains):
            lnLi = c.chain_manager.getLastLnLike()
            lnLikes.append(lnLi)
            #if multichain_path_sampling:
            #    self.parent.path_sample[i].append(lnLi) # DISCRETE PATH SAMPLING
        
        # Only record samples from the current cold chain
        cold_chain = self.parent.mcmc_manager.getColdChain()
        
        # Add line to parameter file if it exists
        if self.parent.paramf:
            self.parent.paramf.write('%d\t' % (cycle + 1))
            if self.parent.opts.doing_path_sampling:
                self.parent.paramf.write('%.5f\t' % (cold_chain.heating_power))
            for lnl in lnLikes:
                self.parent.paramf.write('%.3f\t' % lnl)
            self.parent.paramf.write('%.3f' % cold_chain.tree.edgeLenSum())
            
            self.parent.paramf.write(cold_chain.model.paramReport())
            if self.parent.opts.model.edgelen_hyperprior is not None:
                for p in cold_chain.chain_manager.getEdgeLenHyperparams():
                    self.parent.paramf.write('\t%.5f' % p.getCurrValue())
            if self.parent.opts.model.use_flex_model:
                rates_vector = cold_chain.likelihood.getRateMeans()
                for rr in rates_vector:
                    self.parent.paramf.write('\t%.5f' % rr)
                probs_vector = cold_chain.likelihood.getRateProbs()
                for rp in probs_vector:
                    self.parent.paramf.write('\t%.5f' % rp)
            self.parent.paramf.write('\n')
        
        # Add line to tree file if it exists
        if self.parent.treef:
            self.parent.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, cold_chain.tree.makeNewick()))

