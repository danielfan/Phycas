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

class LikelihoodCore(object):
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
        self.parent                 = parent    # e.g. LikeImpl, MCMCImpl, etc.
        self.model                  = None
        self.likelihood             = None
        self._tree                  = None
        self.r                      = self.parent._getLot()

    def getTree(self):
        if self._tree is None:
            self.setupCore()
        return self._tree
        
    def setTree(self, t):
        if t and not isinstance(t, Phylogeny.Tree):
            tr = Phylogeny.Tree()
            self._tree = t.buildTree(tr)
        else:
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
        if self.parent.opts.model.type == 'codon':
            self.parent.phycassert(not self.parent.opts.model.use_flex_model, 'Cannot currently use flex model within codon model')
            self.parent.phycassert(self.parent.opts.model.num_rates == 1, 'Cannot currently use gamma rate heterogeneity within codon model')
            self.parent.phycassert(not self.parent.opts.model.pinvar_model, 'Cannot currently use invariable sites model within codon model')
            self.model = Likelihood.CodonModel()
            self.model.setKappa(self.parent.opts.model.kappa)
            if self.parent.opts.model.fix_kappa:
                self.model.fixKappa()
            self.model.setOmega(self.parent.opts.model.omega)
            if self.parent.opts.model.fix_omega:
                self.model.fixOmega()
            if not self.parent.opts.model.update_freqs_separately:
                ndimensions = self.parent.opts.model.state_freq_prior.getNParams() 
                self.parent.phycassert(ndimensions == 61, 'state_freq_prior should be 61-dimensional, but instead has %d dimensions' % ndimensions)
            self.parent.phycassert(self.parent.opts.model.state_freqs, 'state_freqs is None, but should be a list containing 61 (unnormalized) relative codon frequencies')
            self.parent.phycassert(len(self.parent.opts.model.state_freqs) == 61, 'state_freqs should be a list containing exactly 61 codon frequencies; instead, it contains %d values' % len(self.parent.opts.model.state_freqs))
            for c in range(61):
                self.parent.phycassert(self.parent.opts.model.state_freqs[c] >= 0.0, 'state_freqs[%d] cannot be negative (%f was specified)' % (c,self.parent.opts.model.state_freqs[c]))
            self.model.setStateFreqsUnnorm(self.parent.opts.model.state_freqs)
            if self.parent.opts.model.fix_freqs:
                self.model.fixStateFreqs()
        elif self.parent.opts.model.type in ['gtr','hky']:
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
            self.parent.phycassert(self.parent.opts.model.state_freqs, 'state_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.parent.phycassert(len(self.parent.opts.model.state_freqs) == 4, 'state_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(self.parent.opts.model.state_freqs))
            self.parent.phycassert(self.parent.opts.model.state_freqs[0] >= 0.0, 'state_freqs[0] cannot be negative (%f was specified)' % self.parent.opts.model.state_freqs[0])
            self.parent.phycassert(self.parent.opts.model.state_freqs[1] >= 0.0, 'state_freqs[1] cannot be negative (%f was specified)' % self.parent.opts.model.state_freqs[1])
            self.parent.phycassert(self.parent.opts.model.state_freqs[2] >= 0.0, 'state_freqs[2] cannot be negative (%f was specified)' % self.parent.opts.model.state_freqs[2])
            self.parent.phycassert(self.parent.opts.model.state_freqs[3] >= 0.0, 'state_freqs[3] cannot be negative (%f was specified)' % self.parent.opts.model.state_freqs[3])
            self.model.setNucleotideFreqs(self.parent.opts.model.state_freqs[0], self.parent.opts.model.state_freqs[1], self.parent.opts.model.state_freqs[2], self.parent.opts.model.state_freqs[3])  #POL should be named setStateFreqs?
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
            
        # Create the object
        self.likelihood = Likelihood.TreeLikelihood(self.model)
        self.likelihood.setLot(self.r)
        self.likelihood.setUFNumEdges(self.parent.opts.uf_num_edges)
        self.likelihood.useUnimap(self.parent.opts.use_unimap)
        if self.parent.data_matrix:
            if POLPY_NEWWAY:
                from phycas import partition
                self.likelihood.copyDataFromDiscreteMatrix(self.parent.data_matrix, partition.getSiteModelVector())
            else:
                self.likelihood.copyDataFromDiscreteMatrix(self.parent.data_matrix)
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
        self.boldness                = 0.0
        self.heating_power           = power
        self.relrate_prior           = cloneDistribution(self.parent.opts.model.relrate_prior)
        self.relrate_param_prior     = cloneDistribution(self.parent.opts.model.relrate_param_prior)
        self.state_freq_prior        = cloneDistribution(self.parent.opts.model.state_freq_prior)
        self.state_freq_param_prior  = cloneDistribution(self.parent.opts.model.state_freq_param_prior)
        self.gamma_shape_prior       = cloneDistribution(self.parent.opts.model.gamma_shape_prior)
        self.external_edgelen_prior  = cloneDistribution(self.parent.opts.model.external_edgelen_prior)
        self.internal_edgelen_prior  = cloneDistribution(self.parent.opts.model.internal_edgelen_prior)
        self.kappa_prior             = cloneDistribution(self.parent.opts.model.kappa_prior)
        self.omega_prior             = cloneDistribution(self.parent.opts.model.omega_prior)
        self.pinvar_prior            = cloneDistribution(self.parent.opts.model.pinvar_prior)
        self.flex_prob_param_prior   = cloneDistribution(self.parent.opts.model.flex_prob_param_prior)
        self.edgelen_hyperprior      = None
        if self.parent.opts.model.edgelen_hyperprior is not None:
            self.edgelen_hyperprior  = cloneDistribution(self.parent.opts.model.edgelen_hyperprior)
        self.chain_manager           = None

        self.setupChain()

    def resetNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the resetNumLikelihoodEvals function of self.likelihood. This
        resets the number of likelihood evaluations performed to 0.

        """
        return self.likelihood.resetNumLikelihoodEvals()
    
    def getNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the getNumLikelihoodEvals function of self.likelihood. This 
        returns the number of likelihood evaluations performed since the last
        call of resetNumLikelihoodEvals.

        """
        return self.likelihood.getNumLikelihoodEvals()
        
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
        if self.parent.opts.doing_steppingstone_sampling:
            paramf.write('Gen\tbeta\tLnL')
        else:
            paramf.write('Gen\tLnL')
            
        if self.parent.opts.fix_topology:
            nbrlens = 2*self.parent.ntax - 3    # will need to be changed if rooted trees allowed
            for i in range(nbrlens):
                paramf.write('\tbrlen%d' % (i+1))
        else:
            paramf.write('\tTL')
        
        paramf.write(self.model.paramHeader())
        if self.parent.opts.model.edgelen_hyperprior is not None:
            if self.parent.opts.model.internal_edgelen_prior is self.parent.opts.model.external_edgelen_prior:
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

    def setBoldness(self, boldness):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the boldness data member and calls the setBoldness method for
        every updater so that all updaters are immediately informed that the
        boldness for this chain has changed. The boldness is a value from 0 to
        100 that specifies the boldness of Metropolis-Hastings moves. Setting
        the boldness has no effect on slice sampling based updaters. Each 
        move class defines what is meant by boldness. The boldness is changed
        during an MCMC run in some circumstances, such as during a path 
        sampling analysis where the target distribution changes during the 
        run.
        
        """
        self.boldness = boldness
        for updater in self.chain_manager.getAllUpdaters():
            updater.setBoldness(boldness)

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
        self.kappa_prior.setLot(self.r)
        self.relrate_prior.setLot(self.r)
        self.state_freq_prior.setLot(self.r)
        self.state_freq_param_prior.setLot(self.r)
        self.flex_prob_param_prior.setLot(self.r)
        self.gamma_shape_prior.setLot(self.r)
        self.pinvar_prior.setLot(self.r)
        if self.edgelen_hyperprior is not None:
            self.edgelen_hyperprior.setLot(self.r)
        if self.parent.opts.model.external_edgelen_prior is not None:
            self.external_edgelen_prior.setLot(self.r)
        if self.parent.opts.model.internal_edgelen_prior is not None:
            self.internal_edgelen_prior.setLot(self.r)
        
        # Define priors for the model parameters
        if self.parent.opts.model.type == 'codon':
            self.model.setKappaPrior(self.kappa_prior)
            self.model.setOmegaPrior(self.omega_prior)
            if self.parent.opts.model.update_freqs_separately:
                self.model.setStateFreqParamPrior(self.state_freq_param_prior)
            else:
                self.model.setStateFreqPrior(self.state_freq_prior)
        elif self.parent.opts.model.type == 'gtr':
            if self.parent.opts.model.update_relrates_separately:
                self.model.setRelRateParamPrior(self.relrate_param_prior)
            else:
                self.model.setRelRatePrior(self.relrate_prior)
            if self.parent.opts.model.update_freqs_separately:
                self.model.setStateFreqParamPrior(self.state_freq_param_prior)
            else:
                self.model.setStateFreqPrior(self.state_freq_prior)
        elif self.parent.opts.model.type == 'hky':
            self.model.setKappaPrior(self.kappa_prior)
            if self.parent.opts.model.update_freqs_separately:
                self.model.setStateFreqParamPrior(self.state_freq_param_prior)
            else:
                self.model.setStateFreqPrior(self.state_freq_prior)

        # If rate heterogeneity is to be assumed, add priors for these model parameters here
        if self.parent.opts.model.use_flex_model:
            self.model.setNumFlexSpacers(self.parent.opts.model.flex_num_spacers)
            self.model.setFLEXProbParamPrior(self.flex_prob_param_prior)
        elif self.parent.opts.model.num_rates > 1:
            self.model.setDiscreteGammaShapePrior(self.gamma_shape_prior)
        if self.parent.opts.model.pinvar_model:
            self.model.setPinvarPrior(self.pinvar_prior)
        
        # Define edge length prior distributions
        separate_edge_len_dists = self.parent.opts.model.internal_edgelen_prior is not self.parent.opts.model.external_edgelen_prior
        self.model.separateInternalExternalEdgeLenPriors(separate_edge_len_dists)
        self.model.setExternalEdgeLenPrior(self.external_edgelen_prior)
        self.model.setInternalEdgeLenPrior(self.internal_edgelen_prior)

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
            #False,                         # separate_edgelen_params (deprecated: always False)
            self.parent.opts.slice_max_units,    # maximum number of slice units allowed
            self.parent.opts.slice_weight)       # weight for each parameter added
            
        # Create a TreeScalerMove object to handle scaling the entire tree to allow faster
        # convergence in edge lengths. This move is unusual in using slice sampling rather
        # than Metropolis-Hastings updates: most "moves" in parent are Metropolis-Hastings.
        if self.parent.opts.tree_scaler_weight > 0:
            self.tree_scaler_move = Likelihood.TreeScalerMove()
            self.tree_scaler_move.setName("Tree scaler move")
            self.tree_scaler_move.setWeight(self.parent.opts.tree_scaler_weight)
            self.tree_scaler_move.setPosteriorTuningParam(self.parent.opts.tree_scaler_lambda)
            self.tree_scaler_move.setPriorTuningParam(self.parent.opts.tree_scaler_lambda0)
            self.tree_scaler_move.setTree(self.tree)
            self.tree_scaler_move.setModel(self.model)
            self.tree_scaler_move.setTreeLikelihood(self.likelihood)
            self.tree_scaler_move.setLot(self.r)
            if self.model.edgeLengthsFixed():
                self.tree_scaler_move.fixParameter()
            self.chain_manager.addMove(self.tree_scaler_move)

        if not self.parent.opts.model.type == 'jc' and not self.parent.opts.model.update_freqs_separately:
            # Create a StateFreqMove to update entire state frequency vector
            self.state_freq_move = Likelihood.StateFreqMove()
            if self.parent.opts.model.type == 'codon': 
                self.state_freq_move.setDimension(61)
            else:
                self.state_freq_move.setDimension(4)
            self.state_freq_move.setName("State freq move")
            self.state_freq_move.setWeight(self.parent.opts.state_freq_weight)
            self.state_freq_move.setPosteriorTuningParam(self.parent.opts.state_freq_psi)
            self.state_freq_move.setPriorTuningParam(self.parent.opts.state_freq_psi0)
            self.state_freq_move.setWeight(self.parent.opts.state_freq_weight)
            self.state_freq_move.setTree(self.tree)
            self.state_freq_move.setModel(self.model)
            self.state_freq_move.setTreeLikelihood(self.likelihood)
            self.state_freq_move.setLot(self.r)
            if self.model.stateFreqsFixed():
                self.state_freq_move.fixParameter()
            self.state_freq_move.setMultivarPrior(self.state_freq_prior)
            self.chain_manager.addMove(self.state_freq_move)
        
        if self.parent.opts.model.type == 'gtr' and not self.parent.opts.model.update_relrates_separately:
            # Create a RelRateMove to update entire relative rates vector
            self.rel_rate_move = Likelihood.RelRatesMove()
            self.rel_rate_move.setName("Relative rates move")
            self.rel_rate_move.setWeight(self.parent.opts.rel_rate_weight)
            self.rel_rate_move.setPosteriorTuningParam(self.parent.opts.rel_rate_psi)
            self.rel_rate_move.setPriorTuningParam(self.parent.opts.rel_rate_psi0)
            self.rel_rate_move.setTree(self.tree)
            self.rel_rate_move.setModel(self.model)
            self.rel_rate_move.setTreeLikelihood(self.likelihood)
            self.rel_rate_move.setLot(self.r)
            #if self.model.relRatesFixed():
            #    self.rel_rate_move.fixParameter()
            self.rel_rate_move.setMultivarPrior(self.relrate_prior)
            self.chain_manager.addMove(self.rel_rate_move)
        
        if self.parent.opts.fix_topology:
            # Create an EdgeMove object to handle Metropolis-Hastings
            # updates to the edge lengths only (does not change the topology)
            self.edge_move = Likelihood.EdgeMove()
            self.edge_move.setName("Edge length move")
            self.edge_move.setWeight(self.parent.opts.edge_move_weight)
            self.edge_move.setPosteriorTuningParam(self.parent.opts.edge_move_lambda)
            self.edge_move.setPriorTuningParam(self.parent.opts.edge_move_lambda0)
            self.edge_move.setTree(self.tree)
            self.edge_move.setModel(self.model)
            self.edge_move.setTreeLikelihood(self.likelihood)
            self.edge_move.setLot(self.r)
            self.edge_move.setLambda(self.parent.opts.edge_move_lambda)
            if self.model.edgeLengthsFixed():
                self.edge_move.fixParameter()
            self.chain_manager.addMove(self.edge_move)

        if self.parent.opts.use_unimap:
            # Create a UnimapFastNNIMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_fast_nni_move = Likelihood.UnimapFastNNIMove()
            self.unimap_fast_nni_move.setName("Unimap Fast NNI move")
            self.unimap_fast_nni_move.setWeight(self.parent.opts.unimap_fast_nni_move_weight)
            self.unimap_fast_nni_move.setTree(self.tree)
            self.unimap_fast_nni_move.setModel(self.model)
            self.unimap_fast_nni_move.setTreeLikelihood(self.likelihood)
            self.unimap_fast_nni_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_fast_nni_move)

            # Create a UnimapSampleAmbigMove
            wt = self.parent.opts.unimap_sample_ambig_move_weight
            self.unimap_sample_ambig_move = Likelihood.UnimapSampleAmbigMove(self.likelihood, self.tree, self.model, wt)
            self.unimap_sample_ambig_move.setName("Unimap Sample Ambig move")
            num_ambig = self.unimap_sample_ambig_move.getNumAmbigNodes()
            if wt > 0.0 and num_ambig > 0:
                self.unimap_sample_ambig_move.setLot(self.r)
                self.chain_manager.addMove(self.unimap_sample_ambig_move)

            # Create a UnimapNNIMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_nni_move = Likelihood.UnimapNNIMove()
            self.unimap_nni_move.setName("Unimap NNI move")
            self.unimap_nni_move.setWeight(self.parent.opts.unimap_nni_move_weight)
            self.unimap_nni_move.setTree(self.tree)
            self.unimap_nni_move.setModel(self.model)
            self.unimap_nni_move.setTreeLikelihood(self.likelihood)
            self.unimap_nni_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_nni_move)

            # Create a UnimapNodeSlideMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_node_slide_move = Likelihood.UnimapNodeSlideMove()
            self.unimap_node_slide_move.setName("Unimap NodeSlide move")
            self.unimap_node_slide_move.setWeight(self.parent.opts.unimap_node_slide_move_weight)
            self.unimap_node_slide_move.setTree(self.tree)
            self.unimap_node_slide_move.setModel(self.model)
            self.unimap_node_slide_move.setTreeLikelihood(self.likelihood)
            self.unimap_node_slide_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_node_slide_move)

            # Create a MappingMove (to refresh the mapping for all sites)
            self.unimapping_move = Likelihood.MappingMove()
            self.unimapping_move.setName("Univent mapping move")
            self.unimapping_move.setWeight(self.parent.opts.mapping_move_weight)
            self.unimapping_move.setTree(self.tree)
            self.unimapping_move.setModel(self.model)
            self.unimapping_move.setTreeLikelihood(self.likelihood)
            self.unimapping_move.setLot(self.r)
            self.chain_manager.addMove(self.unimapping_move)
            
            # Create a UnimapEdgeMove object to handle Metropolis-Hastings
            # updates to the edge lengths only (does not change the topology)
            self.unimap_edge_move = Likelihood.UnimapEdgeMove()
            self.unimap_edge_move.setName("Unimap edge length move")
            self.unimap_edge_move.setWeight(self.parent.opts.unimap_edge_move_weight)
            self.unimap_edge_move.setPosteriorTuningParam(self.parent.opts.unimap_edge_move_lambda)
            self.unimap_edge_move.setPriorTuningParam(self.parent.opts.unimap_edge_move_lambda0)
            self.unimap_edge_move.setTree(self.tree)
            self.unimap_edge_move.setModel(self.model)
            self.unimap_edge_move.setTreeLikelihood(self.likelihood)
            self.unimap_edge_move.setLot(self.r)
            self.unimap_edge_move.setLambda(self.parent.opts.unimap_edge_move_lambda)
            if self.model.edgeLengthsFixed():
                self.unimap_edge_move.fixParameter()
            self.chain_manager.addMove(self.unimap_edge_move)

        elif not self.parent.opts.fix_topology:
            # Create a LargetSimonMove object to handle Metropolis-Hastings
            # updates to the tree topology and edge lengths
            self.larget_simon_move = Likelihood.LargetSimonMove()
            self.larget_simon_move.setName("Larget-Simon move")
            self.larget_simon_move.setWeight(self.parent.opts.ls_move_weight)
            self.larget_simon_move.setPosteriorTuningParam(self.parent.opts.ls_move_lambda)
            self.larget_simon_move.setPriorTuningParam(self.parent.opts.ls_move_lambda0)
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

        if self.parent.opts.use_unimap:
            if num_ambig > 0:
                self.unimap_sample_ambig_move.sampleTipsAsDisconnected()
            self.likelihood.fullRemapping(self.tree, self.r, True) 

        # Calculate relative rates if among-site rate heterogeneity is part of the model
        # Currently, the unnormalized rates are stored in the model; the function call
        # below normalizes the rate means and probabilities so that they will be accurately
        # recorded in the first line of the parameter file
        self.likelihood.recalcRelativeRates()

        # Make sure each updater knows the heating power and heating type
        # and set boldness to 0 for each updater (0% boldness is the default,
        # with more boldness used only in steppingstone sampling where bolder
        # moves are needed as the influence of the likelihood is progressively
        # diminished)
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(self.heating_power)
            updater.setBoldness(0.0)
            if self.parent.opts.ss_heating_likelihood:
                updater.setLikelihoodHeating()
            else:
                updater.setStandardHeating()

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
        self.swap_table = None

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
    
        # Create the chains
        # print '@@@@@ heat vector =',self.parent.heat_vector
        for i, heating_power in enumerate(self.parent.heat_vector):
            markov_chain = MarkovChain(self.parent, heating_power)
            self.chains.append(markov_chain)
            
        n = len(self.parent.heat_vector)
        self.swap_table = [[0]*n for i in range(n)]

    def getNumChains(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of chains being managed by this object
        (i.e. returns the length of the self.chains list).
        
        """
        return len(self.chains)

    def getColdChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the cold chain MarkovChain object. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        i = self.parent.heat_vector.index(1.0)
        return self.chains[0]

    def getColdChainManager(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the chain manager for the cold chain. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        return self.chains[0].chain_manager

    def setChainPower(self, chain_index, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Loops through all updaters in the chain having the specified 
        chain_index, setting the power for each updater to the supplied power.
        
        """
        self.parent.phycassert(len(self.chains) > chain_index, 'chain index specified (%d) too large for number of chains (%d)' % (chain_index, len(self.chains)))
        self.parent.phycassert(len(self.parent.heat_vector) == len(self.chains), 'length of heat vector (%d) not equal to number of chains (%d)' % (len(self.parent.heat_vector), len(self.chains)))
        self.chains[chain_index].setPower(power)

    def resetNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls resetNumLikelihoodEvals method for every MarkovChain in the
        self.chains list.
        
        """
        for c in self.chains:
            c.resetNumLikelihoodEvals()

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
        MarkovChain objects in the self.chains list by calling the 
        getNumLikelihoodEvals method for each chain.
        
        """
        total = 0
        for c in self.chains:
            total += c.getNumLikelihoodEvals()
        return total

    def recordSample(self, cycle = -1):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Records the current tree topology and edge lengths by adding a line to
        the tree file, and records tree length and substitution parameters
        by adding a line to the parameter file.
        
        """
        float_format_str = '%%.%df\t' % self.parent.opts.ndecimals
        float_format_notab_str = '%%.%df' % self.parent.opts.ndecimals
        
        # Gather log-likelihoods, and if path sampling save in path_sample list for later
        lnLikes = []
        for i,c in enumerate(self.chains):
            lnLi = c.chain_manager.getLastLnLike()
            lnLikes.append(lnLi)
        
        # Only record samples from the current cold chain
        cold_chain = self.parent.mcmc_manager.getColdChain()
        
        # Add line to parameter file if it exists
        if self.parent.paramf:
            self.parent.paramf.write('%d\t' % (cycle + 1))
            if self.parent.opts.doing_steppingstone_sampling:
                self.parent.paramf.write(float_format_str % (cold_chain.heating_power))
            for lnl in lnLikes:
                self.parent.paramf.write(float_format_str % lnl)
            
            # TREE_LENGTH
            if self.parent.opts.fix_topology:
                self.parent.paramf.write('%s\t' % '\t'.join([float_format_notab_str % brlen for brlen in cold_chain.tree.edgeLens()]))
            else:
                self.parent.paramf.write(float_format_str % cold_chain.tree.edgeLenSum())
            
            self.parent.paramf.write(cold_chain.model.paramReport(self.parent.opts.ndecimals))
            if self.parent.opts.model.edgelen_hyperprior is not None:
                for p in cold_chain.chain_manager.getEdgeLenHyperparams():
                    self.parent.paramf.write(float_format_str % p.getCurrValue())
            if self.parent.opts.model.use_flex_model:
                rates_vector = cold_chain.likelihood.getRateMeans()
                for rr in rates_vector:
                    self.parent.paramf.write(float_format_str % rr)
                probs_vector = cold_chain.likelihood.getRateProbs()
                for rp in probs_vector:
                    self.parent.paramf.write(float_format_str % rp)
            self.parent.paramf.write('\n')
            self.parent.paramf.flush()
            
        # Add line to tree file if it exists
        if self.parent.treef:
            self.parent.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, cold_chain.tree.makeNewick(self.parent.opts.ndecimals)))
            self.parent.treef.flush()

        # Add line to sitelike file if it exists and if we are saving site-likelihoods
        if self.parent.opts.saving_sitelikes:
            if self.parent.sitelikef:
                cold_chain.likelihood.storeSiteLikelihoods(True)
                cold_chain.likelihood.calcLnL(cold_chain.tree)
                cold_chain.likelihood.storeSiteLikelihoods(False)
                
                patternLnLikes = cold_chain.likelihood.getSiteLikelihoods()
                siteloglikes = [0.0]*self.parent.nchar
                for i,patternLnL in enumerate(patternLnLikes):
                    for site in self.parent.siteIndicesForPatternIndex[i]:
                        siteloglikes[site] = patternLnL
                for siteLnL in siteloglikes:
                    self.parent.sitelikef.write(float_format_str % siteLnL)
                self.parent.sitelikef.write('\n')
                self.parent.sitelikef.flush()
                
    def attemptChainSwap(self, cycle):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Attempts to swap two randomly-chosen chains. Letting
        p_i be posterior density for chain i
        p_j be posterior density for chain j
        theta_i be the parameter vector for chain i
        theta_j be the parameter vector for chain j
        
                       p_i(theta_j)  p_j(theta_i)
        accept_ratio = -----------   ------------
                       p_i(theta_i)  p_j(theta_j)
        
          L(theta_j)^power_i f(theta_j)  L(theta_i)^power_j f(theta_i) 
        = -----------------------------  -----------------------------
          L(theta_i)^power_i f(theta_i)  L(theta_j)^power_j f(theta_j) 
        
        = L(theta_i)^(power_j - power_i) L(theta_j)^(power_i - power_j)
        
        """
        n = len(self.chains)
        if n < 2:
            return
            
        # Figure out which chains to swap
        
        r = self.parent._getLot()
        i = r.sampleUInt(n)
        j = r.sampleUInt(n - 1)
        if j >= i:
            j += 1
        else:
            i,j = j,i   # make sure i < j so attempted swaps are consistently in upper triangle of swap table
        
        # Compute acceptance ratio
        cmi = self.chains[i].chain_manager
        power_i = self.chains[i].heating_power
        cmi.refreshLastLnLike()
        lnLi = cmi.getLastLnLike()
        cmi.refreshLastLnPrior()
        lnPriori = cmi.getLastLnPrior()

        cmj = self.chains[j].chain_manager
        power_j = self.chains[j].heating_power
        cmj.refreshLastLnLike()
        lnLj = cmj.getLastLnLike()
        cmj.refreshLastLnPrior()
        lnPriorj = cmj.getLastLnPrior()

        log_accept_ratio = (power_j - power_i)*(lnLi + lnPriori - lnLj - lnPriorj)
        log_u = math.log(r.uniform())
        swap_accepted = False
        self.swap_table[i][j] += 1  # upper triangle
        if log_u < log_accept_ratio:
            swap_accepted = True
            self.setChainPower(i, power_j)
            self.setChainPower(j, power_i)
            self.chains[i], self.chains[j] = self.chains[j], self.chains[i]
            
            self.swap_table[j][i] += 1  # lower triangle
            
        if cycle % 1 == 0:
            print '\n********** Attempting chain swap **********'
            print 'n =',n
            print 'i =',i
            print 'j =',j
            print 'log_accept_ratio =',log_accept_ratio
            print 'log_u            =',log_u
            if swap_accepted:
                print 'swap_accepted    = True'
            else:
                print 'swap_accepted    = False'
            print 'swap_table:'
            for ii in range(n):
                for jj in range(n):
                    print '%12d' % self.swap_table[ii][jj],
                print
            print '*********************************************\n'
