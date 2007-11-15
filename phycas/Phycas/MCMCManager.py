import os,sys,math
from phycas import *

def cloneDistribution(d):
    if d == None:
        return None
    else:
        return d.clone()

class LikelihoodCore:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Not yet documented.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Not yet documented.
        
        """
        self.phycas                             = phycas
        self.model                              = None
        self.likelihood                         = None
        self.tree                               = Phylogeny.Tree()
        self.r                                  = ProbDist.Lot()
        self.starting_edgelen_dist              = cloneDistribution(self.phycas.starting_edgelen_dist)
    
    def setupCore(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Not yet documented.
        
        """
        # Set seed if user has supplied one
        if self.phycas.random_seed != 0:
            self.r.setSeed(int(self.phycas.random_seed))

        self.starting_edgelen_dist.setLot(self.r)

        # Create a substitution model
        if self.phycas.default_model in ['gtr','hky']:
            if self.phycas.default_model == 'gtr':
                self.model = Likelihood.GTRModel()
                self.model.setRelRates(self.phycas.starting_relrates)
                if self.phycas.fix_relrates:
                    self.model.fixRelRates()
            else:
                self.model = Likelihood.HKYModel()
                self.model.setKappa(self.phycas.starting_kappa)
                if self.phycas.fix_kappa:
                    self.model.fixKappa()
            self.phycas.phycassert(self.phycas.starting_freqs, 'starting_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.phycas.phycassert(len(self.phycas.starting_freqs) == 4, 'starting_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(self.phycas.starting_freqs))
            self.phycas.phycassert(self.phycas.starting_freqs[0] >= 0.0, 'starting_freqs[0] cannot be negative (%f was specified)' % self.phycas.starting_freqs[0])
            self.phycas.phycassert(self.phycas.starting_freqs[1] >= 0.0, 'starting_freqs[1] cannot be negative (%f was specified)' % self.phycas.starting_freqs[1])
            self.phycas.phycassert(self.phycas.starting_freqs[2] >= 0.0, 'starting_freqs[2] cannot be negative (%f was specified)' % self.phycas.starting_freqs[2])
            self.phycas.phycassert(self.phycas.starting_freqs[3] >= 0.0, 'starting_freqs[3] cannot be negative (%f was specified)' % self.phycas.starting_freqs[3])
            self.model.setNucleotideFreqs(self.phycas.starting_freqs[0], self.phycas.starting_freqs[1], self.phycas.starting_freqs[2], self.phycas.starting_freqs[3])  #POL should be named setStateFreqs?
            if self.phycas.fix_freqs:
                self.model.fixStateFreqs()
        else:
            self.model = Likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
        # Note must defer setting up pattern specific rates model until we know number of patterns
        if self.phycas.use_flex_model:
            self.model.setNGammaRates(self.phycas.num_rates)
            self.model.setFlexModel()
            self.model.setFlexRateUpperBound(self.phycas.flex_L)
        elif self.phycas.num_rates > 1:
            self.model.setNGammaRates(self.phycas.num_rates)
            self.model.setPriorOnShapeInverse(self.phycas.use_inverse_shape)    #POL should be named useInverseShape rather than setPriorOnShapeInverse
            self.model.setShape(self.phycas.starting_shape)
            if self.phycas.fix_shape:
                self.model.fixShape()
        else:
            self.model.setNGammaRates(1)
            
        if self.phycas.estimate_pinvar:
            assert not self.phycas.use_flex_model, 'Cannot currently use flex model with pinvar'
            self.model.setPinvarModel()
            self.model.setPinvar(self.phycas.starting_pinvar)
            if self.phycas.fix_pinvar:
                self.model.fixPinvar()
        else:
            self.model.setNotPinvarModel()

        if self.phycas.fix_edgelens:
            self.model.fixEdgeLengths()
            
        # Create the likelihood object (formerly the Phycas.setupLikelihood function)
        self.likelihood = Likelihood.TreeLikelihood(self.model)
        self.likelihood.setUFNumEdges(self.phycas.uf_num_edges)
        if self.phycas.data_source == 'file':
            self.likelihood.copyDataFromDiscreteMatrix(self.phycas.data_matrix)
            self.npatterns = self.likelihood.getNPatterns()
        elif not self.phycas.data_source:
            self.phycas.phycassert(self.phycas.ntax > 0, 'data_source is None, which indicates data will be simulated, but in this case sim_taxon_labels should not be an empty list')

        # Build the starting tree
        if self.phycas.starting_tree == None:
            # Build a random tree
            Phylogeny.TreeManip(self.tree).randomTree(
                self.phycas.ntax,           # number of tips
                self.r,                     # pseudorandom number generator
                self.starting_edgelen_dist, # distribution from which to draw starting edge lengths
                False)                      # Yule tree if True, edge lengths independent if False
            self.phycas.warn_tip_numbers = False
            self.phycas.starting_tree = self.tree.makeNewick()
        else:
            # Build user-specified tree
            self.tree.buildFromString(self.phycas.starting_tree)
            if not self.tree.tipNumbersSetUsingNames():
                self.phycas.warn_tip_numbers = True

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
        self.phycas.phycassert(self.phycas.sim_nchar > 0, 'sim_nchar must be greater than zero in order to perform simulations')
        self.likelihood.simulateFirst(sim_data, self.tree, self.r, self.phycas.sim_nchar)
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
    Not yet documented.
    
    """
    def __init__(self, phycas, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Not yet documented.
        
        """
        LikelihoodCore.__init__(self, phycas)
        
        self.phycas                             = phycas
        self.heating_power                      = power
        self.relrate_prior                      = cloneDistribution(self.phycas.relrate_prior)
        self.base_freq_param_prior              = cloneDistribution(self.phycas.base_freq_param_prior)
        self.gamma_shape_prior                  = cloneDistribution(self.phycas.gamma_shape_prior)
        self.edgelen_hyperprior                 = cloneDistribution(self.phycas.edgelen_hyperprior)
        self.external_edgelen_dist              = cloneDistribution(self.phycas.external_edgelen_dist)
        self.internal_edgelen_dist              = cloneDistribution(self.phycas.internal_edgelen_dist)
        self.kappa_prior                        = cloneDistribution(self.phycas.kappa_prior)
        self.pinvar_prior                       = cloneDistribution(self.phycas.pinvar_prior)
        self.flex_prob_param_prior              = cloneDistribution(self.phycas.flex_prob_param_prior)
        self.chain_manager                      = None

        self.setupChain()

    def resetNEvals(self):
        return self.likelihood.resetNEvals()
    
    def getNEvals(self):
        return self.likelihood.getNEvals()
        
    def paramFileHeader(self, paramf):
        paramf.write('[ID: %d]\n' % self.r.getInitSeed())
        paramf.write(self.model.paramHeader())
        if self.phycas.using_hyperprior:
            if self.phycas.internal_edgelen_dist is self.phycas.external_edgelen_dist:
                paramf.write('\thyper(all)')
            else:
                paramf.write('\thyper(external)')
                paramf.write('\thyper(internal)')
        if self.phycas.use_flex_model:
            paramf.write('\trates_probs')

    def treeFileHeader(self, treef):
        treef.write('#NEXUS\n')
        treef.write('[ID: %d]\n' % self.r.getInitSeed())
        treef.write('begin trees;\n')
        treef.write('   translate\n')
        for i in range(self.phycas.ntax):
            if self.phycas.taxon_labels[i].find(' ') < 0:
                # no spaces found in name
                treef.write('       %d %s%s\n' % (i + 1, self.phycas.taxon_labels[i], i == self.phycas.ntax - 1 and ';' or ','))
            else:
                # at least one space in taxon name, so enclose name in quotes
                treef.write("       %d '%s'%s\n" % (i + 1, self.phycas.taxon_labels[i], i == self.phycas.ntax - 1 and ';' or ','))

    def setPower(self, power):
        self.heating_power = power
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(power)

    def setupChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Not yet documented.
        
        """
        LikelihoodCore.setupCore(self)
        LikelihoodCore.prepareForLikelihood(self)
        
        # Make sure that each is using the same pseudorandom number generator object
        self.relrate_prior.setLot(self.r)
        self.base_freq_param_prior.setLot(self.r)
        self.flex_prob_param_prior.setLot(self.r)
        self.gamma_shape_prior.setLot(self.r)
        self.pinvar_prior.setLot(self.r)
        self.edgelen_hyperprior.setLot(self.r)
        self.external_edgelen_dist.setLot(self.r)
        if self.phycas.internal_edgelen_dist:
            self.internal_edgelen_dist.setLot(self.r)
        
        # define priors for the model parameters
        if self.phycas.default_model == 'gtr':
            self.model.setRelRatePrior(self.relrate_prior)
            self.model.setStateFreqParamPrior(self.base_freq_param_prior)   #POL should be named state_freq_param_prior
        elif self.phycas.default_model == 'hky':
            self.model.setKappaPrior(self.kappa_prior)
            self.model.setStateFreqParamPrior(self.base_freq_param_prior)   #POL should be named state_freq_param_prior

        # If rate heterogeneity is to be assumed, add priors for these model parameters here
        if self.phycas.use_flex_model:
            self.model.setNumFlexSpacers(self.phycas.flex_num_spacers)
            self.model.setFLEXProbParamPrior(self.flex_prob_param_prior)
        elif self.phycas.num_rates > 1:
            self.model.setDiscreteGammaShapePrior(self.gamma_shape_prior)
        if self.phycas.estimate_pinvar:
            self.model.setPinvarPrior(self.pinvar_prior)
        
        # Define edge length prior distributions
        separate_edge_len_dists = self.phycas.internal_edgelen_dist is not self.phycas.internal_edgelen_dist
        self.model.separateInternalExternalEdgeLenPriors(separate_edge_len_dists)
        self.model.setExternalEdgeLenPrior(self.external_edgelen_dist)
        self.model.setInternalEdgeLenPrior(self.internal_edgelen_dist)

        if self.phycas.using_hyperprior:
            #self.edgelen_hyperprior.setMeanAndVariance(1.0, 10.0)
            self.model.setEdgeLenHyperPrior(self.edgelen_hyperprior)
            #todo self.model.starting_edgelen_hyperparam
            if self.phycas.fix_edgelen_hyperparam:
                self.model.fixEdgeLenHyperprior()   #POL should be named fixEdgeLenHyperparam
        else:
            self.model.setEdgeLenHyperPrior(None)

        self.likelihood.replaceModel(self.model)            

        if self.phycas.data_source == None:
            self.likelihood.setNoData() # user wants to run MCMC with no data
        
        # Create an MCMCChainManager object and add all necessary updaters
        self.chain_manager = Likelihood.MCMCChainManager()
        self.chain_manager.addMCMCUpdaters(
            self.model,              # substitution model
            self.tree,               # tree
            self.likelihood,         # likelihood calculation machinery
            self.r,                  # pseudorandom number generator
            #POLPY_NEWWAY False,     # separate_edgelen_params (deprecated: always False)
            self.phycas.slice_max_units,    # maximum number of slice units allowed
            self.phycas.slice_weight)       # weight for each parameter added

        # Create a TreeScalerMove object to handle scaling the entire tree to allow faster
        # convergence in edge lengths. This move is unusual in using slice sampling rather
        # than Metropolis-Hastings updates: most "moves" in Phycas are Metropolis-Hastings.
        if self.phycas.tree_scaler_weight > 0:
            self.tree_scaler_move = Likelihood.TreeScalerMove()
            self.tree_scaler_move.setName("Tree scaler move")
            self.tree_scaler_move.setWeight(self.phycas.tree_scaler_weight)
            self.tree_scaler_move.setTree(self.tree)
            self.tree_scaler_move.setModel(self.model)
            self.tree_scaler_move.setTreeLikelihood(self.likelihood)
            self.tree_scaler_move.setLot(self.r)
            if self.model.edgeLengthsFixed():
                self.tree_scaler_move.fixParameter()
            self.chain_manager.addMove(self.tree_scaler_move)

        # Create a LargetSimonMove object to handle Metropolis-Hastings
        # updates to the tree topology and edge lengths
        self.larget_simon_move = Likelihood.LargetSimonMove()
        self.larget_simon_move.setName("Larget-Simon move")
        self.larget_simon_move.setWeight(self.phycas.ls_move_weight)
        self.larget_simon_move.setTree(self.tree)
        self.larget_simon_move.setModel(self.model)
        self.larget_simon_move.setTreeLikelihood(self.likelihood)
        self.larget_simon_move.setLot(self.r)
        self.larget_simon_move.setLambda(self.phycas.ls_move_lambda)
        if self.model.edgeLengthsFixed():
            self.larget_simon_move.fixParameter()
        self.chain_manager.addMove(self.larget_simon_move)

        # If requested, create an NCatMove object to allow the number of rate categories to change
        if self.phycas.use_flex_model:
            # Create an NCatMove object
            self.ncat_move = Likelihood.NCatMove()
            
            # Set up features specific to NCatMove
            self.ncat_move.setCatProbPrior(self.flex_prob_param_prior)
            self.ncat_move.setL(self.phycas.flex_L)
            self.ncat_move.setS(self.phycas.flex_num_spacers)
            self.ncat_move.setLambda(self.phycas.flex_lambda)
            self.ncat_move.setPhi(self.phycas.flex_phi)

            # Continue setting up NCatMove object
            self.ncat_move.setName("NCat move")
            self.ncat_move.setWeight(self.phycas.flex_ncat_move_weight)
            self.ncat_move.setTree(self.tree)
            self.ncat_move.setModel(self.model)
            self.ncat_move.setTreeLikelihood(self.likelihood)
            self.ncat_move.setLot(self.r)
            
            self.chain_manager.addMove(self.ncat_move)
            
        # If requested, create a BushMove object to allow polytomous trees
        if self.phycas.allow_polytomies:
            # Create a BushMove object
            self.bush_move = Likelihood.BushMove()

            # Set up the topology prior
            self.topo_prior_calculator = self.bush_move.getTopoPriorCalculator()
            self.topo_prior_calculator.chooseUnrooted()
            self.topo_prior_calculator.setC(self.phycas.topo_prior_C)
            if self.phycas.polytomy_prior:
                self.topo_prior_calculator.choosePolytomyPrior()
            else:
                self.topo_prior_calculator.chooseResolutionClassPrior()
                
            # Continue setting up BushMove object
            self.bush_move.setName("Bush move")
            self.bush_move.setWeight(self.phycas.bush_move_weight)
            self.bush_move.setTree(self.tree)
            self.bush_move.setModel(self.model)
            self.bush_move.setTreeLikelihood(self.likelihood)
            self.bush_move.setLot(self.r)
            self.bush_move.setEdgeLenDistMean(self.phycas.bush_move_edgelen_mean)
            self.bush_move.viewProposedMove(self.phycas.bush_move_debug)
            if self.model.edgeLengthsFixed():
                self.bush_move.fixParameter()
            self.bush_move.finalize()
            
            self.chain_manager.addMove(self.bush_move)

        if self.phycas.doing_samc:
            self.samc_move = Likelihood.SamcMove(self.starting_edgelen_dist)

            # Continue setting up SAMC move object
            self.samc_move.setName("SAMC move")
            self.samc_move.setWeight(self.phycas.samc_move_weight)
            self.samc_move.setTree(self.tree)
            self.samc_move.setModel(self.model)
            self.samc_move.setTreeLikelihood(self.likelihood)
            self.samc_move.setLot(self.r)
            if self.model.edgeLengthsFixed():
                self.samc_move.fixParameter()
            #self.samc_move.finalize()
            self.chain_manager.addMove(self.samc_move)

        self.chain_manager.finalize()

        # Calculate relative rates if among-site rate heterogeneity is part of the model
        # Currently, the unnormalized rates are stored in the model; the function call
        # below normalizes the rate means and probabilities so that they will be accurately
        # recorded in the first line of the parameter file
        self.likelihood.recalcRelativeRates()

        # Make sure each updater knows the heating power and heating type
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(self.heating_power)
            if self.phycas.is_standard_heating:
                updater.setStandardHeating()
            else:
                updater.setLikelihoodHeating()

class MCMCManager:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This class manages one to many Markov chains participating in a
    Bayesian MCMC analysis. It creates the chains and manages swapping
    and sampling the chains. If likelihood heating is employed, it also
    has the ability to estimate the marginal likelihood of the model.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Stores the phycas object passed into the constructor and creates an
        empty chains vector.
        
        """
        self.phycas = phycas
        self.chains = []

    def paramFileHeader(self, paramf):
        m = self.getColdChain()
        m.paramFileHeader(paramf)

    def treeFileHeader(self, treef):
        m = self.getColdChain()
        m.treeFileHeader(treef)

    def createChains(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a separate MarkovChain for every element in phycas.heat_vector
        and, after calling the setupChain function for the MarkovChain, adds
        it to the data member chains.
        
        """
        for heating_power in self.phycas.heat_vector:
            markov_chain = MarkovChain(self.phycas, heating_power)
            self.chains.append(markov_chain)

    def getNumChains(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of chains being managed by this object.
        
        """
        return len(self.chains)

    def getColdChainIndex(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the index (0..nchains-1) of the chain corresponding to the
        cold chain (i.e. the first chain whose power = 1.0).
        
        """
        return self.phycas.heat_vector.index(1.0)

    def getColdChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the cold chain MCMCChain object. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        i = self.phycas.heat_vector.index(1.0)
        return self.chains[i]

    def getColdChainManager(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the chain manager for the cold chain. The cold chain is the
        first chain in the data member chains whose power equals 1.0.
        
        """
        i = self.phycas.heat_vector.index(1.0)
        return self.chains[i].chain_manager

    def setChainPower(self, chain_index, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Loops through all updaters in the chain having the specified index,
        and sets the power for each updater to the supplied power.
        
        """
        phycas.phycassert(len(chains) > chain_index, 'chain index specified (%d) too large for number of chains (%d)' % (chain_index, len(chains)))
        phycas.phycassert(len(self.phycas.heat_vector) == len(self.chains), 'length of heat vector (%d) not equal to number of chains (%d)' % (len(self.heat_vector, len(self.nchains))))
        self.phycas.heat_vector[chain_index] = power
        self.chains[chain_index].setPower(power)

    def resetNEvals(self):
        for c in self.chains:
            c.resetNEvals()

    def setRandomSeedAllChains(self, rnseed):
        self.phycas.random_seed = rnseed
        for c in self.chains:
            c.r.setSeed(int(rnseed))

    def getTotalEvals(self):
        total = 0
        for c in self.chains:
            total += c.getNEvals()
        return total

    def recordSample(self, cycle):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Records current tree topology and edge lengths by adding a line to
        the tree file, and records tree length and substitution parameters
        by adding a line to the parameter file. If multiple chains and doing
        path sampling, stores log-likelihoods of all chains.
        
        """
        # Gather log-likelihoods, and if path sampling save in path_sample list for later
        lnLikes = []
        doing_path_sampling = self.phycas.nchains > 1 and not self.phycas.is_standard_heating
        for i,c in enumerate(self.chains):
            lnLi = c.chain_manager.getLastLnLike()
            lnLikes.append(lnLi)
            if doing_path_sampling:
                self.phycas.path_sample[i].append(lnLi)
        
        # Only record samples from the current cold chain
        cold_chain = self.phycas.mcmc_manager.getColdChain()
        
        # Add line to parameter file if it exists
        if self.phycas.paramf:
            # old way: cycle, lnL_1, TL
            #lnL = cold_chain.chain_manager.getLastLnLike()
            #self.phycas.paramf.write('%d\t%.3f\t%.3f' % (cycle + 1, lnL, cold_chain.tree.edgeLenSum()))
            
            # new way: cycle, lnL_1, lnL_2, ..., lnL_nchains, TL
            self.phycas.paramf.write('%d\t' % (cycle + 1))
            for lnl in lnLikes:
                self.phycas.paramf.write('%.3f\t' % lnl)
            self.phycas.paramf.write('%.3f' % cold_chain.tree.edgeLenSum())
            
            self.phycas.paramf.write(cold_chain.model.paramReport())
            if self.phycas.using_hyperprior:
                for p in cold_chain.chain_manager.getEdgeLenHyperparams():
                    self.phycas.paramf.write('\t%.5f' % p.getCurrValue())
            if self.phycas.use_flex_model:
                rates_vector = cold_chain.likelihood.getRateMeans()
                for rr in rates_vector:
                    self.phycas.paramf.write('\t%.5f' % rr)
                probs_vector = cold_chain.likelihood.getRateProbs()
                for rp in probs_vector:
                    self.phycas.paramf.write('\t%.5f' % rp)
            self.phycas.paramf.write('\n')
        
        # Add line to tree file if it exists
        if self.phycas.treef:
            self.phycas.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, cold_chain.tree.makeNewick()))

