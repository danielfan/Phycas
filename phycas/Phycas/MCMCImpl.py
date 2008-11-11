import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from MCMCManager import MCMCManager
from phycas.ProbDist import StopWatch
from phycas.ReadNexus import NexusReader

def check(msg = 'check'):
    raw_input(msg)

class MCMCImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Needs to be written.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes MCMCImpl object by assigning supplied phycas object
        to a data member variable.
        
        """
        CommonFunctions.__init__(self, opts)
        # These copied over from Phycas.py - many are not used and should be weeded out
        self.data_matrix            = None
        self.file_name_trees_stored = None
        self.do_marginal_like       = False
        self.mcmc_manager           = MCMCManager(self)
        self.heat_vector            = None      # Leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using self.opts.nchains and self.heating_lambda
        self.stopwatch              = StopWatch()
        self.sim_model_tree         = None      # Will hold the model tree used by simulateDNA 
        self.starting_tree          = None      # Will contain description of actual starting tree used
        self.warn_tip_numbers       = False     # True only if tip numbers were not able to be created using the tip names in the tree description (always False if starting_tree_source == 'random' because BuildTreeFromString is not called in this case)
        self.ntax                   = 0         # Will hold the actual number of taxa after data file read
        self.nchar                  = 0         # Will hold the actual number of characters after data file has been read
        self.npatterns              = 0         # Will hold the actual number of patterns after data file has been read
        self.taxon_labels           = []        # Will hold taxon labels from data file or default names if self.data_source equals None
        self.paramf                 = None
        self.treef                  = None
        #self.tree_file_name         = ''        # Will hold tree file name (see openParameterAndTreeFiles)
        #self.param_file_name        = ''        # Will hold parameter file name (see openParameterAndTreeFiles)
        #self.tmp_simdata            = SimData()
        self.gg_Pm                  = 0.0       # Penalty component (same for all k)
        self.gg_Gm                  = []        # Vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []        # Vector of overall measures (one for each k in gg_kvect)
        self.reader                 = NexusReader()
        #self._logFileName           = None
        self.addition_sequence      = []        # List of taxon numbers for addition sequence
        self.samc_theta             = []        # Normalizing factors (will have length ntax - 3 because levels with 1, 2 or 3 taxa are not examined)
        self.samc_distance_matrix   = None      # Holds ntax x ntax hamming distance matrix used by SamcMove
        self.stored_tree_defs       = None
        self.psf                    = None
        self.pdf_splits_to_plot     = None
        #self.param_file_name        = None  
        #self.tree_file_name         = None
        self.nsamples               = None
        self.unimap_manager         = None
        self.nsamples               = 0
        self.burnin                 = 0     # same as self.opts.burnin except for path sampling, when it drops to 0 after first beta value
        self.last_adaptation        = 0
        self.next_adaptation        = 0
        self.ps_beta                = 1.0
        self.ps_beta_index          = 0
        self.ps_sampled_betas       = None
        self.ps_sampled_likes       = None
        #self.ps_delta_beta          = 0.0   # can be deleted when new system in place
        self.phycassert(self.opts.ps_nbetavals > 0, 'ps_nbetavals cannot be less than 1')

    def sliceSamplerReport(self, s, nm):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reports status of one slice sampler s that is updating a parameter
        named nm.
        
        """
        avg = float(s.getNumFuncEvals())/float(s.getNumSamples())
        mode = s.getMode()
        return '  mode=%.5f, avgevals=%.3f (%s)\n' % (mode, avg, nm)

    def adaptOneSliceSampler(self, p):
        self.phycassert(p, 'could not adapt slice sampler; parameter non-existant')
        summary = ''
        if p.hasSliceSampler():
            s = p.getSliceSampler()
            nm = p.getName()
            if s.getNumSamples() > 0:
                s.adaptSimple(self.opts.adapt_simple_param)
                #self.total_evals += float(s.getNumFuncEvals())
                summary = self.sliceSamplerReport(s, nm)
                s.resetDiagnostics()
        return summary
        
    def adaptSliceSamplers(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Cycles through all slice samplers and adapts each one. Adaptation of
        a slice sampler involves changing the unit width of segments used to
        construct a slice. If the slice unit width is too small, too many
        likelihood function evaluations are needed to span the full
        conditional posterior density. If the slice unit width is too large,
        too many failed sampling attempts are needed before a valid sample
        can be obtained. Adaptation adjusts the slice unit width of each
        slice sampler in an attempt to bring it closer to the optimum width
        using experience from past sampling attempts.
        
        """
        summary = ''
        # need to adapt all chains, not just the cold one!
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getAllUpdaters():
            summary += self.adaptOneSliceSampler(p)
        
        if self.opts.verbose and summary != '':
            self.output('\nSlice sampler diagnostics:')
            self.output(summary)

    def updateAllUpdaters(self, chain, chain_index, cycle):
        if self.opts.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
        for p in chain.chain_manager.getAllUpdaters():
            w = p.getWeight()
            #print p.getName(), "(weight = %d)" % w
            #if p.getName() == 'State freq move':
            #    raw_input('debug stop')
            for x in range(w):
                if self.opts.debugging:
                    p.setSaveDebugInfo(True)
                p.update()
                if self.opts.debugging:
                    tmpf.write('%s | %s\n' % (p.getName(), p.getDebugInfo()))
        
        if self.opts.debugging:
            tmpf.close()

    def showTopoPriorInfo(self):
        m = self.mcmc_manager.getColdChain()
        self.output('Topology prior:')
        if not self.opts.allow_polytomies:
            self.output('  flat across all fully-resolved tree topologies (polytomies not allowed)')
        else:            
            if m.topo_prior_calculator.isPolytomyPrior():
                self.output('  Prior type: polytomy prior')
            else:
                self.output('  Prior type: resolution class prior')
            self.output('  Prior strength (C): %s' % m.topo_prior_calculator.getC())
            self.output('  Prior probability for each resolution class:')
            self.output('  Note: 0.00000000 does *not* mean that the prior is zero! It simply')
            self.output('        indicates that the prior is less than 0.000000005\n')
            self.output('%20s %20s' % ('internal nodes', 'prior probability'))
            self.output('%20s %20s' % ('--------------', '-----------------'))
            topo_priors = m.topo_prior_calculator.getRealizedResClassPriorsVect()
            for i,v in enumerate(topo_priors):
                if i == 0:
                    denom = v   # first element of vector is log of normalization constant (sum of all other elements)
                else:
                    topo_prior = math.exp(v - denom)
                    self.output('%20d %20.8f' % (i,topo_prior))
            self.output()

    def showParamInfo(self, p):
        self.output('  Parameter name:     %s' % p.getName())
        self.output('  Prior distribution: %s' % p.getPriorDescr())
        if p.isMasterParameter():
            self.output('  Master parameter (no current value)')
        else:
            self.output('  Current value:      %s' % p.getCurrValue())
        self.output('  Prior log-density:  %s' % p.getLnPrior())
        self.output()
                
    def treeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the tree file and writes a translate table.
        
        """
        #self.tree_file_name = self.opts.out.trees
        
        tree_file_spec = self.opts.out.trees
        self.treef = None
        try:
            self.treef = tree_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open tree file (%s) failed.' % self.opts.out.trees.filename

        if self.treef:
            self.mcmc_manager.treeFileHeader(self.treef)

    def treeFileClose(self):
        self.treef.write('end;\n')
        self.treef.close()

    def paramFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the parameter file and writes a header line.
        
        """
        #self.param_file_name = self.opts.out.params
        
        param_file_spec = self.opts.out.params
        self.paramf = None
        try:
            self.paramf = param_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open parameter file (%s) failed.' % self.opts.out.params.filename

        if self.paramf:
            self.mcmc_manager.paramFileHeader(self.paramf)
            self.paramf.write('\n')

    def paramFileClose(self):
        self.paramf.close()

    def openParameterAndTreeFiles(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates parameter and tree file names based on the data file name or the
        user-supplied prefix and opens the files
        
        """
        #prefix = self.getPrefix()
        #self.param_file_name = prefix + '.p'
        #self.tree_file_name = prefix + '.t'

        self.paramFileOpen()
        self.treeFileOpen()
        
    def _loadData(self, matrix):
        self.data_matrix = matrix
        if matrix is None:            
            self.taxon_labels = []
            self.ntax = 0
            self.nchar = 0 # used for Gelfand-Ghosh simulations only
        else:
            self.taxon_labels = matrix.taxa
            self.ntax = self.data_matrix.getNTax()
            self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")
        
    def calcMarginalLikelihood(self):
        marginal_like = 0.0
        if self.opts.doing_path_sampling:
            if self.opts.ps_nbetavals > 1:
                # Calculate marginal likelihood using continuous path sampling
                for i in range(self.opts.ps_nbetavals):
                    n = len(self.ps_sampled_likes[i])
                    self.phycassert(n == self.opts.ncycles//self.opts.sample_every, 'number of sampled likelihoods (%d) does not match the expected number (%d) in path sampling calculation' % (n, self.opts.ncycles//self.opts.sample_every))
                    mean = sum(self.ps_sampled_likes[i])/float(n)
                    if i == 0 or i == self.opts.ps_nbetavals - 1:
                        marginal_like += (mean/2.0)
                    else:
                        marginal_like += mean
                marginal_like /= float(self.opts.ps_nbetavals - 1)
                self.output('Log of marginal likelihood (continuous path sampling method) = %f' % marginal_like)
            else:
                self.output('Log of marginal likelihood not computed (ps_nbetavals must be greater than 1)')
        else:
            # Calculate marginal likelihood using harmonic mean method on cold chain
            nignored = 0
            n = len(self.ps_sampled_likes[0])
            self.phycassert(n == self.opts.ncycles//self.opts.sample_every, 'number of sampled likelihoods (%d) does not match the expected number (%d) in harmonic mean calculation' % (n, self.opts.ncycles//self.opts.sample_every))
            min_lnL = min(self.ps_sampled_likes[0])
            sum_diffs = 0.0
            for lnl in self.ps_sampled_likes[0]:
                diff = lnl - min_lnL
                if diff < 500.0:
                    sum_diffs += math.exp(-diff)
                else:
                    nignored += 1
            log_harmonic_mean = math.log(n) + min_lnL - math.log(sum_diffs)
            if nignored > 0:
                self.warning('ignoring %d sampled log-likelihoods in harmonic mean calculation' % nignored)
            self.output('Log of marginal likelihood (harmonic mean method) = %f' % log_harmonic_mean)

    def getStartingTree(self):
        if self.starting_tree is None:
            if False:
                if self.opts.starting_tree_source == 'file':
                    self.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")
                    tree_defs = self.reader.getTrees()
                    self.phycassert(len(tree_defs) > 0, 'a trees block defining at least one tree must be stored in the nexus data file')
                    # Grab first tree description in the data file
                    # TODO allow some other tree than the first
                    self.starting_tree = tree_defs[0]
                elif self.opts.starting_tree_source == 'usertree':
                    self.starting_tree = Newick(self.opts.tree_topology)
                elif self.opts.starting_tree_source == 'random':
                    self.phycassert(self.ntax > 0, 'expecting ntax to be greater than 0')
                    self.starting_tree = None
                else:
                    self.phycassert(False, "starting_tree_source should equal 'random', 'file', or 'usertree', but instead it was this: %s" % self.starting_tree_source)
            else:
                # If user failed to specify starting_tree_source, get starting tree from randomtree object
                # as it is currently configured
                tr_source = self.opts.starting_tree_source
                if tr_source is None:
                    tr_source = randomtree()
                try:
                    tr_source.setActiveTaxonLabels(self.taxon_labels)
                    i = iter(tr_source)
                    self.starting_tree = i.next()
                except:
                    self.stdout.error("A starting tree could not be obtained from the starting_tree_source")
                    raise
        return self.starting_tree

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
        3) creates an appropriate heat_vector
        4) calls MCMCManager's createChains function to handle setup for each
        individual chain; 
        5) opens the parameter and tree files; and
        6) establishes an output log file name if requested
        
        """
        if self.opts.model.edgelen_prior is not None:
            # set both internal and external edge length priors to edgelen_prior
            self.opts.model.internal_edgelen_prior = self.opts.model.edgelen_prior
            self.opts.model.external_edgelen_prior = self.opts.model.edgelen_prior
        else:
            # Ensure that user has specified both internal and external edge length priors
            self.phycassert(self.opts.model.internal_edgelen_prior is not None, 'internal_edgelen_prior cannot be None if edgelen_prior is None')
            self.phycassert(self.opts.model.external_edgelen_prior is not None, 'external_edgelen_prior cannot be None if edgelen_prior is None')
            
        if self.opts.model.edgelen_hyperprior is not None:
            # Ensure that both internal and external edgelen priors are Exponential
            if self.opts.model.internal_edgelen_prior.getDistName() != 'Exponential':
                self.opts.model.internal_edgelen_prior = Exponential(1.0)
                self.warning('internal_edgelen_prior reset to Exponential because edgelen_hyperprior was specified')
            if self.opts.model.external_edgelen_prior.getDistName() != 'Exponential':
                self.opts.model.external_edgelen_prior = Exponential(1.0)
                self.warning('external_edgelen_prior reset to Exponential because edgelen_hyperprior was specified')

        ds = self.opts.data_source
        mat = ds and ds.getMatrix() or None
        self._loadData(mat)
        
        # Create a tree description to be used for building starting trees
        self.getStartingTree()

        # Determine heating levels if multiple chains
        if self.heat_vector == None:
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
            # User supplied his/her own heat_vector; perform sanity checks
            self.opts.nchains = len(self.heat_vector)
            self.phycassert(self.heat_vector.index(1.0) < self.opts.nchains, 'user-supplied heat_vector does not allow for a cold chain (one power must be 1.0)')

        self.mcmc_manager.createChains()
        self.openParameterAndTreeFiles()
        if self.opts.doing_path_sampling:
            self.ps_beta = self.opts.ps_maxbeta
            cc = self.mcmc_manager.getColdChain()
            cc.setPower(self.ps_beta)
        
    def beyondBurnin(self, cycle):
        c = cycle + 1
        return (c > self.burnin)
        
    def doThisCycle(self, cycle, mod):
        c = cycle + 1
        return ((c % mod) == 0)
        
    def explorePrior(self, cycle):
        chain_index = self.mcmc_manager.getColdChainIndex()
        chain = self.mcmc_manager.getColdChain()
        tm = Phylogeny.TreeManip(chain.tree)
        if self.opts.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
        edgelens_generated = False
        for p in chain.chain_manager.getAllUpdaters():
            w = p.getWeight()
            name = p.getName()
            if name == 'edge length hyperparameter':    # C++ class HyperPriorParam
                # Choose hyperparam, then use it to choose new edge lengths for a newly-created tree
                if chain.model.isSeparateInternalExternalEdgeLenPriors():
                    edgelen_hyperparam = chain.model.getEdgeLenHyperPrior().sample()
                    chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                    chain.model.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                    edgelen_hyperparam = chain.model.getEdgeLenHyperPrior().sample()
                    chain.chain_manager.setEdgeLenHyperparam(1, edgelen_hyperparam)
                    chain.model.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                else:
                    edgelen_hyperparam = chain.model.getEdgeLenHyperPrior().sample()
                    chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                    chain.model.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                    chain.model.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                tm.equiprobTree(chain.tree.getNTips(), chain.r, chain.model.getInternalEdgeLenPrior(), chain.model.getExternalEdgeLenPrior())
                edgelens_generated = True
            elif name == 'trs/trv rate ratio':              # C++ class KappaParam
                new_kappa = chain.model.getKappaPrior().sample()
                chain.model.setKappa(new_kappa)
            elif name == 'nonsynon./synon. rate ratio':      # C++ class OmegaParam
                new_omega = chain.model.getOmegaPrior().sample()
                chain.model.setOmega(new_omega)
            elif name == 'rAC':                    # C++ class GTRRateParam
                new_rAC = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(0, new_rAC) 
            elif name == 'rAG':                    # C++ class GTRRateParam
                new_rAG = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(1, new_rAG) 
            elif name == 'rAT':                    # C++ class GTRRateParam
                new_rAT = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(2, new_rAT) 
            elif name == 'rCG':                    # C++ class GTRRateParam
                new_rCG = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(3, new_rCG) 
            elif name == 'rCT':                    # C++ class GTRRateParam
                new_rCT = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(4, new_rCT) 
            elif name == 'rGT':                    # C++ class GTRRateParam
                new_rGT = chain.model.getRelRatePrior().sample()
                chain.model.setRelRateUnnorm(5, new_rGT) 
            elif name == 'base freq. A':                    # C++ class StateFreqParam
                new_freq_param_A = chain.model.getStateFreqParamPrior().sample()
                chain.model.setStateFreqParam(0, new_freq_param_A) 
            elif name == 'base freq. C':                    # C++ class StateFreqParam
                new_freq_param_C = chain.model.getStateFreqParamPrior().sample()
                chain.model.setStateFreqParam(1, new_freq_param_C) 
            elif name == 'base freq. G':                    # C++ class StateFreqParam
                new_freq_param_G = chain.model.getStateFreqParamPrior().sample()
                chain.model.setStateFreqParam(2, new_freq_param_G) 
            elif name == 'base freq. T':                    # C++ class StateFreqParam
                new_freq_param_T = chain.model.getStateFreqParamPrior().sample()
                chain.model.setStateFreqParam(3, new_freq_param_T) 
            elif name == 'Discrete gamma shape':            # C++ class DiscreteGammaShapeParam
                new_gamma_shape = chain.model.getDiscreteGammaShapePrior().sample()
                chain.model.setShape(new_gamma_shape)
            elif name == 'Proportion of invariable sites':  # C++ class PinvarParam
                new_pinvar = chain.model.getPinvarPrior().sample()
                chain.model.setPinvar(new_pinvar)
            elif name == 'master edge length parameter':
                pass
            elif name == 'State freq move':                # C++ class StateFreqMove
                freq_vector = chain.model.getStateFreqPrior().sample()
                # should use generic function below (i.e. not specific to nucleotide data)
                chain.model.setNucleotideFreqs(freq_vector[0],freq_vector[1],freq_vector[2],freq_vector[3])
            elif name == 'Edge length move':                # C++ class EdgeMove
                pass
            elif name == 'Larget-Simon move':               # C++ class LargetSimonMove
                pass
            elif name == 'Tree scaler move':                # C++ class TreeScalerMove
                pass
            elif name == 'Bush move':                       # C++ class BushMove
                pass    # polytomies handled further down (by randomly pruning fully-resolved equiprobable tree)
            elif name == 'FLEX probs':                      # C++ class FlexProbParam
                self.phycassert(0, 'sampling directly from the prior not yet implemented for FLEX model (workaround: specify mcmc.draw_directly_from_prior = False)')
            elif name == 'FLEX rates':                      # C++ class FlexRateParam
                self.phycassert(0, 'sampling directly from the prior not yet implemented for FLEX model (workaround: specify mcmc.draw_directly_from_prior = False)')
            elif name == 'NCat move':                       # C++ class NCatMove
                self.phycassert(0, 'sampling directly from the prior not yet implemented for FLEX model (workaround: specify mcmc.draw_directly_from_prior = False)')
            elif name == 'Univent mapping move':            # C++ class MappingMove
                self.phycassert(0, 'sampling directly from the prior not yet implemented for unimap analyses (workaround: specify mcmc.draw_directly_from_prior = False)')
            elif name == 'Unimap NNI move':                 # C++ class UnimapNNIMove
                self.phycassert(0, 'sampling directly from the prior not yet implemented for unimap analyses (workaround: specify mcmc.draw_directly_from_prior = False)')
            elif name == 'SAMC move':                       # C++ class SamcMove
                self.phycassert(0, 'sampling directly from the prior not yet implemented for SAMC (workaround: specify mcmc.draw_directly_from_prior = False)')
            else:
                self.phycassert(0, 'model uses an updater (%s) that has not yet been added to MCMCImpl.explorePrior (workaround: specify mcmc.draw_directly_from_prior = False)' % name)

        # If no edge length hyperprior was specified, then build the tree with edge lengths now
        if not edgelens_generated:
            tm.equiprobTree(chain.tree.getNTips(), chain.r, chain.model.getInternalEdgeLenPrior(), chain.model.getExternalEdgeLenPrior())
            
        if self.opts.allow_polytomies:
            # Choose number of internal nodes
            num_internal_nodes = chain.topo_prior_calculator.sample(chain.r)
                    
            # Delete nodes at random from tree to achieve chosen number of internal nodes
            orig_num_internal_nodes = chain.tree.getNInternals()
            num_internals_to_delete = orig_num_internal_nodes - num_internal_nodes
            for i in range(num_internals_to_delete):
                tm.deleteRandomInternalEdge(chain.r)
                
        chain.prepareForLikelihood()
        chain.likelihood.replaceModel(chain.model)
        
        if False:
            # debugging code
            chain.likelihood.storeSiteLikelihoods(True)
            from phycas.Utilities.kappa2tratio import convert
            f = chain.model.getStateFreqs()
            k = convert(chain.model.getKappa(), f[0], f[1], f[2], f[3])
            print 'cycle = %d, model = %s' % (cycle + 1, chain.model.getModelName())
            print '  lset tratio=%.5f basefreq=(%.5f %.5f %.5f) rates=gamma ncat=4 shape=%.5f;' % (k, f[0], f[1], f[2], chain.model.getShape())
            print 'taxon names:', self.opts.data_source.taxon_labels
            chain.tree.rectifyNames(self.opts.data_source.taxon_labels)
            print '  utree one = %s;' % chain.tree.makeNewick()
            print '  sum of edge lengths = %.5f' % chain.tree.edgeLenSum()
            raw_input('stopped before computing likelihood')
        
        # recalculate the likelihood
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        cold_chain_manager.refreshLastLnLike()
        
        if False:
            # debugging code
            counts = chain.likelihood.getPatternCounts()
            sitelikes = chain.likelihood.getSiteLikelihoods()
            print '  lnL = %.6f' % cold_chain_manager.getLastLnLike()
            sumlikes = 0.0
            for sitelike,count in zip(sitelikes, counts):
                if count > 100:
                    print '%6.0f -> %12.5f' % (count, sitelike)
                sumlikes += count*sitelike
            raw_input('check: sum = %.5f' % sumlikes)
        
        if self.opts.debugging:
            tmpf.close()
        
    def mainMCMCLoop(self, explore_prior = False):
        self.last_adaptation = 0
        self.next_adaptation = self.opts.adapt_first
        for cycle in xrange(self.burnin + self.opts.ncycles):
            if explore_prior and self.opts.draw_directly_from_prior:
                self.explorePrior(cycle)
            else:
                for i,c in enumerate(self.mcmc_manager.chains):
                    self.updateAllUpdaters(c, i, cycle)
            if self.opts.verbose and self.doThisCycle(cycle, self.opts.report_every):
                self.stopwatch.normalize()
                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                if self.opts.doing_path_sampling:
                    msg = 'beta = %.5f, cycle = %d, lnL = %.5f (%.5f secs)' % (self.ps_beta, cycle + 1, cold_chain_manager.getLastLnLike(), self.stopwatch.elapsedSeconds())
                else:
                    msg = 'cycle = %d, lnL = %.5f (%.5f secs)' % (cycle + 1, cold_chain_manager.getLastLnLike(), self.stopwatch.elapsedSeconds())
                self.output(msg)
            if self.beyondBurnin(cycle) and self.doThisCycle(cycle - self.burnin, self.opts.sample_every):
                self.mcmc_manager.recordSample(cycle)
                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                sampled_lnL = cold_chain_manager.getLastLnLike()
                self.ps_sampled_likes[self.ps_beta_index].append(sampled_lnL)
                self.stopwatch.normalize()
            if self.doThisCycle(cycle, self.next_adaptation):
                self.adaptSliceSamplers()
                self.next_adaptation += 2*(self.next_adaptation - self.last_adaptation)
                self.last_adaptation = cycle + 1
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the MCMC analysis. 
        
        """
        self.setup()
        
        # If user has set quiet to True, then phycas.output calls will have no effect
        self.quiet = self.opts.quiet
        
        # Tell TreeLikelihood object if user wants to run with no data
        #if not self.data_source:
        #    self.likelihood.setNoData()

        self.nsamples = self.opts.ncycles//self.opts.sample_every
        
        if self.opts.verbose:
            if self.data_matrix == None:
                self.output('Data source:    None (running MCMC with no data to explore prior)')
            else:
                self.output('Data source:    %s' % str_value_for_user(self.opts.data_source))
                all_missing = self.mcmc_manager.getColdChain().likelihood.getListOfAllMissingSites()
                num_excluded = len(all_missing)
                if num_excluded > 0:
                    self.output('*** Note: the following %d sites were automatically excluded because' % num_excluded)
                    self.output('*** they exhibited completely missing data for all taxa:')
                    while len(all_missing) > 0:
                        tmp = all_missing[:10]
                        all_missing = all_missing[10:]
                        self.output('***   '+','.join([str(i+1) for i in tmp]))
            self.output('Starting tree:  %s' % self.starting_tree)
            if self.opts.doing_path_sampling:
                self.output('\nPerforming path sampling MCMC to estimate marginal likelihood.')
                self.output('Likelihood will be raised to the power beta, and beta will be')
                self.output('decremented from 1.0 to 0.0 in a series of steps.')
                self.output('  No. steps:               %s' % self.opts.ps_nbetavals)
                self.output('  No. cycles per step:     %s' % self.opts.ncycles)
                self.output('  Sample every:            %s' % self.opts.sample_every)
                self.output('  No. samples per step:    %s' % self.nsamples)
                self.output('\n')
            else:
                self.output('No. cycles:     %s' % self.opts.ncycles)
                self.output('Sample every:   %s' % self.opts.sample_every)
                self.output('No. samples:    %s' % self.nsamples)
            self.output('Sampled trees will be saved in %s' % str_value_for_user(self.opts.out.trees))
            self.output('Sampled parameters will be saved in %s' % str_value_for_user(self.opts.out.params))
            if self.opts.use_unimap:
                self.output('Using uniformized mapping MCMC')
            else:
                self.output('Using standard MCMC (i.e. no uniformized mapping)')

            if not self.warn_tip_numbers:
                self.output('Tip node numbers were set using the names in the tree description')
            else:
                self.output('Warning: tip node numbers were NOT set using the names in the tree description')

        if self.opts.nchains == 1:
            self.output('Creating one chain (i.e. not using heated chains to improve mixing)')
        else:
            self.output('Creating %d chains with these temperatures:' % (self.opts.nchains))
            for t in self.heat_vector:
                self.output('  %.5f %s' % (t, t == 1.0 and '(cold chain)' or ''))
            
        # Compute the current log-likelihood and log-prior in case first updater 
        # is a move and will thus depend on these quantities being accurate
        for c in self.mcmc_manager.chains:
            c.chain_manager.refreshLastLnLike()
            c.chain_manager.refreshLastLnPrior()
            if c.heating_power == 1.0:
                self.output('Starting log-likelihood = %s' % c.chain_manager.getLastLnLike())
                self.output('Starting log-prior = %s' % c.chain_manager.getLastLnPrior())

        # Show starting parameter info 
        self.output('\nParameter starting values and prior densities:')
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getEdgeLenParams():
            self.showParamInfo(p)
        for p in cold_chain_manager.getEdgeLenHyperparams():
            self.showParamInfo(p)
        for p in cold_chain_manager.getModelParams():
            self.showParamInfo(p)

        # Debugging: show data patterns
        if self.opts.debugging:
            cold_chain = self.mcmc_manager.getColdChain()
            s = cold_chain.likelihood.listPatterns()
            print '\nDebug Info: List of data patterns and their frequencies:'
            print s

        # Show information about topology prior to be used
        self.showTopoPriorInfo()

        self.stopwatch.start()
        self.mcmc_manager.resetNEvals()
        
        if self.opts.doing_path_sampling:
            self.output('\nSampling (%d cycles for each of the %d values of beta)...' % (self.opts.ncycles, self.opts.ps_nbetavals))
        else:
            self.output('\nSampling (%d cycles)...' % self.opts.ncycles)
        if self.opts.verbose:
            print
            
        # POL moved these lines to beginning of mainMCMCLoop so that adaptation cycle
        # starts again each time path sampling beta value is changed
        #self.last_adaptation = 0
        #self.next_adaptation = self.opts.adapt_first

        # Lay down first line in params file (recorded as cycle 0) containing starting values of parameters
        self.mcmc_manager.recordSample()

        if self.opts.doing_path_sampling:
            self.phycassert(self.data_matrix is not None, 'path sampling requires data')
            self.phycassert(self.opts.nchains == 1, 'path sampling requires nchains to be 1')
            chain = self.mcmc_manager.getColdChain()
            if self.opts.ps_nbetavals > 1:
                if False:
                    # old way (constant ps_delta_beta)
                    self.ps_delta_beta = float(self.opts.ps_maxbeta - self.opts.ps_minbeta)/float(self.opts.ps_nbetavals - 1)
                    self.ps_sampled_betas = [self.opts.ps_maxbeta - self.ps_delta_beta*float(i) for i in range(self.opts.ps_nbetavals)]
                else:
                    # new way (ps_delta_beta is a vector of values determined by discretized beta distribution)
                    # Beta distribution will be divided into ps_nbetavals intervals, each of which has an equal area
                    segment_area = 1.0/float(self.opts.ps_nbetavals - 1)
                    cum_area = 0.0
                    lower_boundary = 0.0
                    self.ps_sampled_betas = [self.opts.ps_minbeta]
                    total_extent = float(self.opts.ps_maxbeta - self.opts.ps_minbeta)
                    betadist = ProbDist.Beta(self.opts.ps_shape1, self.opts.ps_shape2)
                    for i in range(self.opts.ps_nbetavals - 1):
                        cum_area += segment_area
                        upper_boundary = betadist.getQuantile(cum_area)
                        scaled_upper_boundary = self.opts.ps_minbeta + total_extent*upper_boundary
                        self.ps_sampled_betas.append(scaled_upper_boundary)
                        lower_boundary = upper_boundary
                        
                    # Reverse the sampled betas so that they start at 1.0 and decrease toward 0.0
                    self.ps_sampled_betas.reverse()
                
                # Output the beta values that will be used
                self.output('%d %s chosen from a discrete\nBeta(%.5f, %.5f) distribution:' % (self.opts.ps_nbetavals, (self.opts.ps_nbetavals == 1 and 'value was' or 'values were'), self.opts.ps_shape1, self.opts.ps_shape2))
                for i,x in enumerate(self.ps_sampled_betas):
                    self.output('%6d %12.5f' % (i+1,x))
                self.output('An MCMC analysis will be performed exploring each of the')
                self.output('power posteriors defined by these values.')
                self.output()
            else:
                self.ps_sampled_betas = [self.opts.ps_minbeta]
            
            self.ps_sampled_likes = []
            for self.ps_beta_index,self.ps_beta in enumerate(self.ps_sampled_betas):
                self.ps_sampled_likes.append([])
                chain.setPower(self.ps_beta)
                if self.ps_beta_index > 0:
                    self.burnin = 0
                else:
                    self.burnin = self.opts.burnin
                if self.ps_beta == 0.0:
                    self.mainMCMCLoop(explore_prior=True)
                else:
                    self.mainMCMCLoop()
        else:
            self.ps_sampled_likes = []
            self.ps_sampled_likes.append([])
            self.ps_beta_index = 0
            self.burnin = self.opts.burnin
            if self.data_matrix is None:
                self.mainMCMCLoop(explore_prior=True)
            else:
                self.mainMCMCLoop()

        self.adaptSliceSamplers()
        total_evals = self.mcmc_manager.getTotalEvals() #self.likelihood.getNEvals()
        total_secs = self.stopwatch.elapsedSeconds()
        self.output('%d likelihood evaluations in %.5f seconds' % (total_evals, total_secs))
        if (total_secs > 0.0):
            self.output('  = %.5f likelihood evaluations/sec' % (total_evals/total_secs))

        if self.treef:
            self.treeFileClose()
        if self.paramf:
            self.paramFileClose()
            
        # Report marginal likelihood calculation (using harmonic mean method, or path sampling if 
        # this mcmc method was called from the ps command)
        self.calcMarginalLikelihood()
