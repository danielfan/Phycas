import os, sys, math, threading
import MCMCManager  # poorly named, as MCMCManager is now only one of many classes within
from phycas.Conversions import *
from phycas.DataMatrix import *
from phycas.Likelihood import *
from phycas.Phylogeny import *
from phycas.ProbDist import *
from phycas.ReadNexus import *

# see http://mail.python.org/pipermail/python-list/2002-January/121376.html
import inspect

class Phycas(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Performs a Bayesian phylogenetic MCMC analysis. The tree topology and
    edge lengths are updated via the Metropolis-Hastings algorithm (using
    the Larget-Simon LOCAL move without a molecular clock). Slice 
    sampling is used to update all model parameters except edge lengths.

    For examples of how to use Phycas, see the Examples folder:
    phycas/Phycas/Phycas.py   <-- you are here
    phycas/Examples           <-- here is the Examples directory

    See the __init__ function below for variables that can be modified
    before Phycas is run.

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the object with default values for all settings. If these
        defaults do not suit, they can be changed before mcmc() is called.
        
        """
        # Variables controlling the MCMC analysis and progress reporting
        self.random_seed            = 0         # determines random number seed (0 means generate seed automatically)
        self.ncycles                = 10000     # the number of update cycles (a cycle is analogous to, but different than, a "generation" in MrBayes)
        self.sample_every           = 100       # a new sample of tree topology and model parameters will be taken after this many cycles
        self.report_every           = 100       # a progress report will be displayed after this many cycles
        self.verbose                = True      # you will get more output if True, less output if False
        self.quiet                  = False     # if True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well
        self.outfile_prefix         = None      # If None, parameter and tree files created will have a name beginning with the name of the data file; if provided, this prefix will form the first part of the parameter and tree file names

        # Variables associated with substitution models (except for edge lengths)
        self.default_model          = 'hky'     # Can be 'jc', 'hky' or 'gtr'

        self.relrate_prior          = ExponentialDist(1.0)
        self.starting_relrates      = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
        self.fix_relrates           = False 

        self.kappa_prior            = ExponentialDist(1.0) 
        self.starting_kappa         = 4.0
        self.fix_kappa              = False 

        self.num_rates              = 1         # default is no rate heterogeneity (i.e. 1 rate)
        self.gamma_shape_prior      = ExponentialDist(1.0)
        self.starting_shape         = 0.5
        self.fix_shape              = False 
        self.use_inverse_shape      = False      # if True, gamma_shape_prior is applied to 1/shape rather than shape

        self.estimate_pinvar        = False
        self.pinvar_prior           = BetaDist(1.0, 1.0)
        self.starting_pinvar        = 0.2
        self.fix_pinvar             = False 

        self.base_freq_param_prior  = ExponentialDist(1.0)
        self.starting_freqs         = [1.0, 1.0, 1.0, 1.0] 
        self.fix_freqs              = False 

        # Variables associated with the source of data        
        self.data_source            = 'file'    # Specify None to explore the joint prior or to simulate data; if 'file', self.data_file_name should be a valid nexus file name
        self.data_file_name         = ''        # will hold actual data file name

        # Variables associated with simulating data (see function simulateDNA below)
        self.sim_file_name          = 'simulated.nex'   # name of file in which to save simulated data
        self.sim_taxon_labels       = ['taxon1', 'taxon2', 'taxon3', 'taxon4']  # names to use for taxa in simulated data set (number of labels defined determines the number of taxa in the simulated dataset)
        self.sim_nchar              = 1000      # number of characters to generate

        # Variables associated with the source of starting tree
        self.starting_tree_source   = 'random'  # source of starting tree topology: can be either 'random' or 'usertree'
        self.tree_topology          = None      # unused unless starting_tree_source is 'usertree'

        # Variables associated with Larget-Simon moves
        self.ls_move_lambda         = 0.2       # The value of the tuning parameter for the Larget-Simon move
        self.ls_move_weight         = 100       # Larget-Simon moves will be performed this many times per cycle
        self.ls_move_debug          = False     # If set to true, TreeViewer will popup on each Larget-Simon move update showing edges affected by the proposed move
        
        # Variables associated with tree scaler move
        self.tree_scaler_weight     = 0         # whole-tree scaling will be performed this many times per cycle

        # Variables associated with Polytomy (Bush) moves
        self.allow_polytomies       = False     # if True, do Bush moves in addition to Larget-Simon moves; if False, do Larget-Simon moves only
        self.polytomy_prior         = True      # if True, use polytomy prior; if False, use resolution class prior
        self.topo_prior_C           = 2.0       # specifies the strength of the prior (C = 1 is flat prior; C > 1 favors less resolved topologies)
        self.bush_move_edgelen_mean = 1.0       # specifies mean of exponential edge length generation distribution used by BushMove when new edges are created
        self.bush_move_weight       = 100       # Bush moves will be performed this many times per cycle if
        self.bush_move_debug        = False     # If set to true, TreeViewer will pop up on each Bush move update showing edges affected by the proposed move

        # Variables associated with SAMC analyses
        self.doing_samc             = False     # if True, using Cheon and Liang "SSAMC" method
        self.samc_move_edgelen_mean = 1.0       # specifies mean of exponential edge length generation distribution used by SamcMove when new edges are created
        self.samc_move_debug        = False     # If set to True, output will be saved to a file named cf.txt (if it doesn't exist) or, if cf.txt already exists, the output will be compared to cf.txt and the program will halt if a discrepency is found
        self.samc_t0                = 10000.0   # samc_gain_factor = samc_t0/max(samc_t0, cycle)
        self.samc_move_weight       = 1         # number of times per cycle that SAMC moves will be performed (currently unused because SAMC moves are not used in standard MCMC analyses)
        self.samc_temperature       = 0.6       # temperature used in extrapolate move to smooth out differences in probabilities of different possible attachment points

        # Variables associated with slice samplers
        self.slice_weight           = 1         # Slice sampled parameters will be updated this many times per cycle
        self.slice_max_units        = 1000      # Max. number of units used in slice sampling
        self.adapt_first            = 100       # Adaptation of slice samplers is performed the first time at cycle adapt_first.
                                                # Subsequent adaptations wait twice the number of cycles as the previous adaptation.
                                                # Thus, adaptation n occurs at cycle adapt_first*(2^n - 1)
                                                # The total number of adaptations that will occur during an MCMC run is
                                                #      [ln(adapt_first + ncycles) - ln(adapt_first)]/ln(2)
        self.adapt_simple_param     = 0.5       # Slice sampler adaptation parameter
        #self.adapt_ycond_param      = 1.3
        #self.adapt_ycond_from_ends  = 0.25

        # Variables associated with edge length prior distributions
        # If using_hyperprior is False, external_edgelen_dist and internal_edgelen_dist will govern
        # external and internal edge length priors, respectively. If either external_edgelen_dist or
        # internal_edgelen_dist is None, then whichever prior is defined will govern both internal 
        # and external edges. If using_hyperprior is True, then hyperparameters will be used to
        # govern the mean of all edge length priors defined. The same hyperprior will be used for both
        # hyperparameters (should two hyperparameters be used)
        self.using_hyperprior       = True      
        self.edgelen_hyperprior     = InverseGammaDist(2.1, 1.0/1.1)
        self.starting_edgelen_hyperparam = 0.05 #POL doesn't do anything! currently ignored
        self.fix_edgelen_hyperparam = False 

        # Note that edgelen_dist is a property: its setter function (setEdgelenDist) sets
        # both internal_edgelen_dist and external_edgelen_dist to the value specified
        self.edgelen_dist           = ExponentialDist(2.0)
        self.fix_edgelens           = False 
        
        # Variables associated with initializing the MCMC sampler
        self.starting_edgelen_dist  = ExponentialDist(10.0) # Used to select the starting edge lengths when starting_tree_source is 'random'

        # Variables associated with Metropolis coupling (heated chains)
        self.heating_lambda         = 0.2
        self.nchains                = 4
        self.is_standard_heating    = True
        
        # Variables associated with Gelfand-Ghosh calculation
        self.gg_outfile             = 'gg.txt'  # file in which to save gg results (use None to not save results)
        self.gg_nreps               = 1         # the number of replicate simulations to do every MCMC sample
        self.gg_kvect               = [1.0]     # vector of k values to use when computing Gm and Dm
        self.gg_save_postpreds      = False     # if True, all posterior predictive data sets will be saved
        self.gg_postpred_prefix     = 'pp'      # prefix to use for posterior predictive dataset filenames (only used if gg_save_postpreds is True)
        self.gg_burnin              = 1         # number of starting samples to skip when computing Gelfand-Ghosh measures
        self.gg_pfile               = None      # name of parameter file to use for Gelfand-Ghosh calculations
        self.gg_tfile               = None      # name of tree file to use for Gelfand-Ghosh calculations
        self.gg_bin_patterns        = False     # if True, patterns will be classified into 7 bins, corresponding to 'A only', 'C only', 'G only', 'T only', 'any 2 states', 'any 3 states' and 'any 4 states'. Gelfand-Ghosh statistics will be computed on this vector of counts instead of the complete vector of pattern counts. Can only be used for DNA/RNA data.

        # Variables associated with the FLEXCAT model
        self.use_flex_model         = False
        self.flex_ncat_move_weight  = 1         # number of times each cycle to attempt an ncat move
        self.flex_num_spacers       = 1         # number of fake rates between each adjacent pair of real rates
        self.flex_phi               = 0.25      # proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)
        self.flex_L                 = 1.0       # upper bound of interval used for unnormalized relative rate parameter values
        self.flex_lambda            = 1.0       # parameter of Poisson prior on the number of extra categories
        self.flex_prob_param_prior  = ExponentialDist(1.0)

        # Variables associated with underflow protection                
        self.uf_num_edges           = 50        # number of edges to traverse before taking action to prevent underflow
        
        # ***** IT IS BEST NOT TO CHANGE ANYTHING BELOW HERE *****
        self.debugging              = False      # if set to True expect lots of debug output (e.g. data pattern table)
        self.data_matrix            = None
        self.file_name_data_stored  = None
        self.do_marginal_like       = False
        self.mcmc_manager           = MCMCManager.MCMCManager(self)
        self.heat_vector            = None      # leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using self.nchains and self.heating_lambda
        self.stopwatch              = StopWatch()
        self.sim_model_tree         = None      # will hold the model tree used by simulateDNA 
        self.starting_tree          = None      # will contain description of actual starting tree used
        self.warn_tip_numbers       = False     # True only if tip numbers were not able to be created using the tip names in the tree description (always False if starting_tree_source == 'random' because BuildTreeFromString is not called in this case)
        self.ntax                   = 0         # will hold the actual number of taxa after data file read
        self.nchar                  = 0         # will hold the actual number of characters after data file has been read
        self.npatterns              = 0         # will hold the actual number of patterns after data file has been read
        self.taxon_labels           = []        # will hold taxon labels from data file or default names if self.data_source equals None
        self.paramf                 = None
        self.treef                  = None
        self.tmp_simdata            = SimData()
        self.gg_Pm                  = 0.0       # penalty component (same for all k)
        self.gg_Gm                  = []        # vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []        # vector of overall measures (one for each k in gg_kvect)
        self.reader                 = NexusReader()
        self.logf                   = None
        self._logFileName           = None
        #self.use_tree_viewer        = False    # popup graphical TreeViewer to show trees during run POLPY_NEWWAY
        self.addition_sequence      = []        # list of taxon numbers for addition sequence
        self.samc_theta             = []        # normalizing factors (will have length ntax - 3 because levels with 1, 2 or 3 taxa are not examined)
        self.samc_distance_matrix   = None      # holds ntax x ntax hamming distance matrix used by SamcMove
        self.path_sample            = None
        
    # see http://mail.python.org/pipermail/python-list/2002-January/121376.html
    def source_line():
        return inspect.getouterframes(inspect.currentframe())[1][2]
    
    def phycassert(self, assumption, msg):
        if not assumption:
            if Phycas.PhycassertRaisesException:
                raise AssertionError(msg)
            sys.exit('Error: ' + msg)

    def setEdgelenDist(self, dist):
        self.internal_edgelen_dist = self.external_edgelen_dist = dist
        
    def getEdgelenDist(self):
        self.phycassert(self.internal_edgelen_dist is self.external_edgelen_dist, "There are separate distributions for internal and external edge lengths")
        return self.internal_edgelen_dist

    edgelen_dist = property(getEdgelenDist, setEdgelenDist)
    
    def shutdown(self):
        if self.logf:
            self.logf.close()

    def output(self, msg = None):
        if not self.quiet:
            if msg:
                print msg
            else:
                print
        if self.logf:
            if not msg:
                msg = ''
            self.logf.write('%s\n' % (msg))
            self.logf.flush()

    def showParamInfo(self, p):
        self.output('  Parameter name:     %s' % p.getName())
        self.output('  Prior distribution: %s' % p.getPriorDescr())
        if p.isMasterParameter():
            self.output('  Master parameter (no current value)')
        else:
            self.output('  Current value:      %s' % p.getCurrValue())
        self.output('  Prior log-density:  %s' % p.getLnPrior())
        self.output()
                
    def showTopoPriorInfo(self):
        m = self.mcmc_manager.getColdChain()
        self.output('Topology prior:')
        if not self.allow_polytomies:
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

    def treeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the tree file and writes a translate table.
        
        """
        self.treef = file(self.tree_file_name, 'w')
        self.mcmc_manager.treeFileHeader(self.treef)

    def treeFileClose(self):
        self.treef.write('end;\n')
        self.treef.close()

    def paramFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the parameter file and writes a header line.
        
        """
        self.paramf = file(self.param_file_name, 'w')
        self.mcmc_manager.paramFileHeader(self.paramf)
        self.paramf.write('\n')

    def paramFileClose(self):
        self.paramf.close()

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
                s.adaptSimple(self.adapt_simple_param)
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
        
        if self.verbose and summary != '':
            self.output('\nSlice sampler diagnostics:')
            self.output(summary)

    def updateAllUpdaters(self, chain, chain_index, cycle):
        if self.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
        for p in chain.chain_manager.getAllUpdaters():
            w = p.getWeight()
            for x in range(w):
                if self.debugging:
                    p.setSaveDebugInfo(True)
                p.update()
                if self.debugging:
                    tmpf.write('%s | %s\n' % (p.getName(), p.getDebugInfo()))
        
        if self.debugging:
            tmpf.close()

    def runMCMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the MCMC analysis. Assumes setupMCMC() has already been
        called.
        
        """
        # Tell TreeLikelihood object if user wants to run with no data
        #if not self.data_source:
        #    self.likelihood.setNoData()

        self.nsamples = self.ncycles//self.sample_every
        
        if self.verbose:
            if self.data_source == None:
                self.output('Data source:    None (running MCMC with no data to explore prior)')
            elif self.data_source == 'file':
                self.output('Data source:    %s' % self.data_file_name)
            else:
                self.output("Data source:    Unknown (something other than 'file' or None was specified for data_source)")
            self.output('No. cycles:     %s' % self.ncycles)
            self.output('Sample every:   %s' % self.sample_every)
            self.output('Starting tree:  %s' % self.starting_tree)
            self.output('No. samples:    %s' % self.nsamples)
            self.output('Sampled trees will be saved in %s' % self.tree_file_name)
            self.output('Sampled parameters will be saved in %s' % self.param_file_name)

            if not self.warn_tip_numbers:
                self.output('Tip node numbers were set using the names in the tree description')
            else:
                self.output('Warning: tip node numbers were NOT set using the names in the tree description')

        self.output('Creating %d chains with these temperatures:' % (self.nchains))
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
        self.output('Parameter starting values and prior densities:')
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getEdgeLenParams():
            self.showParamInfo(p)
        for p in cold_chain_manager.getEdgeLenHyperparams():
            self.showParamInfo(p)
        for p in cold_chain_manager.getModelParams():
            self.showParamInfo(p)

        # Debugging: show data patterns
        if self.debugging:
            cold_chain = self.mcmc_manager.getColdChain()
            s = cold_chain.likelihood.listPatterns()
            print '\nDebug Info: List of data patterns and their frequencies:'
            print s

        # Show information about topology prior to be used
        self.showTopoPriorInfo()

        self.stopwatch.start()
        self.mcmc_manager.resetNEvals()

        self.output('\nSampling (%d cycles)...' % self.ncycles)
        if self.verbose:
            print
        self.mcmc_manager.recordSample(0)
        last_adaptation = 0
        next_adaptation = self.adapt_first

        #if self.use_tree_viewer:
        #    self.tree_mutex = threading.Lock()
        #    self.tree_viewer = TreeViewer.TreeViewer(tree=self.tree, mutex=self.tree_mutex) # POLPY_NEWWAY
        #    self.tree_viewer.start() # POLPY_NEWWAY    
    
        for cycle in xrange(self.ncycles):
            for i,c in enumerate(self.mcmc_manager.chains):
                self.updateAllUpdaters(c, i, cycle)
            if self.verbose and (cycle + 1) % self.report_every == 0:
                self.stopwatch.normalize()
                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                msg = 'cycle = %d, lnL = %.5f (%.5f secs)' % (cycle + 1, cold_chain_manager.getLastLnLike(), self.stopwatch.elapsedSeconds())
                self.output(msg)
            if (cycle + 1) % self.sample_every == 0:
                self.mcmc_manager.recordSample(cycle)
                self.stopwatch.normalize()
            if (cycle + 1) % next_adaptation == 0:
                self.adaptSliceSamplers()
                next_adaptation += 2*(next_adaptation - last_adaptation)
                last_adaptation = cycle + 1

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

        # If we have been path sampling, compute marginal likelihood using lnL values
        # for each chain stored in self.path_sample
        if self.nchains > 1 and not self.is_standard_heating:
            # Calculate marginal likelihood using path sampling
            C = self.nchains
            marginal_like = 0.0
            self.output('\nCalculation of marginal likelihood:')
            self.output('%12s%12s' % ('chain', 'avg. lnL'))
            for i,v in enumerate(self.path_sample):
                n = len(v)
                avg = sum(v)/float(n)
                self.output('%12d%12.5f' % (i, avg))
                if (i == 0) or (i == C - 1):
                    avg /= 2.0
                marginal_like += avg
            marginal_like /= float(C)
            self.output('  Marginal likelihood (path sampling method) = %f' % marginal_like)
            
            # Calculate marginal likelihood using harmonic mean method on cold chain
            sample_size = len(self.path_sample[0])
            min_lnL = min(self.path_sample[0])
            sum_diffs = 0.0
            for lnl in self.path_sample[0]:
                diff = lnl - min_lnL
                if diff < 500.0:
                    sum_diffs += math.exp(-diff)
                else:
                    self.output('warning: ignoring large diff (%f) in harmonic mean calculation' % diff)
            log_harmonic_mean = math.log(sample_size) + min_lnL - math.log(sum_diffs)
            self.output('  Marginal likelihood(harmonic mean method)= %f' % log_harmonic_mean)

    def runSAMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the SSAMC analysis.
        
        """
        if self.samc_move_debug:
            stop_at_cycle = -1  # set this to the cycle you want to examine in detail, or -1 to ignore

        current_level = 0
        prev_level = 0
        addseq = self.addition_sequence
        chain = self.mcmc_manager.getColdChain()
        mgr = self.mcmc_manager.getColdChainManager()
        ls = chain.larget_simon_move
        samc_move = chain.samc_move
        max_level = len(addseq) - 4
        num_levels = len(addseq) - 3
        counts = [0]*(max_level + 1)
        mgr.refreshLastLnPrior()
        mgr.refreshLastLnLike()
        ln_proposal_ratio = 0.0 # always 0.0 because \tilde{T}_{m,m-1} = \tilde{T}_{m,m+1} = 1/3
        ls_accepted = [0]*num_levels
        ls_tried = [0]*num_levels
        if self.samc_move_debug:
            outvect = None
            if os.path.exists('cf.txt'):
                outvect = file('cf.txt', 'r').readlines()
            else:
                outfile = file('cf.txt', 'w')

        # Start of the main SAMC loop        
        for cycle in xrange(self.ncycles):

            if self.samc_move_debug:
                chain.tree.debugMode(False)
                if cycle == stop_at_cycle:
                    chain.tree.debugMode(True)
                    raw_input('Stopped at beginning of cycle %d. Press enter to continue...' % stop_at_cycle)

            # proposal for changing current level
            u = chain.r.uniform()
            proposed_level = current_level
            if u < 1.0/3.0:
                if current_level > 0:
                    proposed_level = current_level - 1
            elif u > 2.0/3.0:
                if current_level < max_level:
                    proposed_level = current_level + 1
                    
            if proposed_level == current_level:
                c = "*"
                samc_move.setSaveDebugInfo(True)
                if ls.update():
                    ls_accepted[current_level] += 1
                    ls_tried[current_level] += 1
                    c = "*"
                else:
                    ls_tried[current_level] += 1
                    c = "* :("
            elif proposed_level > current_level:
                leaf_num = self.addition_sequence[proposed_level + 3]
                theta_diff = self.samc_theta[current_level] - self.samc_theta[proposed_level]
                if samc_move.extrapolate(leaf_num, theta_diff, ln_proposal_ratio):
                    current_level = proposed_level
                    c = "+"
                else:
                    c = "+ :("                    
            else:
                leaf_num = self.addition_sequence[current_level + 3]
                theta_diff = self.samc_theta[current_level] - self.samc_theta[proposed_level]
                if samc_move.project(leaf_num, theta_diff, ln_proposal_ratio):
                    current_level = proposed_level
                    c = "-"
                else:
                    c = "- :("                    

            # Bump up the normalizing constant for the current level
            gain_factor = self.samc_t0/max([self.samc_t0, float(cycle)])
            self.samc_theta[current_level] += gain_factor
            lnL = chain.calcLnLikelihood()

            self.mcmc_manager.recordSample(cycle)

            outstr = "cycle = %d, level = %d: lnL =%*f %s" % (cycle, current_level, 15*current_level, lnL, c)
            print outstr

            if self.samc_move_debug:
                #if chain.debugCheckForUncachedCLAs():
                #    sys.exit('cached CLAs found at cycle %d' % cycle)
                if outvect:
                    if outvect[cycle].strip() != outstr:
                        print outvect[cycle].strip(),' <-- expected'
                        sys.exit('output differs from expected at cycle %d' % cycle)
                else:
                    outfile.write('%s\n' % outstr)

            prev_level = current_level
            counts[current_level] += 1
            
        if self.samc_move_debug and not outvect:
            outfile.close()
            
        print "ls accept pct =", [n > 0 and (100.0*float(a)/float(n) or 0) for a,n in zip(ls_accepted,ls_tried)]
        print "theta vector = ", self.samc_theta
        print "counts = ", counts
        avg_count = float(sum(counts))/float(len(counts))
        normalized_counts = [100.0*float(c)/avg_count for c in counts]
        print "normalized_counts = [%.1f" % normalized_counts[0],
        for c in normalized_counts[1:]:
            print ', %.1f' % c,
        print ']'
    
    def readDataFromFile(self):
        if not self.file_name_data_stored or (self.data_file_name != self.file_name_data_stored):
            self.reader.readFile(self.data_file_name)
            self.taxon_labels = self.reader.getTaxLabels()
            self.data_matrix = getDiscreteMatrix(self.reader, 0)
            self.ntax = self.data_matrix.getNTax()
            self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
            self.file_name_data_stored = self.data_file_name    # prevents rereading same data file later

    def distanceToTaxaSet(self, leaf, taxon_set, d):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Finds the minimum distance between a leaf and the members of a set of
        taxa.
        """
        
        min_dist = float('1e3000')
        for i in taxon_set:
            d_i = d[leaf][i]
            if d_i < min_dist:
                min_dist = d_i
        return min_dist
        
    def getSamcAdditionSequence(self):
        
        n = self.ntax
        d = self.calcDistances()
        self.samc_distance_matrix = d

        # Compute the addition order.  First find the most distant pair.
        max_dist = -1.0
        for i in range(n):
            for j in range(i+1,n):
                dij = d[i][j]
                if dij > max_dist:
                    max_dist = dij
                    self.addition_sequence = [i,j]

        # level_set = [self.addition_sequence]
        available_set = range(n)
        available_set.remove(self.addition_sequence[0])
        available_set.remove(self.addition_sequence[1])
        # print "level 2 self.addition_sequence =>", self.addition_sequence

        for level in range(2,n):
            max_dist = -1.0
            for i in available_set:
                d_to_set = self.distanceToTaxaSet(i, self.addition_sequence, d)
                if d_to_set > max_dist:
                    max_dist = d_to_set
                    next_i = i
            self.addition_sequence.append(next_i)
            available_set.remove(next_i)
            print "level", level, " self.addition_sequence =>", self.addition_sequence
            
    def getSamcStartingTree(self):
        #POL TEMP: not the best way to deal with this. See LikelihoodCore.__init__()
        r = Lot()
        r.setSeed(int(self.random_seed))
        self.starting_edgelen_dist.setLot(r)

        addseq = self.addition_sequence
        brlens = [self.starting_edgelen_dist.sample() for i in range(5)]
        self.tree_topology = "((%d:%.5f,%d:%.5f):%.5f,%d:%.5f,%d:%.5f)" % (
                             addseq[0] + 1, brlens[0], addseq[1] + 1, brlens[1], brlens[4], addseq[2] + 1,
                             brlens[2], addseq[3] + 1, brlens[3])
        self.starting_tree_source = 'usertree'
        self.starting_tree = self.tree_topology
        print "starting tree = ", self.tree_topology
        
    def setupSAMC(self):
        # Read the data
        self.phycassert(self.data_source == 'file', "This only works for data read from a file (specify data_source == 'file')")
        self.readDataFromFile()
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")

        if self.nchains != 1:
            print "Note: nchains reset to 1 for SSAMC"
            self.nchains = 1

        if not self.addition_sequence:
            self.getSamcAdditionSequence()
            
        if not self.tree_topology:
            self.getSamcStartingTree()

        # Initialize vector of normalizing constant parameters
        # Subtract 3 because never examine levels in which number of taxa < 4
        nlevels = len(self.addition_sequence) - 3
        self.samc_theta = [0.0]*(nlevels)
        
        # Set up MCMC chain
        self.heat_vector = [1.0]
        self.mcmc_manager.createChains()

        # TODO: createChains, when it calls prepareForLikelihood, should add all tips
        # and all internal nodes to vectors (not stacks) and use pointers to the appropriate
        # elements in the tree itself

        # Create tip nodes not already in tree and add them to tree's
        # tip node storage vector in reverse order
        internal_node_number = self.ntax
        m = self.mcmc_manager.getColdChain()
        for row in self.addition_sequence[-1:3:-1]:
            m.likelihood.addOrphanTip(m.tree, row, '%d' % (row+1))
            m.likelihood.addDecoratedInternalNode(m.tree, internal_node_number)
            internal_node_number += 1

        # samc_move needs a copy of the distance matrix for purposes of selecting
        # a node to which to join the next taxon in the addition sequence in
        # extrapolate moves
        for row in self.samc_distance_matrix:
            m.samc_move.setDistanceMatrixRow(row)

        # Set the number of taxa and temperature
        m.samc_move.setNTax(self.ntax)
        m.samc_move.setTemperature(self.samc_temperature)

        self.openParameterAndTreeFiles()
        
    def setLogFile(self, filename):
        # Open a log file if requested
        if self.logf:
            self.logf.close()
            self.logf = None
        if not filename:
            self._logFileName = None
        else:
            # TODO check first to see if it exists before blindly overwriting
            self.logf = file(filename, 'w')
            self._logFileName = filename

    def getLogFile(self):
        return self._logFileName
        
    log_file_name = property(getLogFile, setLogFile)
    
    def openParameterAndTreeFiles(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates parameter and tree file names based on the data file name or the
        user-supplied prefix and opens the files
        
        """
        prefix = os.path.abspath(self.data_file_name) #os.path.basename(self.data_file_name)
        if self.outfile_prefix:
            prefix = self.outfile_prefix
        self.param_file_name = prefix + '.p'
        self.tree_file_name = prefix + '.t'

        self.paramFileOpen()
        self.treeFileOpen()

    
    def setupMCMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is for parts of the setup that should occur right before
        runMCMC() is called. Setup is deferred until this point to give the
        user a chance to change the default settings before the call to
        runMCMC(). A call to the MCMCManager's createChains function is the
        last thing done by this function.
        
        """

        # Read the data
        #self.phycassert(self.data_source == 'file', "This only works for data read from a file (specify data_source == 'file')")
        if self.data_source == 'file':
            self.readDataFromFile()
                # Store (and create, if necessary) list of taxon labels
        elif (len(self.taxon_labels) != self.ntax):
            for i in range(self.ntax):
                s = 'taxon_%d' % (i + 1)
                self.taxon_labels.append(s)

        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")

        # Create a tree description to be used for building starting trees (formerly Phycas.setupTree function)
        if self.starting_tree_source == 'file':
            self.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")

            # Grab first tree description in the data file
            # TODO allow some other tree than the first
            newicks = []
            for t in self.reader.getTrees():
                newicks.append(t.newick)
            self.phycassert(len(newicks) > 0, 'a trees block defining at least one tree must be stored in the nexus data file')
            self.starting_tree = newicks[0]
        elif self.starting_tree_source == 'usertree':
            self.starting_tree = self.tree_topology
        elif self.starting_tree_source == 'random':
            self.phycassert(self.ntax > 0, 'expecting ntax to be greater than 0')
            self.starting_tree = None
        else:
            self.phycassert(False, "starting_tree_source should equal 'random', 'file', or 'usertree', but instead it was this: %s" % self.starting_tree_source)
        
        # Determine heating levels if multiple chains
        if self.heat_vector == None:
            if self.nchains == 1:
                self.heat_vector = [1.0]
            else:
                # Create a list for each chain to hold sampled lnL values
                self.path_sample = []
                for i in range(self.nchains):
                    self.path_sample.append([])
                    
                self.heat_vector = []
                if self.is_standard_heating:
                    for i in range(self.nchains):
                        # Standard heating 
                        # 0 1.000 = 1/1.0 cold chain explores posterior
                        # 1 0.833 = 1/1.2
                        # 2 0.714 = 1/1.4
                        # 3 0.625 = 1/1.6
                        temp = 1.0/(1.0 + float(i)*self.heating_lambda)
                        self.heat_vector.append(temp)
                else:
                    for i in range(self.nchains):
                        self.do_marginal_like = True
                        # Likelihood heating for thermodynamic integration
                        # 0 1.000 = (3-0)/3 cold chain explores posterior
                        # 1 0.667 = (3-1)/3
                        # 2 0.333 = (3-2)/3
                        # 3 0.000 = (3-3)/3 hottest chain explores prior
                        temp = 1.0
                        denom = float(self.nchains - 1)
                        temp = float(self.nchains - i - 1)/denom
                        self.heat_vector.append(temp)
        else:
            # User supplied heat_vector; check to make sure they allowed for a cold chain
            self.nchains = len(self.heat_vector)
            self.phycassert(self.heat_vector.index(1.0) < self.nchains, 'no cold chain was specified')

        # Call MCMCManager's createChains function
        self.mcmc_manager.createChains()
        
        self.openParameterAndTreeFiles()
        
    def mcmc(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs an MCMC analysis.
        
        """
        self.setupMCMC()
        self.runMCMC()

    def samc(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a sequential stochastic Markov chain analysis.
        
        """
        self.doing_samc = True;
        self.setupSAMC()
        self.runSAMC()

    def simulateDNA(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simulates a DNA dataset and stores it in NEXUS format in the supplied
        filename.
        
        """
        self.starting_tree = self.tree_topology
        self.phycassert(self.data_source == None, 'set data_source to None before calling simulateDNA')
        self.ntax = len(self.sim_taxon_labels)
        core = MCMCManager.LikelihoodCore(self)
        core.setupCore()
        if not core.tree.hasEdgeLens():
            tm = TreeManip(core.tree)
            tm.setRandomEdgeLengths(core.starting_edgelen_dist)
        self.sim_model_tree = core.tree
        core.prepareForSimulation()
        sim_data = core.simulate()
        sim_data.saveToNexusFile(self.sim_file_name, self.sim_taxon_labels, 'dna', ('a','c','g','t'))
        
    def likelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current tree and current
        model.
        
        """
        self.starting_tree = self.tree_topology
        self.phycassert(self.data_source == 'file', "set data_source to 'file' and specify data_file_name before calling the likelihood function")
        self.readDataFromFile()
        core = MCMCManager.LikelihoodCore(self)
        core.setupCore()
        core.prepareForLikelihood()
        return core.calcLnLikelihood()
        
    def gg(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes Gelfand-Ghosh on a pre-existing MCMC sample defined in the
        files self.gg_pfile and self.gg_tfile.
        
        """
        self.phycassert(self.gg_pfile, 'gg_pfile cannot be None if gg function called')
        self.phycassert(self.gg_tfile, 'gg_pfile cannot be None if gg function called')
        import GGImpl
        gelfand_ghosh = GGImpl.GelfandGhosh(self)
        self.gg_Pm, self.gg_Gm, self.gg_Dm = gelfand_ghosh.run()
        return (self.gg_Pm, self.gg_Gm, self.gg_Dm)

    def calcDistances(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes a matrix of pairwise distances between sequences.
        
        """
        self.phycassert(self.data_source == 'file', "data_source variable must equal 'file'")

        # matrix holds the data matrix
        matrix = []
        for i in range(self.ntax):
            matrix.append(self.data_matrix.getRow(i))

        # distance holds the JC69 distance matrix
        # Note: currently assumes no missing or ambiguous data (i.e. A vs. ? treated as difference)
        distance = []
        for i in range(self.ntax):
            distance_i = []
            for j in range(self.ntax):
                diff = 0.0
                if i < j:
                    diff = 0.0
                    for k in range(self.nchar):
                        xik = matrix[i][k]
                        xjk = matrix[j][k]
                        if xik != xjk:
                            diff += 1.0
                    p = diff/float(self.nchar)
                    if p > 0.74999:
                        p = 0.74999
                    jc = -0.75*math.log(1.0 - 4.0*p/3.0)
                    distance_i.append(jc)
                    #print 'saving distance[%d][%d] = %f' % (i,j,jc)
                elif i > j:
                    #print 'accessing [%d][%d] = %f' % (j,i,distance[j][i])
                    distance_i.append(distance[j][i])
                else:
                    distance_i.append(0.0)
            distance.append(distance_i)
            
#         for i in range(self.ntax):
#             for j in range(self.ntax):
#                 print "%10f" % distance[i][j],
#             print

        return distance

    # by default phycassert sys.exit.
    # When debugging, it is nice to set this to True so that you can see the stack trace
    PhycassertRaisesException = False
    CPPCompiledInDebug = False

if __name__ == '__main__':
    print "The Phycas.py file should be imported, not run directly. To import it,"
    print "create a python script (i.e. a file with a name ending in .py) with the"
    print "following text:"
    print
    print "from phycas import *"
    print "myphycas = Phycas()"
    print "..."
    print
    print "See examples in the Examples and Tests folders for more information, or"
    print "consult the PDF manual."
    print
    raw_input('Press the enter key when finished reading the lecture')
