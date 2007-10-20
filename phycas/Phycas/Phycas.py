import os, sys, math, threading
import MCMCManager  # poorly named, as MCMCManager is now only one of many classes within
from phycas.Conversions import *
from phycas.DataMatrix import *
from phycas.Likelihood import *
from phycas.Phylogeny import *
from phycas.ProbDist import *
from phycas.ReadNexus import *

class Phycas:
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
        self.quiet                  = False     # if True, output will only be send to self.logfile, if defined (see below); if False, output will be sent to the console as well
        self.logfile                = None      # specify a filename to save output to file 
        self.outfile_prefix         = None      # If None, parameter and tree files created will have a name beginning with the name of the data file; if provided, this prefix will form the first part of the parameter and tree file names
        
        # Variables associated with substitution models (except for edge lengths)
        self.default_model          = 'hky'     # Can be 'jc', 'hky' or 'gtr'

        self.relrate_prior          = ExponentialDist(1.0)
        self.starting_relrates      = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
        self.fix_relrates           = False # new

        self.kappa_prior            = ExponentialDist(1.0) # new
        self.starting_kappa         = 4.0
        self.fix_kappa              = False # new

        self.num_rates              = 1         # default is no rate heterogeneity (i.e. 1 rate)
        self.gamma_shape_prior      = ExponentialDist(1.0)
        self.starting_shape         = 0.5
        self.fix_shape              = False # new
        self.use_inverse_shape      = False      # if True, gamma_shape_prior is applied to 1/shape rather than shape

        self.estimate_pinvar        = False
        self.pinvar_prior           = BetaDist(1.0, 1.0)
        self.starting_pinvar        = 0.2
        self.fix_pinvar             = False # new

        self.base_freq_param_prior  = ExponentialDist(1.0)
        self.starting_freqs         = [1.0, 1.0, 1.0, 1.0] # new
        self.fix_freqs              = False # new

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
        self.starting_edgelen_hyperparam = 0.05 # new  #POL doesn't do anything! currently ignored
        self.fix_edgelen_hyperparam = False # new
        
        self.external_edgelen_dist  = ExponentialDist(2.0)
        self.internal_edgelen_dist  = None
        self.fix_edgelens           = False # new
        
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
        #self.r                      = Lot()     #POL still being used?
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
        self.separate_int_ext_edgelen_priors = False
        #self.use_tree_viewer        = False    # popup graphical TreeViewer to show trees during run POLPY_NEWWAY

    def phycassert(self, assumption, msg):
        if not assumption:
            sys.exit('Error: ' + msg)

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

    #def readNexusFile(self, fn):            
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Sets self.data_file_name to supplied filename fn, calls the
    #    self.reader.readFile function, storing the data matrix in 
    #    self.data_matrix. Also sets self.ntax and self.nchar accordingly.
    #    
    #    """
    #    assert self.data_source == 'file', "set data_source to 'file' before calling readNexusFile"
    #    if not self.data_file_name == fn:
    #        self.data_file_name = fn
    #    self.reader.readFile(self.data_file_name)
    #    self.data_matrix = getDiscreteMatrix(self.reader, 0)
    #    self.ntax = self.data_matrix.getNTax()
    #    self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
    #    # used to avoid next two lines if self.data_source was None, but we have to
    #    # copy over the data from data_matrix to set self.npatterns, so without the
    #    # next two lines, we would not be able to run the pattern-specific rates model
    #    # without data because we would not know how many rate parameters to create
    #    assert self.likelihood, 'call Phycas.setupLikelihood before calling Phycas.readNexusFile'
    #    self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
    #    self.npatterns = self.likelihood.getNPatterns()

    #def setupTree(self, source = None):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    The purpose of setupTree is to create a string description of the
    #    starting tree that can be handed to each chain for purposes of
    #    building the starting tree. If starting_tree_source equals 'file',
    #    then the first newick tree description that appears in the data file
    #    is used. In this case, the program will abort if there are no trees
    #    stored in the specified data file. If starting_tree_source is equal
    #    to 'random', then a random starting tree will be used (with edge
    #    lengths drawn from self.starting_edgelen_dist. If starting_tree_source
    #    is 'usertree', it is assumed that self.tree_topology contains a
    #    valid newick tree description.
    #    
    #    """
    #    if source:
    #        self.starting_tree_source = source
    #        
    #    # If user requested using a tree from the data file, grab the first one stored there
    #    self.tree = Tree()
    #    if self.starting_tree_source == 'file':
    #        assert self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)"
    #        
    #        # Grab first tree description in the data file
    #        newicks = []
    #        for t in self.reader.getTrees():
    #            newicks.append(t.newick)
    #        assert len(newicks) > 0, 'Error: a trees block defining at least one tree must be stored in the nexus data file'
    #        self.starting_tree = newicks[0]
    #
    #        # Build a Tree object from the description stored in the data file
    #        self.tree.buildFromString(self.starting_tree)
    #        if not self.tree.tipNumbersSetUsingNames():
    #            self.warn_tip_numbers = True
    #        
    #    elif self.starting_tree_source == 'usertree':
    #        # self.tree_topology should already be created
    #        self.starting_tree = self.tree_topology
    #
    #        # Build a Tree object from the description stored in self.starting_tree
    #        self.tree.buildFromString(self.starting_tree)
    #        if not self.tree.tipNumbersSetUsingNames():
    #            self.warn_tip_numbers = True
    #        
    #    elif self.starting_tree_source == 'random':
    #        assert self.ntax > 0, 'expecting ntax to be greater than 0'
    #        
    #        # Build a random tree
    #        self.starting_edgelen_dist.setLot(self.r)
    #        TreeManip(self.tree).randomTree(
    #            self.ntax,     # number of tips
    #            self.r,        # pseudorandom number generator
    #            self.starting_edgelen_dist, # distribution from which to draw starting edge lengths
    #            False)         # Yule tree if True, edge lengths independent if False
    #        self.starting_tree = self.tree.makeNewick()
    #        self.warn_tip_numbers = False
    #        
    #    else:
    #        # throw exception
    #        assert False, 'starting_tree_source should equal random, file, or usertree, but instead it was this: %s' % self.starting_tree_source
    #
    #    # Make sure that names of tips equal the string equivalent of the tip node number plus 1
    #    # This means that when tree topologies are sampled, the tree definitions that are output
    #    # to the .t file will be similar to MrBayes output
    #    assert self.ntax > 0, 'expecting ntax to be greater than 0'
    #    self.nedges = 2*self.ntax - 3
    #    taxNames = []
    #    for i in range(self.ntax):
    #        taxNames.append(str(i+1))
    #    self.tree.rectifyNames(taxNames)

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

    #def recordSample(self, cycle, lnL = 0.0):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Records current tree topology and edge lengths by adding a line to
    #    the tree file, and records tree length and substitution parameters
    #    by adding a line to the parameter file.
    #    
    #    """
    #
    #    # Add line to parameter file if it exists
    #    if self.paramf:
    #        self.paramf.write('%d\t%.3f\t%.3f' % (cycle + 1, lnL, self.tree.edgeLenSum()))
    #        self.paramf.write(self.model.paramReport())
    #        if self.using_hyperprior:
    #            self.chain_manager = self.mcmc_manager.getColdChainManager()
    #            for p in self.chain_manager.getEdgeLenHyperparams():
    #                #p = self.chain_manager.getEdgeLenHyperparam()
    #                self.paramf.write('\t%.5f' % p.getCurrValue())
    #        if self.use_flex_model:
    #            rates_vector = self.likelihood.getRateMeans()
    #            for rr in rates_vector:
    #                self.paramf.write('\t%.5f' % rr)
    #                #if cycle == 0:
    #                #    print '  rate = %f' % rr 
    #            probs_vector = self.likelihood.getRateProbs()
    #            for rp in probs_vector:
    #                self.paramf.write('\t%.5f' % rp)
    #                #if cycle == 0:
    #                #    print '  prob = %f' % rp 
    #        self.paramf.write('\n')
    #    
    #    # Add line to tree file if it exists
    #    if self.treef:
    #        self.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, self.tree.makeNewick()))

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
        #assert self.metropolis_weight == 0, 'deprecated variable metropolis_weight was used'
        
        if self.verbose:
            if self.data_source == None:
                self.output('Data source:    None (running MCMC with no data to explore prior)')
            elif self.data_source == 'file':
                self.output('Data source:    %s' % self.data_file_name)
            #elif self.data_source == 'memory':
            #    self.output('Data source:    Data already in memory')
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
        #self.likelihood.resetNEvals()

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

        for cycle in range(self.ncycles):
            for i,c in enumerate(self.mcmc_manager.chains):
                #tmpf = file('debug_info.txt', 'a') # comment out for release
                tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,i))
                for p in c.chain_manager.getAllUpdaters():
                    w = p.getWeight()
                    for x in range(w):
                        #p.setSaveDebugInfo(True)  # comment out for release
                        if cycle == 3 and p.getName() == 'Discrete gamma shape':
                            #raw_input('cycle 3, Discrete gamma shape')
                            self.mcmc_manager.getColdChain().likelihood.setDebug(True)
                        # print p.getName()
                        p.update()
                        #tmpf.write('%s | %s\n' % (p.getName(), p.getDebugInfo()))  # comment out for release
                #tmpf.close()  # comment out for release
                # sys.exit('debug kill')  # comment out for release
            if self.verbose and (cycle + 1) % self.report_every == 0:
                self.stopwatch.normalize()
                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                msg = 'cycle = %d, lnL = %.5f (%.5f secs)' % (cycle + 1, cold_chain_manager.getLastLnLike(), self.stopwatch.elapsedSeconds())
                #if self.use_flex_model:
                #    bytes_per_cla = self.likelihood.bytesPerCLA()
                #    ncreated = self.likelihood.numCLAsCreated()
                #    megabytes = ncreated*bytes_per_cla//1048576
                #    msg += ', %d rates, %d MB in %d CLAs' % (self.model.getNGammaRates(),megabytes,ncreated)
                self.output(msg)
                #if self.use_tree_viewer and self.tree_viewer:
                #    self.tree_viewer.refresh('Cycle %d' % (cycle + 1))
            if (cycle + 1) % self.sample_every == 0:
                self.mcmc_manager.recordSample(cycle)
                self.stopwatch.normalize()
            if (cycle + 1) % next_adaptation == 0:
                self.adaptSliceSamplers()
                next_adaptation += 2*(next_adaptation - last_adaptation)
                last_adaptation = cycle + 1

        #if self.use_tree_viewer and self.tree_viewer:
        #    self.tree_viewer.close() # POLPY_NEWWAY

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

    def readDataFromFile(self):
        if not self.file_name_data_stored or (self.data_file_name != self.file_name_data_stored):
            self.reader.readFile(self.data_file_name)
            self.taxon_labels = self.reader.getTaxLabels()
            self.data_matrix = getDiscreteMatrix(self.reader, 0)
            self.ntax = self.data_matrix.getNTax()
            self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
            self.file_name_data_stored = self.data_file_name    # prevents rereading same data file later

    def setupMCMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is for parts of the setup that should occur right before
        run() is called. Setup is deferred until this point to give the user
        a chance to change the default settings before the call to run(). A
        call to the MCMCManager's createChains function is the last thing done
        by this function.
        
        """
        # Open a logfile if requested
        if self.logfile:
            # TODO check first to see if it exists before blindly overwriting
            self.logf = file(self.logfile, 'w')

        # Determine the number of edge length priors
        self.phycassert(self.external_edgelen_dist or self.internal_edgelen_dist, 'external_edgelen_dist and internal_edgelen_dist are both None; a probability distribution must be assigned to at least one of these')
        self.separate_int_ext_edgelen_priors = True
        if not self.external_edgelen_dist:
            # if only one edgelen prior distribution specified, make sure it is external_edgelen_dist
            self.external_edgelen_dist = self.internal_edgelen_dist
            self.internal_edgelen_dist = None
            self.separate_int_ext_edgelen_priors = False
        elif not self.internal_edgelen_dist:
            self.separate_int_ext_edgelen_priors = False

        #if self.data_source == 'memory':
        #    assert self.__dict__.get('likelihood') is not None, "data_source == 'memory' implies that likelihood object exists"
        #    self.likelihood.replaceModel(self.model)
        #else:
        #    self.setupLikelihood() # now gone, moved to MCMCChain.setupModel

        # Read the data
        if self.data_source == 'file':
            self.readDataFromFile()
            #self.reader.readFile(self.data_file_name)
            #self.taxon_labels = self.reader.getTaxLabels()
            #self.data_matrix = getDiscreteMatrix(self.reader, 0)
            #self.ntax = self.data_matrix.getNTax()
            #self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
            #self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
        #elif self.data_source == 'memory':
        #    assert self.ntax > 0
        #    assert self.nchar > 0
            
        #self.npatterns = self.likelihood.getNPatterns()

        # Create a tree description to be used for building starting trees (formerly Phycas.setupTree function)
        if self.starting_tree_source == 'file':
            self.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")

            # Grab first tree description in the data file
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
            self.phycassert(False, 'starting_tree_source should equal random, file, or usertree, but instead it was this: %s' % self.starting_tree_source)
        
        #self.likelihood.prepareForLikelihood(self.tree)

        # Determine heating levels if multiple chains
        if self.heat_vector == None:
            if self.nchains == 1:
                self.heat_vector = [1.0]
            else:
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
                        self.self.do_marginal_like = True
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
        
        # Create parameter and tree file names based on the data file name or the user-supplied prefix
        prefix = os.path.abspath(self.data_file_name) #os.path.basename(self.data_file_name)
        if self.outfile_prefix:
            prefix = self.outfile_prefix
        self.param_file_name = prefix + '.p'
        self.tree_file_name = prefix + '.t'

        # Open the parameter file
        self.paramFileOpen()
        
        # Store (and create, if necessary) list of taxon labels
        if (not self.data_source) or (len(self.taxon_labels) != self.ntax):
            for i in range(self.ntax):
                s = 'taxon_%d' % (i + 1)
                self.taxon_labels.append(s)

        # Open the tree file
        self.treeFileOpen()

    def mcmc(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs an MCMC analysis.
        
        """
        self.setupMCMC()
        self.runMCMC()

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
        import GGImpl
        gelfand_ghosh = GGImpl.GelfandGhosh(self)
        self.gg_Pm, self.gg_Gm, self.gg_Dm = gelfand_ghosh.run()
        return (self.gg_Pm, self.gg_Gm, self.gg_Dm)
