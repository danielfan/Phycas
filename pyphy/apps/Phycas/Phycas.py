import os, sys, time, math
import os, sys, time, math
import DataMatrix
import ReadNexus
import Phylogeny
import Likelihood
import ProbDist

class Phycas:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Performs a phylogenetic MCMC analysis. The tree topology and edge
    lengths are updated via the Metropolis-Hastings algorithm (using the
    Larget-Simon LOCAL move without a molecular clock). Slice sampling is
    used to update all model parameters except edge lengths. A series of
    Larget-Simon updates are periodically interrupted by a round of
    slice sampling in which each model parameter is updated. The period
    is determined by the variables metropolis_weight and slice_weight. If
    the default values 9 and 1 are used for metropolis_weight and
    slice_weight, respectively, then a single cycle will consist of 9
    Larget-Simon metropolis updates followed by one round of slice
    sampling.

    The following example reads a nexus data file, specifies a tree
    topology and runs a single Markov chain sampler for 100 cycles using
    the HKY model.
    
    >>> import Phycas
    >>> mcmc = Phycas.Phycas()
    >>> mcmc.data_file_name = 'nyldna4.nex'
    >>> mcmc.metropolis_weight = 9
    >>> mcmc.slice_weight = 1
    >>> mcmc.ncycles = 100
    >>> mcmc.sample_every = 100
    >>> mcmc.report_every = 10
    >>> mcmc.adapt_first = 25
    >>> mcmc.random_seed = '13579'
    >>> mcmc.model_type = 'hky'
    >>> mcmc.verbose = True
    >>> mcmc.gg_do = False
    >>> mcmc.run() # doctest:+ELLIPSIS
    Data file:      nyldna4.nex
    Prior:          hyper
    No. cycles:     100
    Sample every:   100
    Tree topology:  (0:0.1,1:0.1,(2:0.1,3:0.1):0.1)
    No. samples:    1
    Sampled trees will be saved in nyldna4.nex.t
    Sampled parameters will be saved in nyldna4.nex.p
    Tip node numbers were set using the names in the tree description
    <BLANKLINE>
    Sampling (100 cycles)...
    <BLANKLINE>
    cycle = 10, lnL = -7129.98700
    cycle = 20, lnL = -7128.98328
    <BLANKLINE>
    Slice sampler diagnostics:
      mode=0.05148, avgevals=10.440 (edge length for node 1)
      mode=0.06327, avgevals=10.120 (edge length for node 1)
      mode=0.00804, avgevals=10.080 (edge length for node 0)
      mode=0.03692, avgevals=11.520 (edge length for node 2)
      mode=0.06462, avgevals=8.880 (edge length for node 3)
      mode=0.14666, avgevals=6.920 (edge length hyperprior)
      mode=1.90334, avgevals=7.360 (trs/trv rate ratio)
      mode=10.13449, avgevals=8.560 (base freq. A)
      mode=6.47742, avgevals=6.640 (base freq. C)
      mode=7.40993, avgevals=7.000 (base freq. G)
      mode=11.88085, avgevals=8.840 (base freq. T)
    <BLANKLINE>
    cycle = 30, lnL = -7128.54884
    cycle = 40, lnL = -7126.97586
    cycle = 50, lnL = -7128.50749
    cycle = 60, lnL = -7128.66949
    cycle = 70, lnL = -7129.21719
    <BLANKLINE>
    Slice sampler diagnostics:
      mode=0.04988, avgevals=6.080 (edge length for node 1)
      mode=0.06189, avgevals=6.120 (edge length for node 1)
      mode=0.01083, avgevals=8.120 (edge length for node 0)
      mode=0.04150, avgevals=5.960 (edge length for node 2)
      mode=0.06003, avgevals=6.220 (edge length for node 3)
      mode=0.17108, avgevals=6.040 (edge length hyperprior)
      mode=1.92832, avgevals=6.620 (trs/trv rate ratio)
      mode=4.17365, avgevals=6.440 (base freq. A)
      mode=2.46462, avgevals=6.180 (base freq. C)
      mode=3.15126, avgevals=6.160 (base freq. G)
      mode=4.97198, avgevals=6.500 (base freq. T)
    <BLANKLINE>
    cycle = 80, lnL = -7130.20819
    cycle = 90, lnL = -7126.20187
    cycle = 100, lnL = -7131.13567
    ...

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the object with some default values. If these defaults do
        not suit, they can be changed before run() is called.

        Starting Tree. The variable starting_tree_source is equal to 'random'
        by default, which mean that the program will generate a random tree
        topology and populate it with exponentially distributed (mean
        edgelen_prior_mean) edge lengths. If starting_tree_source is 'file',
        then the program expects to find at least one tree topology defined
        in the nexus data file. If starting_tree_source is 'usertree', then
        the program expects the variable tree_topology to contain a correct
        newick tree description.
        
        Edge Length Hyperprior. The variable using_hyperprior is True by
        default, which means that the mean of the edge length prior
        distribution is governed by a hyperparameter. If using_hyperprior is
        changed to False before run() is called, then edgelen_prior_mean
        will be used as the mean of an exponential prior distribution for
        edge lengths.

        Substitution Model. The HKY model is used by default; set model_type
        to 'jc' to use the JC model instead.

        Pseudorandom Number Generation. The variable random_seed is 'auto'
        by default, which means that the pseudorandom number generator is
        initialized using the system clock. You can, however, set random_seed
        to a particular random number seed to repeat a previous analysis
        exactly.

        Other quantities that typically need to be specified include ncycles,
        sample_every, report_every, adapt_first and data_file_name.
        
        """
        self.starting_tree_source   = 'random'  # source of starting tree topology
        self.starting_tree          = ''        # will contain description of actual starting tree used
        self.tree_topology          = ''        # unused unless starting_tree_source is 'usertree'
        self.using_hyperprior       = True      # Hyperprior for edge length prior used if True
        self.edgelen_prior_mean     = 0.5       # Prior mean of exponential edge length distribution
        self.model_type             = 'hky'     # HKY model used if True, JC used if False
        self.verbose                = True      # more output if True
        self.random_seed            = 'auto'    # determines random number seed
        self.metropolis_weight      = 9         # Metropolis moves will be performed this many times per cycle
        self.slice_weight           = 1         # Slice sampled parameters will be updated this many times per cycle
        self.slice_max_units        = 1000      # Max. number of units used in slice sampling

        # Related to the data file        
        self.data_file_name         = ''        # will hold actual data file name
        self.ntax                   = 0         # will hold the actual number of taxa after data file read
        self.nchar                  = 0         # will hold the actual number of characters after data file has been read
        self.npatterns              = 0         # will hold the actual number of patterns after data file has been read
        self.taxon_labels           = []        # will hold taxon labels from data file or default names if self.no_data is True

        self.no_data                = False     # if True, Phycas will pretend that there is no data during MCMC runs (useful for exploring the prior)

        # Settings related to the model and prior distributions
        self.num_rates              = 1         # default is rate homogeneity (1 rate)
        self.relrate_prior          = ProbDist.ExponentialDist(1.0)
        self.base_freq_param_prior  = ProbDist.ExponentialDist(1.0)
        self.gamma_shape_prior      = ProbDist.ExponentialDist(1.0)
        self.use_inverse_shape      = True  # if True, gamma_shape_prior applied to 1/shape rather than shape
        self.estimate_pinvar        = False
        self.pinvar_prior           = ProbDist.BetaDist(1.0, 1.0)
        self.use_flex_model         = False
        self.num_flex_spacers       = 1
        
        # MCMC settings (used by run function)
        self.ncycles                = 10000
        self.sample_every           = 100
        self.report_every           = self.ncycles//100
        self.paramf                 = None
        self.treef                  = None

        # Bush move settings
        self.allow_polytomies       = False     # if True, do Bush moves in addition to Larget-Simon moves; if False, do Larget-Simon moves only
        self.polytomy_prior         = True      # if True, use polytomy prior; if False, use resolution class prior
        self.topo_prior_C           = 2.0       # specifies the strength of the prior (C = 1 is flat prior; C > 1 favors less resolved topologies)
        self.bush_move_edgelen_mean = 1.0       # specifies mean of exponential edge length generation distribution used by BushMove when new edges are created
        
        # Simulation settings (used by simulate function)
        self.sim_nreps              = 0
        self.sim_outfile            = 'simout.nex'

        # Gelfand-Ghosh settings and variables
        self.gg_do                  = False         # gather GG statistics during MCMC run if True
        self.gg_outfile             = 'ggout.txt'   # file in which to save gg results (use None to not save results)
        self.gg_nreps               = 1             # the number of replicate simulations to do every MCMC sample
        self.gg_y                   = Likelihood.SimData()  # observed dataset
        self.gg_mu                  = Likelihood.SimData()  # mean of all posterior predictive datasets
        self.gg_a                   = []            # vector of compromise actions (one for each k in gg_kvect)
        self.gg_npatterns           = []            # vector containing the number of patterns in each posterior predictive dataset
        self.gg_t                   = []            # vector of t values computed from posterior predictive datasets
        self.gg_kvect               = [1.0]         # vector of k values to use when computing Gm and Dm
        self.gg_total               = 0
        self.gg_t_y                 = 0.0           # t for original dataset
        self.gg_t_mean              = 0.0           # mean of t over all posterior predictive datasets
        self.gg_t_mu                = 0.0           # t of mean over all posterior predictive datasets
        self.gg_t_a                 = []            # vector of t values computed from compromise action (one for each k in gg_kvect)
        self.gg_Gm                  = []            # vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Pm                  = 0.0           # penalty component (same for all k)
        self.gg_Dm                  = []            # vector of overall measures (one for each k in gg_kvect)
        self.gg_save_postpreds      = False         # if True, all posterior predictive data sets will be saved
        self.gg_postpred_prefix     = 'postpred'    # prefix to use for posterior predictive dataset filenames (only used if gg_save_postpreds is True)
        self.gg_save_spectra        = False         # adds all 256 counts for posterior predictive simulated data sets to a file named spectra.txt, with counts separated by tabs (only use for four-taxon problems)
        self.gg_spectrum            = Likelihood.SimData()  # workspace used if gg_save_spectra is True
        self.gg_spectrum_points     = ''            # used for creating surface plot in Maple for spectrum
        self.gg_spectrum_row        = 0

        # Adaptation of slice samplers is performed the first time at cycle adapt_first.
        # Subsequent adaptations wait twice the number of cycles as the previous adaptation.
        # Thus, adaptation n occurs at cycle adapt_first*(2^n - 1)
        # The total number of adaptations that will occur during an MCMC run is
        #      [ln(adapt_first + ncycles) - ln(adapt_first)]/ln(2)
        self.adapt_first            = 100

        # Slice sampler adaptation parameters
        #
        self.adapt_simple_param     = 0.5
        #self.adapt_ycond_param      = 1.3
        #self.adapt_ycond_from_ends  = 0.25
        
        # Create a pseudorandom number generator
        self.r = ProbDist.Lot()

        # Create a Nexus file reader
        self.reader = ReadNexus.NexusReader()

    def setupModel(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Defines prior distributions, creates parameters and adds these to the
        ParamManager object.
        
        """
        # Create a substitution model and define priors for the model parameters
        if self.model_type == 'gtr':
            self.model = Likelihood.GTRModel()
            self.model.setRelRates([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
            self.model.setRelRatePrior(self.relrate_prior)
            self.model.setBaseFreqParamPrior(self.base_freq_param_prior)
        elif self.model_type == 'hky':
            self.model = Likelihood.HKYModel()
            self.model.setKappa(4.0)
            self.model.setKappaPrior(self.relrate_prior)
            self.model.setBaseFreqParamPrior(self.base_freq_param_prior)
        else:
            self.model = Likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
        # Note must defer setting up pattern specific rates model until we know number of patterns
        if self.num_rates > 1 and self.use_flex_model:
            self.model.setNGammaRates(self.num_rates)
            self.model.setFlexModel()
            self.model.setNumFlexSpacers(self.num_flex_spacers)
            self.model.setFLEXProbParamPrior(ProbDist.ExponentialDist(1.0))
        elif self.num_rates > 1:
            self.model.setNGammaRates(self.num_rates)
            self.model.setPriorOnShapeInverse(self.use_inverse_shape)
            self.model.setShape(0.5)
            self.model.setDiscreteGammaShapePrior(self.gamma_shape_prior)
        else:
            self.model.setNGammaRates(1)
            
        if self.estimate_pinvar:
            assert not self.use_flex_model, 'Cannot currently use flex model with pinvar'
            self.model.setPinvarModel()
            self.model.setPinvar(0.2)
            self.model.setPinvarPrior(self.pinvar_prior)
        else:
            self.model.setNotPinvarModel()
        
        # Define an edge length prior distribution
        assert self.edgelen_prior_mean > 0.0, 'edgelen_prior_mean must be a positive, non-zero number'
        v = 1.0/float(self.edgelen_prior_mean)
        self.master_edgelen_dist = ProbDist.ExponentialDist(v);
        self.model.setEdgeLenPrior(self.master_edgelen_dist)
        if self.using_hyperprior:
            # Edge length prior governed by a InverseGamma-distributed hyperprior
            d = ProbDist.InverseGammaDist(2.1, 0.909)
            d.setMeanAndVariance(1.0, 10.0)
            self.model.setEdgeLenHyperPrior(d)
        else:
            # Edge length prior distribution is not hierarchical
            self.model.setEdgeLenHyperPrior(None)

        # Create a TreeLikelihood object. This can be used to compute the likelihood
        # for any tree based on the supplied data matrix and using the supplied model.

    def createChain(self):
        # Create a list of parameters for updating quantities such as kappa and the
        # base frequencies. 
        self.chain_manager = Likelihood.MCMCChainManager()
        self.chain_manager.addMCMCUpdaters(self.model,              # substitution model
                                           self.tree,               # tree
                                           self.likelihood,         # likelihood calculation machinery
                                           self.r,                  # pseudorandom number generator
                                           False,                   # separate_edgelen_params
                                           self.slice_max_units,    # weight for each parameter added
                                           self.slice_weight)       # weight for each parameter added

        # Create a LargetSimonMove object to handle Metropolis-Hastings
        # updates to the tree topology and edge lengths
        self.larget_simon_move = Likelihood.LargetSimonMove()
        self.larget_simon_move.setName("Larget-Simon move")
        self.larget_simon_move.setTree(self.tree)
        self.larget_simon_move.setModel(self.model)
        self.larget_simon_move.setTreeLikelihood(self.likelihood)
        self.larget_simon_move.setLot(self.r)
        self.larget_simon_move.setLambda(0.2) # should be a user setting
        if self.model.edgeLengthsFixed():
            self.larget_simon_move.fixParameter()
        self.chain_manager.addMove(self.larget_simon_move)

        # If requested, create a BushMove object to allow polytomous trees
        if self.allow_polytomies:
            # If allowing polytomies, do about half and half Larget Simon moves vs. Bush moves
            w = self.metropolis_weight//2
            self.larget_simon_move.setWeight(self.metropolis_weight - w)
            
            # Create a BushMove object
            self.bush_move = Likelihood.BushMove()

            # Set up the topology prior
            self.topo_prior_calculator = self.bush_move.getTopoPriorCalculator()
            self.topo_prior_calculator.chooseUnrooted()
            self.topo_prior_calculator.setC(self.topo_prior_C)
            if self.polytomy_prior:
                self.topo_prior_calculator.choosePolytomyPrior()
            else:
                self.topo_prior_calculator.chooseResolutionClassPrior()
                
            # Continue setting up BushMove object
            self.bush_move.setName("Bush move")
            self.bush_move.setWeight(w)
            self.bush_move.setTree(self.tree)
            self.bush_move.setModel(self.model)
            self.bush_move.setTreeLikelihood(self.likelihood)
            self.bush_move.setLot(self.r)
            self.bush_move.setEdgeLenDistMean(self.bush_move_edgelen_mean)
            if self.model.edgeLengthsFixed():
                self.bush_move.fixParameter()
            self.bush_move.finalize()
            self.chain_manager.addMove(self.bush_move)
        else:
            # Only Larget Simon moves if not allowing polytomies
            self.larget_simon_move.setWeight(self.metropolis_weight)

        self.chain_manager.finalize()
        
    def setupTree(self, source = None):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If starting_tree_source equals 'file', then a tree is built using the
        first newick tree description that appears in the data file. In this
        case, the program will abort if there are no trees stored in the
        specified data file. If starting_tree_source is equal to the string
        'random', then a random starting tree will be used (with edge lengths
        drawn from an exponential distribution having mean edgelen_prior_mean.
        If starting_tree_source is 'usertree', it is assumed that
        tree_topology contains a valid newick tree description, and a tree
        is built from this description.
        
        """
        if source:
            self.starting_tree_source = source
            
        # If user requested using a tree from the data file, grab the first one stored there
        self.tree = Phylogeny.Tree()
        if self.starting_tree_source == 'file':
            assert not self.no_data, 'Specified starting_tree_source = file when no_data = True (file was not read)'
            
            # Grab first tree description in the data file
            newicks = []
            for t in self.reader.getTrees():
                newicks.append(t.newick)
            assert len(newicks) > 0, 'Error: a trees block defining at least one tree must be stored in the nexus data file'
            self.starting_tree = newicks[0]

            # Build a Tree object from the description stored in the data file
            self.tree.buildFromString(self.starting_tree)
            
        elif self.starting_tree_source == 'random':
            assert self.ntax > 0, 'expecting ntax to be greater than 0'
            
            # Build a random tree
            edge_dist_param = 1.0/self.edgelen_prior_mean
            edge_len_dist = ProbDist.ExponentialDist(edge_dist_param)
            edge_len_dist.setLot(self.r)
            Phylogeny.TreeManip(self.tree).randomTree(
                self.ntax,     # number of tips
                self.r,        # pseudorandom number generator
                edge_len_dist, # distribution from which to draw edge lengths
                False)         # Yule tree if True, edge lengths independent if False
            self.starting_tree = self.tree.makeNewick()
            
        elif self.starting_tree_source == 'usertree':
            # self.tree_topology should already be created
            self.starting_tree = self.tree_topology

            # Build a Tree object from the description stored in self.tree_topology
            self.tree.buildFromString(self.starting_tree)
            
        else:
            # throw exception
            assert False, 'starting_tree_source should equal random, file, or usertree, but instead it was this: %s' % self.starting_tree_source

        # Make sure that names of tips equal the string equivalent of the tip node number plus 1
        # This means that when tree topologies are sampled, the tree definitions that are output
        # to the .t file will be similar to MrBayes output
        assert self.ntax > 0, 'expecting ntax to be greater than 0'
        self.nedges = 2*self.ntax - 3
        taxNames = []
        for i in range(self.ntax):
            taxNames.append(str(i+1))
        self.tree.rectifyNames(taxNames)

    def showParamInfo(self, p):
        print '  Parameter name:    ', p.getName()
        print '  Prior distribution:', p.getPriorDescr()
        if p.isMasterParameter():
            print '  Master parameter (no current value)'
        else:
            print '  Current value:     ', p.getCurrValue()
        print '  Prior log-density: ', p.getLnPrior()
        print
                
    def showTopoPriorInfo(self):
        print 'Topology prior:'
        if not self.allow_polytomies:
            print '  flat across all fully-resolved tree topologies (polytomies not allowed)'
        else:            
            print '  Prior type:',
            if self.topo_prior_calculator.isPolytomyPrior():
                print 'polytomy prior'
            else:
                print 'resolution class prior'
            print '  Prior strength (C):',self.topo_prior_calculator.getC()
            print '  Expected prior probability for each resolution class:'
            print '   class        prior'
            print '  -------------------'
            topo_priors = self.topo_prior_calculator.getRealizedResClassPriorsVect()
            for i,v in enumerate(topo_priors):
                if i == 0:
                    denom = v
                else:
                    print '%8d %12.5f' % (i,math.exp(v - denom))
            print
            #raw_input('stopped after outputting topo prior table')

    def setupLikelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Create a TreeLikelihood object. This can be used to compute the
        likelihood for any tree based on the supplied data matrix and using
        the supplied model.
        
        """
        assert self.model, 'create Phycas.model before calling Phycas.setupLikelihood'
        self.likelihood = Likelihood.TreeLikelihood(self.model)

    def readNexusFile(self, fn):            
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets self.data_file_name to supplied filename fn, calls the
        self.reader.readFile function, storing the data matrix in 
        self.data_matrix. Also sets self.ntax and self.nchar accordingly.
        
        """
        assert not self.no_data, 'set no_data to False before calling readNexusFile'
        if not self.data_file_name == fn:
            self.data_file_name = fn
        self.reader.readFile(self.data_file_name)
        self.data_matrix = ReadNexus.getDiscreteMatrix(self.reader, 0)
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
        # used to avoid next two lines if self.no_data was true, but we have to
        # copy over the data from data_matrix to set self.npatterns, so without the
        # next two lines, we would not be able to run the pattern-specific rates model
        # without data because we would not know how many rate parameters to create
        assert self.likelihood, 'call Phycas.setupLikelihood before calling Phycas.readNexusFile'
        self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
        self.npatterns = self.likelihood.getNPatterns()
                
    def setup(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is for parts of the setup that should occur right before
        run() is called. In fact, calling setup() is the first thing done
        inside the run() function, and thus this function should not be
        called at all - just call run(). Some things, such as setting the
        random number seed, are deferred until this point to give the user of
        the class a chance to change the defaults between class construction
        and the call to run().
        """
        self.nsamples = self.ncycles//self.sample_every

        # Set seed if user has supplied one
        if self.random_seed != 'auto':
            self.r.setSeed(int(self.random_seed))

        # Read the data file and store the number of taxa and number of characters
        if not self.no_data:
            #print '*** In Phycas.setup(), reading data file %s, ntax = %d' % (self.data_file_name, self.ntax)  # TEMPORARY
            self.reader.readFile(self.data_file_name)
            self.data_matrix = ReadNexus.getDiscreteMatrix(self.reader, 0)
            self.ntax = self.data_matrix.getNTax()
            self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only

        self.setupTree()        
        self.setupModel()        
        self.setupLikelihood()

        if not self.no_data:
            self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
            self.npatterns = self.likelihood.getNPatterns()
            
        # Add data structures to the nodes of the tree to allow likelihood calculations
        # The structures added to tips allow the tips to store data and transition probability matrices
        # The structures added to the internal nodes allow for storage of transition probability matrices
        # as well as conditional likelihood arrays
        self.likelihood.prepareForLikelihood(self.tree)

        # Create parameter and tree file names based on the data file name
        prefix = os.path.basename(self.data_file_name)        
        self.param_file_name = prefix + '.p'
        self.tree_file_name = prefix + '.t'

        # Open the parameter file
        self.paramFileOpen()
        
        # Store (and create, if necessary) list of taxon labels
        self.taxon_labels = []
        if self.no_data:
            for i in range(self.ntax):
                s = 'taxon_%d' % (i + 1)
                self.taxon_labels.append(s)
        else:
            labels_in_datafile = self.reader.getTaxLabels()
            for i in range(self.ntax):
                self.taxon_labels.append(labels_in_datafile[i])

        # Open the tree file
        self.treeFileOpen()
        
        if self.verbose:
            if self.no_data:
                print 'Data file:     ', 'None (running MCMC with no data to explore prior)'
            else:
                print 'Data file:     ', self.data_file_name
            print 'Prior:         ', self.edgelen_prior_mean
            print 'No. cycles:    ', self.ncycles
            print 'Sample every:  ', self.sample_every
            print 'Starting tree: ', self.starting_tree
            print 'No. samples:   ', self.nsamples
            print 'Sampled trees will be saved in', self.tree_file_name
            print 'Sampled parameters will be saved in', self.param_file_name

            if self.tree.tipNumbersSetUsingNames():
                print 'Tip node numbers were set using the names in the tree description'
            else:
                print 'Warning: tip node numbers were NOT set using the names in the tree description'

    def treeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the tree file and writes a translate table.
        
        """
        self.treef = file(self.tree_file_name, 'w')
        self.treef.write('#NEXUS\n')
        self.treef.write('[ID: %d]\n' % self.r.getInitSeed())
        self.treef.write('begin trees;\n')
        self.treef.write('   translate\n')
        for i in range(self.ntax):
            if self.taxon_labels[i].find(' ') < 0:
                # no spaces found in name
                self.treef.write('       %d %s%s\n' % (i + 1, self.taxon_labels[i], i == self.ntax - 1 and ';' or ','))
            else:
                # at least one space in taxon name, so enclose name in quotes
                self.treef.write("       %d '%s'%s\n" % (i + 1, self.taxon_labels[i], i == self.ntax - 1 and ';' or ','))

    def treeFileClose(self):
        self.treef.write('end;\n')
        self.treef.close()

    def paramFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the parameter file and writes a header line.
        
        """
        self.paramf = file(self.param_file_name, 'w')
        self.paramf.write('[ID: %d]\n' % self.r.getSeed())
        self.paramf.write(self.model.paramHeader())
        if self.using_hyperprior:
            self.paramf.write('\thyper')
        self.paramf.write('\n')

    def paramFileClose(self):
        self.paramf.close()

    def recordSample(self, cycle, lnL = 0.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Records current tree topology and edge lengths by adding a line to
        the tree file, and records tree length and substitution parameters
        by adding a line to the parameter file. If Gelfand-Ghosh statistics
        were requested, one or more simulation replicates are performed to
        this end.
        
        """
        # Add line to parameter file if it exists
        if self.paramf:
            self.paramf.write('%d\t%.3f\t%.3f' % (cycle + 1, lnL, self.tree.edgeLenSum()))
            self.paramf.write(self.model.paramReport())
            if self.using_hyperprior:
                p = self.chain_manager.getEdgeLenHyperparam()
                self.paramf.write('\t%.5f' % p.getCurrValue())
            self.paramf.write('\n')
        
        # Add line to tree file if it exists
        if self.treef:
            self.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, self.tree.makeNewick()))

        # Perform posterior predictive simulations if Gelfand-Ghosh statistics requested
        if self.gg_do and cycle > 0:
            for j in range(self.gg_nreps):
                # Simulate from the posterior
                sim_data = Likelihood.SimData()
                self.likelihood.simulate(sim_data, self.tree, self.r, self.nchar)

                if self.gg_save_spectra:
                    self.gg_spectrum.zeroCounts()
                    sim_data.addDataTo(self.gg_spectrum, 1.0)
                    self.gg_spectrum_points += ','
                    self.gg_spectrum_points += self.gg_spectrum.createMapleTuples(self.gg_spectrum_row, 100)
                    self.gg_spectrum_row += 1
                    self.gg_spectrum.appendCountsToFile('spectra.txt', False)

                # Save the simulated data set if desired
                if self.gg_save_postpreds:
                    fn = '%s_cycle%d_rep%d.nex' % (self.gg_postpred_prefix, cycle, j)
                    sim_data.saveToNexusFile(fn, self.taxon_labels, 'dna', ('a','c','g','t'))
                    simf = file(fn, 'a')
                    simf.write('\nbegin trees;\n')
                    simf.write('  translate\n')
                    for num,name in enumerate(self.taxon_labels):
                        simf.write("  %d '%s'%s\n" % (num + 1, name, num == self.ntax - 1 and ';' or ','))
                    simf.write('  ;\n')
                    simf.write('  utree one = %s;\n' % self.tree.makeNewick())
                    simf.write('end;\n')
                    simf.write('\nbegin paup;\n')
                    simf.write('  log file=%s.log start replace;\n' % fn)
                    simf.write('  lset nst=6 basefreq=estimate rmatrix=estimate pinvar=0.0 rates=gamma shape=estimate;\n')
                    simf.write('  lscores 1;\n')
                    simf.write('  log stop;\n')
                    simf.write('  quit;\n')
                    simf.write('end;\n')
                    simf.write('[\n')
                    simf.write('\n')
                    simf.write('cycle                  = %d\n' % cycle)
                    simf.write('lnL                    = %f\n' % lnL)
                    simf.write('TL                     = %f\n' % self.tree.edgeLenSum())
                    if self.using_hyperprior:
                        p = self.chain_manager.getEdgeLenHyperparam()
                        simf.write('edge length hyperparam = %f\n' % p.getCurrValue())
                    simf.write('param headers          = %s\n' % self.model.paramHeader())
                    simf.write('param values           = %s\n' % self.model.paramReport())
                    simf.write(']\n')
                    simf.close()

                # Compute the t function for the simulated dataset                
                curr_t = sim_data.calct(4)

                # Add this value of t to the list (later the mean t will be computed)                
                self.gg_t.append(curr_t)

                # Add the number of patterns in sim_data to the gg_npatterns list
                self.gg_npatterns.append(sim_data.getNUniquePatterns())

                # Add the pattern counts for this data set to gg_mu (later the mean counts will be computed)
                sim_data.addDataTo(self.gg_mu, 1.0)

                # Increment count of the total number of simulated datasets created
                # This value is used to later compute the mean t for all simulated datasets
                # and the mean counts for all simulated data sets
                self.gg_total += 1

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
        assert p, 'could not adapt slice sampler; parameter non-existant'
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
        self.timer_stop = time.clock()
        self.elapsed_secs += self.timer_stop - self.timer_start

        summary = ''
        for p in self.chain_manager.getAllUpdaters():
            summary += self.adaptOneSliceSampler(p)
        
        if self.verbose and summary != '':
            print '\nSlice sampler diagnostics:'
            print summary
            
        self.timer_start = time.clock()

    def fillSpectrum(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Only used for 4-taxon problems where there are 256 possible patterns
        for DNA data. Creates a SimData object with 256 patterns. This object
        is used as a workspace for saving pattern spectra if gg_save_spectra
        is True.
        
        """
        self.gg_spectrum_row = 0
        self.gg_spectrum.resetPatternLength(4)
        for i in range(4):        
            for j in range(4):        
                for k in range(4):        
                    for m in range(4):
                        self.gg_spectrum.setState(0, i)
                        self.gg_spectrum.setState(1, j)
                        self.gg_spectrum.setState(2, k)
                        self.gg_spectrum.setState(3, m)
                        self.gg_spectrum.insertPattern(1.0)
            
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls setup() to finish preparations and then runs the Markov chain.
        Delaying the call to setup() until now allows some parameters to be
        changed after the MCMCSimple object is constructed but before the
        analysis commences.
        
        """
        print '*** Just entered Phycas.run(), ntax =',self.ntax
        if self.gg_do:
            assert not self.no_data, 'cannot set gg_do and no_data both to True'
            assert self.nchar > 0, 'nchar not set, required for Gelfand-Ghosh calculation'

            # Clear gg_y and let it contain the observed pattern counts            
            self.gg_y.clear()
            self.likelihood.addDataTo(self.gg_y)

            # If saving spectra, save the spectrum from the original data set
            if self.gg_save_spectra:
                assert self.ntax == 4, 'gg_save_spectra is designed for 4-taxon problems only (i.e. ntax = 4); ntax is %d in this case' % self.ntax
                self.fillSpectrum()
                self.gg_spectrum.zeroCounts()
                self.gg_y.addDataTo(self.gg_spectrum, 1.0)
                yobs_row = self.gg_spectrum.createMapleTuples(self.gg_spectrum_row, 100)
                self.gg_spectrum_points += yobs_row
                self.gg_spectrum_row += 1
                self.gg_spectrum.appendCountsToFile('spectra.txt', False)
                for z in range(9):
                    self.gg_spectrum_points += ','
                    self.gg_spectrum_points += yobs_row
                    self.gg_spectrum_row += 1
                    self.gg_spectrum.appendCountsToFile('spectra.txt', False)

            # Clear the other quantities that will depend on posterior simulations            
            self.gg_t_mean = 0.0
            self.gg_t_mu = 0.0
            self.gg_t = []
            self.gg_a = []
            self.gg_t_a = []
            self.gg_Pm = 0.0
            self.gg_Gm = []
            self.gg_Dm = []
            self.gg_mu.clear()
            self.gg_total = 0

        # Create a single Markov chain and add the parameters needed by the
        # model (as well as Metropolis moves to modify the topology)
        self.createChain()

        # Tell TreeLikelihood object if user wants to run with no data
        if self.no_data:
            self.likelihood.setNoData()

        # Compute the current log-likelihood and log-prior in case first updater 
        # is a move and will thus depend on these quantities being accurate
        self.chain_manager.refreshLastLnLike()
        print 'Starting log-likelihood =',self.chain_manager.getLastLnLike()
        self.chain_manager.refreshLastLnPrior()
        print 'Starting log-prior =',self.chain_manager.getLastLnPrior()

        # Show starting parameter info 
        print 'Parameter starting values and prior densities:'
        for p in self.chain_manager.getEdgeLenParams():
            self.showParamInfo(p)
        if self.using_hyperprior:
            p = self.chain_manager.getEdgeLenHyperparam()
            self.showParamInfo(p)
        for p in self.chain_manager.getModelParams():
            self.showParamInfo(p)

        # Show information about topology prior to be used
        self.showTopoPriorInfo()
            
        self.elapsed_secs = 0.0
        self.likelihood.resetNEvals()

        print '\nSampling (%d cycles)...' % self.ncycles
        if self.verbose:
            print
        self.recordSample(0)
        self.timer_start = time.clock()
        last_adaptation = 0
        next_adaptation = self.adapt_first
        for cycle in range(self.ncycles):
            for p in self.chain_manager.getAllUpdaters():
                w = p.getWeight()
                #print '*** Updating %s...' % p.getName()
                for x in range(w):
                    p.update()
            if self.verbose and (cycle + 1) % self.report_every == 0:
                print 'cycle = %d, lnL = %.5f' % (cycle + 1, self.chain_manager.getLastLnLike())
            if (cycle + 1) % self.sample_every == 0:
                self.recordSample(cycle, self.chain_manager.getLastLnLike())
            if (cycle + 1) % next_adaptation == 0:
                self.adaptSliceSamplers()
                next_adaptation += 2*(next_adaptation - last_adaptation)
                last_adaptation = cycle + 1
        self.adaptSliceSamplers()
        total_evals = self.likelihood.getNEvals()
        print '%d likelihood evaluations in %.5f seconds' % (total_evals, self.elapsed_secs)
        if (self.elapsed_secs > 0.0):
            print '  = %.5f likelihood evaluations/sec' % (total_evals/self.elapsed_secs)

        if self.treef:
            self.treeFileClose()
        if self.paramf:
            self.paramFileClose()

        if self.gg_do:
            two_n = 2.0*float(self.nchar)

            # Compute the t function for the observed dataset
            self.gg_t_y = self.gg_y.calct(4)

            # Divide the counts stored in gg_mu by gg_total so that gg_mu becomes
            # the mean of all the posterior predictive datasets. Also compute the
            # t function for the mean dataset
            self.gg_mu.divideBy(float(self.gg_total))
            self.gg_t_mu = self.gg_mu.calct(4)

            # If saving spectra, save the spectrum from the original data set
            if self.gg_save_spectra:
                self.gg_spectrum.zeroCounts()
                self.gg_mu.addDataTo(self.gg_spectrum, 1.0)
                self.gg_spectrum_points += ','
                self.gg_spectrum_points += self.gg_spectrum.createMapleTuples(self.gg_spectrum_row, 100)
                self.gg_spectrum_row += 1
                self.gg_spectrum.appendCountsToFile('spectra.txt', False)

                # get patterns
                patterns = self.gg_spectrum.getPatterns(['A','C','G','T'])
                spectf = file('spectra.txt', 'a')
                for i in range(4):
                    for j,p in enumerate(patterns):
                        entry = '%s%s' % ((j == 0 and '' or '\t'), p[i])
                        #if j == 0:
                        #    raw_input('check entry: |%s|' % entry)
                        spectf.write(entry)
                    spectf.write('\n')
                spectf.close()

                # write out the maple_commands file now
                maplef = file('maple_commands', 'w')
                maplef.write('with(linalg);\n')
                maplef.write('with(plots);\n')
                maplef.write('points := [')
                maplef.write(self.gg_spectrum_points)
                maplef.write('];\n')
                maplef.write('surfdata(points, style=patchnogrid, axes=framed, labels=["sim", "pattern", "freq"]);\n')
                maplef.close()

            # Compute the mean of the t values computed for individual posterior
            # predictive datasets
            self.gg_t_mean = sum(self.gg_t)/float(self.gg_total)

            # Compute the penalty term. Guaranteed to be positive by Jensen's
            # inequality and the convexity of the t function.
            self.gg_Pm = two_n*(self.gg_t_mean - self.gg_t_mu)

            # Loop over k values, computing Gm and Dm for each k value in gg_kvect
            for k in self.gg_kvect:
                # Create a dataset representing the compromise "action"
                a = Likelihood.SimData()
                self.gg_mu.addDataTo(a, 1.0)
                self.gg_y.addDataTo(a, k)
                a.divideBy(k + 1.0)
                t_a = a.calct(4)
                self.gg_t_a.append(t_a)
                self.gg_a.append(a)

                # Compute the goodness-of-fit term
                Gkm = (float(k) + 1.0)*two_n*((self.gg_t_mu + k*self.gg_t_y)/(k + 1.0) - t_a)
                self.gg_Gm.append(Gkm)

                # Compute the overall measure            
                Dkm = self.gg_Pm + Gkm
                self.gg_Dm.append(Dkm)

            if self.gg_outfile:
                ggf = file(self.gg_outfile,'w')
                ggf.write('# Pm = %f\n' % self.gg_Pm)
                for i,k in enumerate(self.gg_kvect):
                    ggf.write('# k = %f:\n' % k)
                    ggf.write('#   Gm = %f\n' % self.gg_Gm[i])
                    ggf.write('#   Dm = %f\n' % self.gg_Dm[i])
                ggf.write('\n')

                ggf.write('# no. patterns in original dataset                        = %d\n' % self.gg_y.getNUniquePatterns())
                ggf.write('# no. patterns in mean over posterior preditive datasets  = %d\n' % self.gg_mu.getNUniquePatterns())
                sum_npat = 0.0
                for npat in self.gg_npatterns:
                    sum_npat += float(npat)
                ggf.write('# mean no. patterns over posterior preditive datasets     = %f\n' % (sum_npat/float(len(self.gg_npatterns))))

                ggf.write('# t calculated for original dataset                       = %f\n' % self.gg_t_y)
                ggf.write('# t calculated for mean over posterior preditive datasets = %f\n' % self.gg_t_mu)
                ggf.write('# mean of t over posterior preditive datasets             = %f\n' % self.gg_t_mean)
                ttotal = len(self.gg_t)
                assert ttotal == self.gg_total, 'mismatch between self.gg_total and len(self.gg_t)'
                tsumsq = 0.0
                for t in self.gg_t:
                    tsumsq += t*t
                tvar = tsumsq - float(ttotal)*self.gg_t_mean*self.gg_t_mean
                ggf.write('# std. dev. of t over posterior preditive datasets        = %f\n' % math.sqrt(tvar))
                for i,k in enumerate(self.gg_kvect):
                    ggf.write('# t of compromise action for k = %6f                   = %f\n' % (k, self.gg_t_a[i]))

                ggf.write('\n# GnuPlot commands for making a plot of t values:\n')
                ggf.write('set title "%s"\n' % (self.gg_outfile))                                                           
                ggf.write('\n')
                arrow_number = 1
                ggf.write('# t_y\n')
                ggf.write('set arrow %d from %f,0.95 to %f,1.0 lw 2 lt %d\n' % (arrow_number, self.gg_t_y, self.gg_t_y, arrow_number))
                ggf.write('set label %d "t_y" at %f,0.925\n' % (arrow_number, self.gg_t_y))
                ggf.write('\n')
                arrow_number = 2
                ggf.write('# t_mean\n')
                ggf.write('set arrow %d from %f,0.95 to %f,1.0 lw 1 lt %d\n' % (arrow_number, self.gg_t_mean, self.gg_t_mean, arrow_number))
                ggf.write('set label %d "t_mean" at %f,0.9125\n' % (arrow_number, self.gg_t_mean))
                ggf.write('\n')
                arrow_number = 3
                ggf.write('# t_mu\n')
                ggf.write('set arrow %d from %f,0.95 to %f,1.0 lw 2 lt %d\n' % (arrow_number, self.gg_t_mu, self.gg_t_mu, arrow_number))
                ggf.write('set label %d "t_mu" at %f,0.925\n' % (arrow_number, self.gg_t_mu))
                ggf.write('\n')
                tmin = min(self.gg_t)
                tmax = max(self.gg_t)
                xtremes = [tmin, tmax, self.gg_t_y, self.gg_t_mu, self.gg_t_mean] 
                for i,k in enumerate(self.gg_kvect):
                    arrow_number = 4 + i
                    ggf.write('# t_a (k=%f)\n' % (k))
                    ggf.write('#set arrow %d from %f,0.95 to %f,1.0 lw 2 lt %d\n' % (arrow_number, self.gg_t_a[i], self.gg_t_a[i], arrow_number))
                    ggf.write('#set label %d "t_a" at %f,0.925\n' % (arrow_number, self.gg_t_a[i]))
                    xtremes.append(self.gg_t_a[i])
                ggf.write('\n')
                ggf.write('set nokey\n')
                ggf.write('plot [%f:%f][0.8:1.2] "%s" using 2:3 with points\n' % (min(xtremes), max(xtremes),self.gg_outfile))
                ggf.write('pause -1 "Press return to continue..."\n')
                ggf.write('\n# GnuPlot data:\n')
                for i,t in enumerate(self.gg_t):
                    ggf.write('%d\t%f\t1.0\n' % (i, t))
                ggf.close()
                
                #ggf = file(self.gg_outfile+'.patterns.txt','w')
                #ggf.write('*** yobs:\n')
                #ggf.write(self.gg_y.patternTable(['A', 'C', 'G', 'T']))
                #ggf.write('\n*** mu:\n')
                #ggf.write(self.gg_mu.patternTable(['A', 'C', 'G', 'T']))
                #for i,k in enumerate(self.gg_kvect):
                #    ggf.write('\n*** a (k = %f):\n' % k)
                #    ggf.write(self.gg_a[i].patternTable(['A', 'C', 'G', 'T']))
                #ggf.close()
        
if __name__ == '__main__':
    mcmc = Phycas()
    
    mcmc.data_file_name = '../../pyphy/nyldna4.nex'
    mcmc.starting_tree_source = 'usertree'
    mcmc.tree_topology = '(0:0.1,1:0.1,(2:0.1,3:0.1):0.1)'
    mcmc.ncycles = 1500
    mcmc.sample_every = 15
    mcmc.adapt_first = 60
    mcmc.random_seed = '13579'
    mcmc.model_type = 'hky'
    mcmc.using_hyperprior = True
    mcmc.edgelen_prior_mean = 0.1
    mcmc.verbose = True
    mcmc.metropolis_weight = 300
    mcmc.slice_weight = 1
    mcmc.gg_do = False
    mcmc.slice_max_units = 0
    mcmc.use_inverse_shape = False

    mcmc.setup()
    mcmc.run()
    
