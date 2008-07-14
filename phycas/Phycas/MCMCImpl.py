import os,sys,math,random
from phycas import *
from phycas.Phycas.PhycasCommand import *
from MCMCManager import MCMCManager
from phycas.ProbDist import StopWatch
from phycas.ReadNexus import NexusReader

def check(msg = 'check'):
    raw_input(msg)

class MCMCImpl(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Needs to be written.
    
    """
    def __init__(self, phycas, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes MCMCRunner object by assigning supplied phycas object
        to a data member variable.
        
        """
        self.phycas = phycas
        self.opts = opts
        
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
        self.tree_file_name         = ''        # Will hold tree file name (see openParameterAndTreeFiles)
        self.param_file_name        = ''        # Will hold parameter file name (see openParameterAndTreeFiles)
        #self.tmp_simdata            = SimData()
        self.gg_Pm                  = 0.0       # Penalty component (same for all k)
        self.gg_Gm                  = []        # Vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []        # Vector of overall measures (one for each k in gg_kvect)
        self.reader                 = NexusReader()
        self.logf                   = None
        self._logFileName           = None
        self.addition_sequence      = []        # List of taxon numbers for addition sequence
        self.samc_theta             = []        # Normalizing factors (will have length ntax - 3 because levels with 1, 2 or 3 taxa are not examined)
        self.samc_distance_matrix   = None      # Holds ntax x ntax hamming distance matrix used by SamcMove
        self.stored_tree_defs       = None
        self.ps_delta_beta          = 0.0
        self.doing_path_sampling    = False
        self.path_sample            = None
        self.psf                    = None
        self.pdf_splits_to_plot     = None
        self.param_file_name        = None  
        self.tree_file_name         = None
        self.nsamples               = None
        self.ps_beta                = 1.0
        self.wangang_sampled_betas  = None
        self.wangang_sampled_likes  = None
        self.unimap_manager         = None
        self.nsamples               = 0

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
        self.phycas.phycassert(p, 'could not adapt slice sampler; parameter non-existant')
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
            self.phycas.output('\nSlice sampler diagnostics:')
            self.phycas.output(summary)

    def updateAllUpdaters(self, chain, chain_index, cycle):
        if self.opts.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
        for p in chain.chain_manager.getAllUpdaters():
            w = p.getWeight()
            #print p.getName(), "weight =", w
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
        self.phycas.output('Topology prior:')
        if not self.opts.allow_polytomies:
            self.phycas.output('  flat across all fully-resolved tree topologies (polytomies not allowed)')
        else:            
            if m.topo_prior_calculator.isPolytomyPrior():
                self.phycas.output('  Prior type: polytomy prior')
            else:
                self.phycas.output('  Prior type: resolution class prior')
            self.phycas.output('  Prior strength (C): %s' % m.topo_prior_calculator.getC())
            self.phycas.output('  Prior probability for each resolution class:')
            self.phycas.output('  Note: 0.00000000 does *not* mean that the prior is zero! It simply')
            self.phycas.output('        indicates that the prior is less than 0.000000005\n')
            self.phycas.output('%20s %20s' % ('internal nodes', 'prior probability'))
            self.phycas.output('%20s %20s' % ('--------------', '-----------------'))
            topo_priors = m.topo_prior_calculator.getRealizedResClassPriorsVect()
            for i,v in enumerate(topo_priors):
                if i == 0:
                    denom = v   # first element of vector is log of normalization constant (sum of all other elements)
                else:
                    topo_prior = math.exp(v - denom)
                    self.phycas.output('%20d %20.8f' % (i,topo_prior))
            self.phycas.output()

    def showParamInfo(self, p):
        self.phycas.output('  Parameter name:     %s' % p.getName())
        self.phycas.output('  Prior distribution: %s' % p.getPriorDescr())
        if p.isMasterParameter():
            self.phycas.output('  Master parameter (no current value)')
        else:
            self.phycas.output('  Current value:      %s' % p.getCurrValue())
        self.phycas.output('  Prior log-density:  %s' % p.getLnPrior())
        self.phycas.output()
                
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

    def getPrefix(self):
        prefix = os.path.abspath(self.opts.data_file_name) #os.path.basename(self.data_file_name)
        if self.opts.outfile_prefix:
            prefix = self.opts.outfile_prefix
        return prefix
    
    def openParameterAndTreeFiles(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates parameter and tree file names based on the data file name or the
        user-supplied prefix and opens the files
        
        """
        prefix = self.getPrefix()
        self.param_file_name = prefix + '.p'
        self.tree_file_name = prefix + '.t'

        self.paramFileOpen()
        self.treeFileOpen()

    def readDataFromFile(self):
        self.reader.readFile(self.opts.data_file_name)
        self.taxon_labels = self.reader.getTaxLabels()
        self.data_matrix = self.reader.getLastDiscreteMatrix(True)
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only

    def calcMarginalLikelihood(self):
        marginal_like = 0.0
        if self.doing_path_sampling:
            # Calculate marginal likelihood using continuous path sampling
            # The path_sampling vector is 1-dimensional, each element holds average for one increment
            C = len(self.path_sample) - 1
            marginal_like += (self.path_sample[0] + self.path_sample[-1])/2.0
            for v in self.path_sample[1:-1]:
                marginal_like += v
            marginal_like /= float(C)
            if self.ps_filename:
                self.psf.write('%s\t%.5f\n' % ('-->',marginal_like))
                self.psf.close()
            self.phycas.output('Marginal likelihood (continuous path sampling method) = %f' % marginal_like)
        elif self.opts.nchains > 1 and not self.opts.is_standard_heating:
            # Calculate marginal likelihood using discrete path sampling
            # The path_sampling vector is 2-dimensional, each element holds samples from one chain
            C = self.opts.nchains - 1
            self.phycas.output('\nCalculation of marginal likelihood:')
            self.phycas.output('%12s%12s' % ('chain', 'avg. lnL'))
            for i,v in enumerate(self.path_sample):
                n = len(v)
                avg = sum(v)/float(n)
                self.phycas.output('%12d%12.5f' % (i, avg))
                if (i == 0) or (i == C):
                    avg /= 2.0
                marginal_like += avg
            marginal_like /= float(C)
            self.phycas.output('  Marginal likelihood (discrete path sampling method) = %f' % marginal_like)
            
            # Calculate marginal likelihood using harmonic mean method on cold chain
            sample_size = len(self.path_sample[0])
            min_lnL = min(self.path_sample[0])
            sum_diffs = 0.0
            for lnl in self.path_sample[0]:
                diff = lnl - min_lnL
                if diff < 500.0:
                    sum_diffs += math.exp(-diff)
                else:
                    self.phycas.output('warning: ignoring large diff (%f) in harmonic mean calculation' % diff)
            log_harmonic_mean = math.log(sample_size) + min_lnL - math.log(sum_diffs)
            self.phycas.output('  Marginal likelihood(harmonic mean method)= %f' % log_harmonic_mean)

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
        # REVISIT LATER
        # Hack until PhycasCommand implements properties
        if self.opts.internal_edgelen_dist is None:
            self.opts.internal_edgelen_dist = self.opts.edgelen_dist
        if self.opts.external_edgelen_dist is None:
            self.opts.external_edgelen_dist = self.opts.edgelen_dist
        
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
            if self.opts.nchains == 1:
                self.heat_vector = [1.0]
            else:
                # DISCRETE PATH SAMPLING
                # Create a list for each chain to hold sampled lnL values
                #@ shouldn't this only be done if self.is_standard_heating is False?
                self.path_sample = []
                for i in range(self.opts.nchains):
                    self.path_sample.append([])
                
                # Determine vector of powers for each chain
                self.heat_vector = []
                if self.is_standard_heating:
                    # Goal is to improve mixing; both likelihood and prior raised to power
                    for i in range(self.opts.nchains):
                        # Standard heating 
                        # 0 1.000 = 1/1.0 cold chain explores posterior
                        # 1 0.833 = 1/1.2
                        # 2 0.714 = 1/1.4
                        # 3 0.625 = 1/1.6
                        temp = 1.0/(1.0 + float(i)*self.heating_lambda)
                        self.heat_vector.append(temp)
                else:
                    # DISCRETE PATH SAMPLING
                    # Goal is path sampling, not mixing, so these powers are applied only to likelihood
                    for i in range(self.opts.nchains):
                        self.do_marginal_like = True
                        # Likelihood heating for thermodynamic integration
                        # 0 1.000 = (3-0)/3 cold chain explores posterior
                        # 1 0.667 = (3-1)/3
                        # 2 0.333 = (3-2)/3
                        # 3 0.000 = (3-3)/3 hottest chain explores prior
                        temp = 1.0
                        denom = float(self.opts.nchains - 1)
                        temp = float(self.opts.nchains - i - 1)/denom
                        self.heat_vector.append(temp)
        else:
            # User supplied his/her own heat_vector; perform sanity checks
            self.opts.nchains = len(self.heat_vector)
            self.phycas.phycassert(self.heat_vector.index(1.0) < self.opts.nchains, 'user-supplied heat_vector does not allow for a cold chain (one power must be 1.0)')

        self.mcmc_manager.createChains()
        self.openParameterAndTreeFiles()
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the MCMC analysis. 
        
        """
        self.setup()
        
        # If user has set quiet to True, then phycas.output calls will have no effect
        self.phycas.quiet = self.opts.quiet
        
        # Tell TreeLikelihood object if user wants to run with no data
        #if not self.data_source:
        #    self.likelihood.setNoData()

        self.nsamples = self.opts.ncycles//self.opts.sample_every
        
        if self.opts.verbose:
            if self.opts.data_source == None:
                self.phycas.output('Data source:    None (running MCMC with no data to explore prior)')
            elif self.opts.data_source == 'file':
                self.phycas.output('Data source:    %s' % self.opts.data_file_name)
            else:
                self.phycas.abort("Only 'file' or None are allowed for data_source")
            self.phycas.output('No. cycles:     %s' % self.opts.ncycles)
            self.phycas.output('Sample every:   %s' % self.opts.sample_every)
            self.phycas.output('Starting tree:  %s' % self.starting_tree)
            self.phycas.output('No. samples:    %s' % self.nsamples)
            self.phycas.output('Sampled trees will be saved in %s' % self.tree_file_name)
            self.phycas.output('Sampled parameters will be saved in %s' % self.param_file_name)
            if self.opts.use_unimap:
                self.phycas.output('Using uniformized mapping MCMC')
            else:
                self.phycas.output('Using standard MCMC (i.e. no uniformized mapping)')

            if not self.warn_tip_numbers:
                self.phycas.output('Tip node numbers were set using the names in the tree description')
            else:
                self.phycas.output('Warning: tip node numbers were NOT set using the names in the tree description')

        if self.opts.nchains == 1:
            self.phycas.output('Creating one chain (i.e. not using heated chains to improve mixing)')
        else:
            self.phycas.output('Creating %d chains with these temperatures:' % (self.opts.nchains))
            for t in self.heat_vector:
                self.phycas.output('  %.5f %s' % (t, t == 1.0 and '(cold chain)' or ''))
            
        # Compute the current log-likelihood and log-prior in case first updater 
        # is a move and will thus depend on these quantities being accurate
        for c in self.mcmc_manager.chains:
            c.chain_manager.refreshLastLnLike()
            c.chain_manager.refreshLastLnPrior()
            if c.heating_power == 1.0:
                self.phycas.output('Starting log-likelihood = %s' % c.chain_manager.getLastLnLike())
                self.phycas.output('Starting log-prior = %s' % c.chain_manager.getLastLnPrior())

        # Show starting parameter info 
        self.phycas.output('\nParameter starting values and prior densities:')
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
        
        self.phycas.output('\nSampling (%d cycles)...' % self.opts.ncycles)
        if self.opts.verbose:
            print
        self.mcmc_manager.recordSample(0)
        last_adaptation = 0
        next_adaptation = self.opts.adapt_first

        # REVISIT LATER
        #if self.doing_path_sampling:
        #    self.path_sample = []
        #    chain = self.mcmc_manager.chains[0]
        #    ps_Qsum   = 0.0
        #    ps_Qtotal = 0.0
        #    ps_Qnum   = 0
        #    if self.ps_toward_posterior:
        #        self.ps_beta = self.ps_minbeta
        #    else:
        #        self.ps_beta = self.ps_maxbeta
        #    chain.setPower(self.ps_beta)
        #    if self.ps_filename:
        #        self.psf = open(self.ps_filename,'w')
        #        self.psf.write('beta\tavglnL\n')
        #    self.wangang_sampled_betas = [self.ps_beta]
        #    self.wangang_sampled_likes = []
        #    self.wangang_sampled_likes.append([])
        #    beta_index = 0
            
        for cycle in xrange(self.opts.ncycles):
            for i,c in enumerate(self.mcmc_manager.chains):
                self.updateAllUpdaters(c, i, cycle)
            if self.opts.verbose and (cycle + 1) % self.opts.report_every == 0:
                self.stopwatch.normalize()
                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                msg = 'cycle = %d, lnL = %.5f (%.5f secs)' % (cycle + 1, cold_chain_manager.getLastLnLike(), self.stopwatch.elapsedSeconds())
                self.phycas.output(msg)
            # REVISIT LATER
            #if self.doing_path_sampling and cycle + 1 > self.ps_burnin:
            #    sampled_lnL = cold_chain_manager.getLastLnLike()
            #    ps_Qnum += 1
            #    if (ps_Qnum % self.ps_sample_every) == 0:
            #        ps_Qsum += sampled_lnL
            #        ps_Qtotal += 1.0
            #        self.wangang_sampled_likes[beta_index].append(sampled_lnL)
            #    if ps_Qnum == self.ps_Q:
            #        avg = ps_Qsum/ps_Qtotal
            #        self.path_sample.append(avg)
            #        if self.ps_filename:
            #            self.psf.write('%.3f\t%.5f\n' % (self.ps_beta,avg))
            #        if self.ps_toward_posterior:
            #            self.ps_beta += self.ps_delta_beta
            #        else:
            #            self.ps_beta -= self.ps_delta_beta
            #        chain.setPower(self.ps_beta)
            #        ps_Qsum   = 0.0
            #        ps_Qtotal = 0.0
            #        ps_Qnum   = 0
            #        self.wangang_sampled_betas.append(self.ps_beta)
            #        self.wangang_sampled_likes.append([])
            #        beta_index += 1
            if (cycle + 1) % self.opts.sample_every == 0:
                self.mcmc_manager.recordSample(cycle)
                self.stopwatch.normalize()
            if (cycle + 1) % next_adaptation == 0:
                self.adaptSliceSamplers()
                next_adaptation += 2*(next_adaptation - last_adaptation)
                last_adaptation = cycle + 1

        self.adaptSliceSamplers()
        total_evals = self.mcmc_manager.getTotalEvals() #self.likelihood.getNEvals()
        total_secs = self.stopwatch.elapsedSeconds()
        self.phycas.output('%d likelihood evaluations in %.5f seconds' % (total_evals, total_secs))
        if (total_secs > 0.0):
            self.phycas.output('  = %.5f likelihood evaluations/sec' % (total_evals/total_secs))

        if self.treef:
            self.treeFileClose()
        if self.paramf:
            self.paramFileClose()
            
        # Delete the last (extra) entry in both self.wangang_sampled_betas and self.wangang_sampled_likes
        if self.doing_path_sampling:
            self.wangang_sampled_betas = self.wangang_sampled_betas[:-1]
            self.wangang_sampled_likes = self.wangang_sampled_likes[:-1]

        # If we have been path sampling, compute marginal likelihood using lnL values
        # for each chain stored in self.path_sample
        self.calcMarginalLikelihood()