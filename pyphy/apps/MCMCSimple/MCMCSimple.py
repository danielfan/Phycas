import os, sys, time
import DataMatrix
import ReadNexus
import Phylogeny
import Likelihood
import ProbDist

class MCMCSimple:
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Performs the simplest possible phylogenetic MCMC analysis. The tree
    topology remains fixed and slice sampling is used to update all
    parameters, including all edge lengths. A cycle is defined as one
    round of updating in which each slice sampler is given the chance to
    update its particular parameter.

    The following example reads a nexus data file, specifies a tree
    topology and runs a single Markov chain sampler for 100 cycles using
    the HKY model.
    
    >>> from MCMCSimple import *
    >>> mcmc = MCMCSimple()
    >>> mcmc.data_file_name = 'nyldna4.nex'
    >>> mcmc.tree_topology = '(0:0.1,1:0.1,(2:0.1,3:0.1):0.1)'
    >>> mcmc.ncycles = 100
    >>> mcmc.sample_every = 100
    >>> mcmc.report_every = 10
    >>> mcmc.adapt_first = 25
    >>> mcmc.random_seed = '13579'
    >>> mcmc.model_type = 'hky'
    >>> mcmc.verbose = True
    >>> mcmc.run() # doctest:+ELLIPSIS
    Data file:      nyldna4.nex
    Prior:          0.5
    No. cycles:     100
    Sample every:   100
    Tree topology:  (0:0.1,1:0.1,(2:0.1,3:0.1):0.1)
    No. samples:    1
    Sampled trees will be saved in nyldna4.nex.t
    Sampled parameters will be saved in nyldna4.nex.p
    Tip node numbers were set using the names in the tree description
    Starting log-likelihood = -7574.46605305
    Starting log-prior = 0.621406941325
    Parameter starting values and prior densities:
      Parameter name                    Value   Log density  Prior Distribution
      edge length for node 5          0.10000       2.46574  ExponentialDist(2.00000)
      edge length for node 1          0.10000       2.46574  ExponentialDist(2.00000)
      edge length for node 4          0.10000       2.46574  ExponentialDist(2.00000)
      edge length for node 2          0.10000       2.46574  ExponentialDist(2.00000)
      edge length for node 3          0.10000       2.46574  ExponentialDist(2.00000)
      edge length hyperprior          0.10000      -3.70727  InverseGammaDist(2.10000, 0.90909)
      trs/trv rate ratio              4.00000      -4.00000  ExponentialDist(1.00000)
      base freq. A                    1.00000      -1.00000  ExponentialDist(1.00000)
      base freq. C                    1.00000      -1.00000  ExponentialDist(1.00000)
      base freq. G                    1.00000      -1.00000  ExponentialDist(1.00000)
      base freq. T                    1.00000      -1.00000  ExponentialDist(1.00000)
    <BLANKLINE>
    Sampling (100 cycles)...
    <BLANKLINE>
    cycle = 10, lnL = -7128.28380
    cycle = 20, lnL = -7127.99376
    <BLANKLINE>
    Slice sampler diagnostics:
      mode=0.04763, avgevals=10.880 (edge length for node 5)
      mode=0.06883, avgevals=10.280 (edge length for node 1)
      mode=0.00916, avgevals=11.840 (edge length for node 4)
      mode=0.03791, avgevals=9.080 (edge length for node 2)
      mode=0.06370, avgevals=9.880 (edge length for node 3)
      mode=0.07340, avgevals=7.560 (edge length hyperprior)
      mode=1.79212, avgevals=7.800 (trs/trv rate ratio)
      mode=9.74287, avgevals=7.120 (base freq. A)
      mode=6.08929, avgevals=6.920 (base freq. C)
      mode=6.93802, avgevals=6.960 (base freq. G)
      mode=10.56806, avgevals=7.480 (base freq. T)
    <BLANKLINE>
    cycle = 30, lnL = -7126.93853
    cycle = 40, lnL = -7134.88050
    cycle = 50, lnL = -7131.68925
    cycle = 60, lnL = -7126.69735
    cycle = 70, lnL = -7129.42156
    <BLANKLINE>
    Slice sampler diagnostics:
      mode=0.05399, avgevals=6.180 (edge length for node 5)
      mode=0.06311, avgevals=6.060 (edge length for node 1)
      mode=0.00806, avgevals=6.740 (edge length for node 4)
      mode=0.04009, avgevals=7.940 (edge length for node 2)
      mode=0.06927, avgevals=7.500 (edge length for node 3)
      mode=0.07340, avgevals=6.440 (edge length hyperprior)
      mode=1.88792, avgevals=6.300 (trs/trv rate ratio)
      mode=5.51514, avgevals=5.620 (base freq. A)
      mode=3.43802, avgevals=5.660 (base freq. C)
      mode=3.59260, avgevals=5.760 (base freq. G)
      mode=5.80438, avgevals=6.440 (base freq. T)
    <BLANKLINE>
    cycle = 80, lnL = -7127.18869
    cycle = 90, lnL = -7127.06131
    cycle = 100, lnL = -7134.28293
    <BLANKLINE>
    Slice sampler diagnostics:
      mode=0.05399, avgevals=5.680 (edge length for node 5)
      mode=0.06311, avgevals=6.120 (edge length for node 1)
      mode=0.00806, avgevals=5.880 (edge length for node 4)
      mode=0.04009, avgevals=6.120 (edge length for node 2)
      mode=0.06927, avgevals=5.760 (edge length for node 3)
      mode=0.07340, avgevals=6.120 (edge length hyperprior)
      mode=1.88792, avgevals=6.400 (trs/trv rate ratio)
      mode=3.16380, avgevals=5.400 (base freq. A)
      mode=1.98346, avgevals=6.520 (base freq. C)
      mode=2.36322, avgevals=5.960 (base freq. G)
      mode=3.77085, avgevals=6.400 (base freq. T)
    <BLANKLINE>
    ...

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the object with some default values. If these defaults do
        not suit, they can be changed before run() is called.

        tree_topology is equal to 'file' by default, which means that the
        program will expect to find at least one tree topology defined in the
        nexus data file. tree_topology can alternatively be set to a newick
        tree description if no tree topology is provided in the data file.
        
        prior_param is 'hyper' by default, which means that the mean of the
        edge length prior distribution is governed by a hyperparameter. If
        some number (>0) is specified for prior_param, this number is used
        as the hazard parameter for an exponential distribution that serves
        as the edge length prior distribution.

        The JC model is used by default; set model_type to 'hky' to use the
        HKY model instead, to 'gtr' to use the GTR model.

        random_seed is 'auto' by default, which means that the pseudorandom
        number generator is initialized using the system clock. You can
        however set random_seed to a particular random number seed to repeat
        a previous analysis exactly.

        Other quantities that typically need to be specified include ncycles,
        sample_every, report_every, adapt_first and data_file_name.
        
        """
        self.tree_topology          = 'file'        
        self.using_hyperprior       = True
        self.edgelen_prior_mean     = 0.5
        self.model_type             = 'jc'
        self.verbose                = True
        self.random_seed            = 'auto'
        self.data_file_name         = ''

        # MCMC settings
        self.ncycles                = 10000
        self.sample_every           = 100
        self.report_every           = self.ncycles//100

        # Adaptation of slice samplers is performed the first time at cycle adapt_first.
        # Subsequent adaptations wait twice the number of cycles as the previous adaptation.
        # Thus, adaptation n occurs at cycle adapt_first*(2^n - 1)
        # The total number of adaptations that will occur during an MCMC run is
        #      [ln(adapt_first + ncycles) - ln(adapt_first)]/ln(2)
        self.adapt_first            = 100

        # Experiment varying adapt_first, keeping adapt_simple_param = 0.5
        # adapt_first    secs       evals   evals/sec    no. adapations
        # -------------------------------------------------------------
        #       1        30.52022   658350  21570.94363       13.29
        #       1 (best) 30.49364   658350  21589.74683       13.29
        #      50        30.57674   661482  21633.50307        7.65
        #      50        30.57056   661482  21637.87608        7.65
        #     125        30.68982   669347  21810.06437        6.34
        #     125        30.69365   669347  21807.34360        6.34
        #     250        30.54196   665038  21774.56734        5.36
        #     250        30.56474   665038  21758.33690        5.36
        #     500        30.81820   677774  21992.65503        4.39
        #     500        30.80202   677774  22004.20396        4.39
        #    1000        30.96054   693046  22384.81389        3.46
        #    1000        30.99763   693046  22358.03118        3.46
        #    2000        31.52377   725290  23007.71467        2.58
        #    2000        31.47046   725290  23046.69260        2.58
        #    4000        33.15166   798984  24100.87779        1.81
        #    4000        33.20341   798984  24063.30939        1.81
        #    8000        35.59382   925273  25995.32634        1.17
        #    8000        35.45414   908367  25620.89939        1.17
        #   16000        36.78948   968326  26320.72926        0.70
        #   16000        36.76348   968326  26339.34726        0.70

        # Experiment varying adapt_simple_param, keeping adapt_every = 2500
        # adapt_simple_param    secs       evals   evals/sec
        # ----------------------------------------------------
        #         0.25          33.46301   786597  23506.46617
        #         0.25          33.81276   786597  23263.32027 
        #         0.50          31.05870   739171  23799.15908
        #         0.50          30.97201   739171  23865.77161
        #         0.67          30.90181   744284  24085.45295
        #         0.67 (best)   30.85918   744284  24118.72416
        #         0.75          31.16250   753193  24169.85160
        #         0.75          31.12405   753193  24199.70955
        #         1.00          31.74238   775730  24438.30412
        #         1.00          31.75773   775730  24426.49109
        #

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
        MCMCChainManager object.
        
        """
        # Create a substitution model and define priors for the model parameters
        if self.model_type == 'gtr':
            self.model = Likelihood.GTRModel()
            self.model.setRelRates([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
            self.model.setRelRatePrior(ProbDist.ExponentialDist(1.0))
            self.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
        elif self.model_type == 'hky':
            self.model = Likelihood.HKYModel()
            self.model.setKappa(4.0)
            self.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
            self.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
        else:
            self.model = Likelihood.JCModel()
        
        # Define an edge length prior distribution
        assert self.edgelen_prior_mean > 0.0, 'edgelen_prior_mean must be a positive, non-zero number'
        v = 1.0/float(self.edgelen_prior_mean)
        if self.using_hyperprior:
            # Edge length prior governed by a InverseGamma-distributed hyperprior
            d = ProbDist.InverseGammaDist(2.1, 0.909)
            d.setMeanAndVariance(1.0, 10.0)
            self.model.setEdgeLenHyperPrior(d)
            self.model.setEdgeLenPrior(ProbDist.ExponentialDist(v))
        else:
            # Edge length prior distribution is not hierarchical
            self.model.setEdgeLenHyperPrior(None)
            self.model.setEdgeLenPrior(ProbDist.ExponentialDist(v))

        # Create a TreeLikelihood object. This can be used to compute the likelihood
        # for any tree based on the supplied data matrix and using the supplied model.
        # self.likelihood = Likelihood.TreeLikelihood(self.model, self.data_matrix) #oldway
        self.likelihood = Likelihood.TreeLikelihood(self.model) #newway
        self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix) #newway

        # Create a list of parameters for updating the edge lengths of the tree
        # (and to update kappa and the base frequencies if self.model_type is 'hky')
        self.chain_manager = Likelihood.MCMCChainManager()
        self.chain_manager.addMCMCUpdaters(self.model,  # substitution model
                                      self.tree,        # tree
                                      self.likelihood,  # likelihood calculation machinery
                                      self.r,           # pseudorandom number generator
                                      True,             # separate_edgelen_params
                                      0,                # max slice units (0 means use largest unsigned int)
                                      1)                # weight for each parameter added
        self.chain_manager.finalize()

    def setupTree(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If tree_topology equals 'file', then a tree is built using the first
        newick tree description that appears in the data file. In this case,
        the program will abort if there are no trees stored in the specified
        data file. If tree_topology is not equal to the string 'file', then
        it is assumed that it contains instead a valid newick tree
        description, and a tree is built from this description.
        
        """
        # If user requested using a tree from the data file, grab the first one stored there
        if self.tree_topology == 'file':        
            newicks = []
            for t in self.reader.getTrees():
                newicks.append(t.newick)
            assert len(newicks) > 0, 'Error: a trees block defining at least one tree must be stored in the nexus data file'
            self.tree_topology = newicks[0]

        # Build a Tree object from the description stored in self.tree_topology
        self.tree = Phylogeny.Tree()
        self.tree.buildFromString(self.tree_topology)
            
        # Make sure that names of tips equal the string equivalent of the tip node number plus 1
        # This means that when tree topologies are sampled, the tree definitions that are output
        # to the .t file will be similar to MrBayes output
        self.ntax = self.data_matrix.getNTax()
        self.nedges = 2*self.ntax - 3
        taxNames = []
        for i in range(self.ntax):
            taxNames.append(str(i+1))
        self.tree.rectifyNames(taxNames)

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
        self.total_evals = 0.0
        self.elapsed_secs = 0.0

        # Set seed if user has supplied one
        if self.random_seed != 'auto':
            self.r.setSeed(int(self.random_seed))

        # Read the data file
        self.reader.readFile(self.data_file_name)
        self.data_matrix = ReadNexus.getDiscreteMatrix(self.reader, 0)

        self.setupTree()        
        self.setupModel()        

        # Add data structures to the nodes of the tree to allow likelihood calculations
        # The structures added to tips allow the tips to store data and transition probability matrices
        # The structures added to the internal nodes allow for storage of transition probability matrices
        # as well as conditional likelihood arrays
        self.likelihood.prepareForLikelihood(self.tree)

        # Compute the current log-likelihood and log-prior
        self.chain_manager.refreshLastLnLike()
        self.chain_manager.refreshLastLnPrior()

        # Create parameter and tree file names based on the data file name
        prefix = os.path.basename(self.data_file_name)        
        self.param_file_name = prefix + '.p'
        self.tree_file_name = prefix + '.t'

        # Open the parameter and tree files
        self.paramFileOpen()
        self.treeFileOpen()
        
        if self.verbose:
            print 'Data file:     ', self.data_file_name
            print 'Prior:         ', self.edgelen_prior_mean
            print 'No. cycles:    ', self.ncycles
            print 'Sample every:  ', self.sample_every
            print 'Tree topology: ', self.tree_topology
            print 'No. samples:   ', self.nsamples
            print 'Sampled trees will be saved in', self.tree_file_name
            print 'Sampled parameters will be saved in', self.param_file_name

            if self.tree.tipNumbersSetUsingNames():
                print 'Tip node numbers were set using the names in the tree description'
            else:
                print 'Warning: tip node numbers were NOT set using the names in the tree description'

            print 'Starting log-likelihood =',self.chain_manager.getLastLnLike()
            print 'Starting log-prior =',self.chain_manager.getLastLnPrior()
            print 'Parameter starting values and prior densities:'
            print '  %-25s  %12s  %12s  %s' % ('Parameter name','Value','Log density','Prior Distribution')
            for p in self.chain_manager.getEdgeLenParams():
                self.showParamInfo(p)
            p = self.chain_manager.getEdgeLenHyperparam()
            self.showParamInfo(p)
            for p in self.chain_manager.getModelParams():
                self.showParamInfo(p)
            print
            #raw_input('stopped')

    def showParamInfo(self, p):
        if p.isMasterParameter():
            print '  %-25s             -  %12.5f  %s' % (p.getName(),p.getLnPrior(),p.getPriorDescr())
        else:
            print '  %-25s  %12.5f  %12.5f  %s' % (p.getName(),p.getCurrValue(),p.getLnPrior(),p.getPriorDescr())
                
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
        tax_labels = self.reader.getTaxLabels()
        for i in range(self.ntax - 1):
            self.treef.write('       %d %s,\n' % (i + 1, tax_labels[i]))
        self.treef.write('       %d %s;\n' % (self.ntax, tax_labels[self.ntax - 1]))

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
        s = p.getSliceSampler()
        nm = p.getName()
        if s.getNumSamples() > 0:
            s.adaptSimple(self.adapt_simple_param)
            self.total_evals += float(s.getNumFuncEvals())
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

    def recordSample(self, cycle, lnL = 0.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Records current tree topology and edge lengths by adding a line to
        the tree file, and records tree length and substitution parameters
        by adding a line to the parameter file.
        
        """
        self.paramf.write('%d\t%.3f\t%.3f' % (cycle + 1, lnL, self.tree.edgeLenSum()))
        self.paramf.write(self.model.paramReport())
        if self.using_hyperprior:
            p = self.chain_manager.getEdgeLenHyperparam()
            self.paramf.write('\t%.5f' % p.getCurrValue())
        self.paramf.write('\n')
        self.treef.write('   tree rep.%d = %s;\n' % (cycle + 1, self.tree.makeNewick()))

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls setup() to finish preparations and then runs the Markov chain.
        Delaying the call to setup() until now allows some parameters to be
        changed after the MCMCSimple object is constructed but before the
        analysis commences.
        
        """
        self.setup()
        self.recordSample(0)
        self.total_evals = 0.0

        print 'Sampling (%d cycles)...' % self.ncycles
        if self.verbose:
            print
        self.timer_start = time.clock()
        last_adaptation = 0
        next_adaptation = self.adapt_first
        for cycle in range(self.ncycles):
            for p in self.chain_manager.getAllUpdaters():
                w = p.getWeight()
                #self.chain_manager.refreshLastLnPrior()
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
        print '%d likelihood evaluations in %.5f seconds' % (self.total_evals, self.elapsed_secs)
        if (self.elapsed_secs > 0.0):
            print '  = %.5f likelihood evaluations/sec' % (self.total_evals/self.elapsed_secs)

        self.treeFileClose()
        self.paramFileClose()
        
if __name__ == '__main__':
    mcmc = MCMCSimple()
    
    mcmc.data_file_name = '../../pyphy/nyldna4.nex'
    #mcmc.tree_topology = '(0:0.1,1:0.1,(2:0.1,3:0.1):0.1)'
    mcmc.tree_topology = '(0:0.1,3:0.1,(2:0.1,1:0.1):0.1)'
    mcmc.ncycles = 500
    mcmc.sample_every = 10
    mcmc.adapt_first = 10
    mcmc.random_seed = '13579'
    mcmc.model_type = 'hky'
    mcmc.using_hyperprior = True
    mcmc.edgelen_prior_mean = 0.1
    mcmc.verbose = True

    mcmc.run()
    
