import os,sys,math
from phycas import *

class GelfandGhosh(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Computes Gelfand-Ghosh on a pre-existing MCMC sample.
    
    """
    def __init__(self, phycas):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Gathers information relevant to Gelfand-Ghosh calculations from the
        Phycas object phycas.
        
        """
        # copy relevant setup information from phycas object
        self.phycas                 = phycas
        self.datafname              = phycas.data_file_name
        self.paramfname             = phycas.gg_pfile
        self.treefname              = phycas.gg_tfile
        self.gg_nreps               = phycas.gg_nreps
        self.gg_kvect               = phycas.gg_kvect
        self.gg_burnin              = phycas.gg_burnin
        self.rnseed                 = phycas.random_seed
        self.num_rates              = phycas.num_rates
        self.gg_bin_patterns        = phycas.gg_bin_patterns
        self.gg_bincount_filename   = phycas.gg_bincount_filename
        self.gg_save_postpreds      = phycas.gg_save_postpreds
        self.gg_postpred_prefix     = phycas.gg_postpred_prefix

        # initialize quantities used in Gelfand-Ghosh calculations
        self.gg_simdata             = Likelihood.SimData()     # temporary container used to hold nascent posterior predictive simulation results until they have been analyzed
        self.gg_binned_simdata      = [0.0]*7                  # if self.gg_bin_patterns is True, this vector of 7 floats is used to summarize the counts in self.gg_simdata
        self.gg_y                   = Likelihood.SimData()     # observed dataset
        self.gg_binned_y            = [0.0]*7                  # if self.gg_bin_patterns is True, this vector of 7 floats is used instead of self.gg_y
        self.gg_mu                  = Likelihood.SimData()     # mean of all posterior predictive datasets
        self.gg_binned_mu           = [0.0]*7                  # if self.gg_bin_patterns is True, this vector of 7 floats is used instead of self.gg_mu
        self.gg_a                   = []            # vector of compromise actions (one for each k in gg_kvect)
        self.gg_npatterns           = []            # vector containing the number of patterns in each posterior predictive dataset
        self.gg_t                   = []            # vector of t values computed from posterior predictive datasets
        self.gg_t_y                 = 0.0           # t for original dataset
        self.gg_t_mean              = 0.0           # mean of t over all posterior predictive datasets
        self.gg_t_mu                = 0.0           # t of mean over all posterior predictive datasets
        self.gg_t_a                 = []            # vector of t values computed from compromise action (one for each k in gg_kvect)
        self.gg_Pm                  = 0.0           # penalty component (same for all k)
        self.gg_Gm                  = []            # vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []            # vector of overall measures (one for each k in gg_kvect)
        self.gg_num_post_pred_reps  = 0.0           # counts total number of posterior predictive simulations performed
        self.gg_total               = 0

        self.lot = ProbDist.Lot()
        if self.rnseed != 0:
            self.lot.setSeed(self.rnseed)
            self.rnseed = self.lot.getSeed()
            
        self.pfile = None
        self.tfile = None
        self.model = None
        self.model_type = None
        self.data_matrix = None
        self.likelihood = None
        self.ntax = 0
        self.nchar = 0
        self.npatterns = 0
        self.ntrees = 0
        self.params = None
        self.trees = None
        self.is_invariable_sites_model = False
        self.is_discrete_gamma_model = False
        
        self.num_edgelen_hyperparams = 0
        self.taxon_labels = []

        self.ok = True
        if self.ok:
            if os.path.exists(self.paramfname):
                self.pfile = file(self.paramfname, 'r')
            else:
                self.ok = False
                print 'Parameter file specified (%s) does not exist' % self.paramfname
                sys.exit()
        if self.ok:
            if os.path.exists(self.paramfname):
                self.tfile = file(self.treefname, 'r')
            else:
                self.ok = False
                print 'Tree file specified (%s) does not exist' % self.treefname
                sys.exit()
                
        self.outputHeader()
        self.setupModel()
        self.getData()
        self.getTrees()

    def outputHeader(self):
        self.phycas.output('Gelfand-Ghosh Analysis:')
        self.phycas.output('  data file:          %s' % self.datafname)
        self.phycas.output('  parameter file:     %s' % self.paramfname)
        self.phycas.output('  tree file:          %s' % self.treefname)
        self.phycas.output('  sims/sample:        %d' % self.gg_nreps)
        self.phycas.output('  k values:')
        for kvalue in self.gg_kvect:
            self.phycas.output('    %f' % kvalue)
        self.phycas.output('  samples skipped:    %d' % self.gg_burnin)
        self.phycas.output('  random number seed: %d' % self.rnseed)
        self.phycas.output()

    def outputDataInfo(self):
        self.phycas.output('  Information about the data:')
        self.phycas.output('    number of taxa:         %d' % self.ntax)
        self.phycas.output('    number of characters:   %d' % self.nchar)
        self.phycas.output('    number of patterns:     %d' % self.npatterns)
        self.phycas.output()
        
    def outputTreesInfo(self):
        self.phycas.output('  Information about the sampled trees:')
        self.phycas.output('    number of trees skipped:  %d' % self.gg_burnin)
        self.phycas.output('    number of trees included: %d' % self.ntrees)
        self.phycas.output('    total trees in file:      %d' % (self.gg_burnin + self.ntrees))
        self.phycas.output()

    def outputModelInfo(self):
        pinvar_str = self.is_invariable_sites_model and '+I' or ''
        dgamma_str = self.is_discrete_gamma_model and '+G' or ''
        
        model_str = '%s%s%s' % (self.model_type, pinvar_str, dgamma_str)
        self.phycas.output('  Information about the substitution model:')
        self.phycas.output('    model name:  %s' % model_str)
        self.phycas.output()

    def getData(self):
        self.data_matrix = phycas.readData(self.datafname)[0]
        self.taxon_labels = data_reader.getTaxLabels()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar()
        assert self.likelihood, 'call setupModel before calling getData'
        self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
        self.npatterns = self.likelihood.getNPatterns()
        self.outputDataInfo()
        
    def getTrees(self):
        tree_reader = ReadNexus.NexusReader()
        tree_reader.readFile(self.treefname)
        self.trees = tree_reader.getTrees()
        self.ntrees = len(self.trees) - self.gg_burnin
        self.outputTreesInfo()

    def setupModel(self):
        # Read second line of the parameter file, which should contain the
        # header line, from which the model can be deduced
        self.params = self.pfile.readlines()
        headers = self.params[1].split()

        # Determine model
        are_basefreqs = ('freqA' in headers)
        are_relrates = ('rAC' in headers)
        self.is_invariable_sites_model = ('pinvar' in headers)
        self.is_discrete_gamma_model = ('shape' in headers)
        
        if 'hyper(all)' in headers:
            self.num_edgelen_hyperparams = 1
        elif 'hyper(external)' in headers:
            assert 'hyper(internal)' in headers
            self.num_edgelen_hyperparams = 2
        is_kappa = ('kappa' in headers)
        if are_basefreqs and are_relrates:
            self.model_type = 'GTR'
            self.model = Likelihood.GTRModel()
        elif are_basefreqs and is_kappa:
            self.model_type = 'HKY'
            self.model = Likelihood.HKYModel()
        else:
            self.model_type = 'JC'
            self.model = Likelihood.JCModel()
        self.outputModelInfo()

        # Create a likelihood object to orchestrate both simulations and likelihood calculations
        self.likelihood = Likelihood.TreeLikelihood(self.model)
                
    def parameterizeModelForSimulation(self, i):
        # i is the index of the MCMC sample being processed
        # Example lines from params file:
        # HKY model: 
        #     Gen        LnL     TL    kappa    freqA    freqC    freqG    freqT  hyper(all)
        #   20000  -5364.230  0.333  2.85931  0.21916  0.29672  0.29253  0.19159     0.31617
        # HKY+gamma:
        #   Gen          LnL     TL    kappa    freqA    freqC    freqG    freqT    shape hyper(all)
        #   20000  -5258.755  0.541  4.41408  0.20816  0.29738  0.27837  0.21609  0.16325    0.21676
        parameters = self.params[i+2].split()
        cycle = int(parameters[0])
        if self.model_type == 'HKY':
            kappa = float(parameters[3])
            piA = float(parameters[4])
            piC = float(parameters[5])
            piG = float(parameters[6])
            piT = float(parameters[7])
            self.model.setKappa(kappa)
            self.model.setNucleotideFreqs(piA, piC, piG, piT)
        elif self.model_type == 'GTR':
            print 'GTR model forthcoming'
            sys.exit()
        #@POL need to handle codon model
        if self.is_invariable_sites_model and self.is_discrete_gamma_model:
            self.model.setPinvar(float(parameters[8]))
            self.model.setShape(float(parameters[9]))
            self.model.setNGammaRates(self.num_rates)
        elif self.is_invariable_sites_model:
            self.model.setPinvar(float(parameters[8]))
        elif self.is_discrete_gamma_model:
            self.model.setShape(float(parameters[8]))
            self.model.setNGammaRates(self.num_rates)

    def calcBinnedT(self, v):
        # Calculate t for the vector v, which should have 7 elements
        assert len(v) == 7, 'vector supplied to the GGImpl.calcBinnedT function should have 7 elements, in this case it had %d elements' % len(v)
        t = 0.0
        n = float(sum(v))
        lnn = math.log(n + 1.0)
        eps = 1.0/7.0
        for x in v:
            term1 = (x + eps)/(n + 1.0)
            term2 = math.log(x + eps) - lnn
            t += term1*term2
        return t
        
    def ggCalculate(self):
        two_n = 2.0*float(self.nchar)

        # Compute the t function for the observed dataset
        if self.gg_bin_patterns:
            self.gg_t_y = self.calcBinnedT(self.gg_binned_y)
        else:                
            self.gg_t_y = self.gg_y.calct(4)

        # gg_mu (or gg_binned_mu) represents the mean over all posterior predictive datasets
        # Compute the t function for this mean dataset
        if self.gg_bin_patterns:
            self.gg_t_mu = self.calcBinnedT(self.gg_binned_mu)
        else:
            self.gg_t_mu = self.gg_mu.calct(4)
        
        # Compute the mean of the t values computed for individual posterior
        # predictive datasets
        self.gg_t_mean = sum(self.gg_t)/float(self.gg_total)

        # Compute the penalty term. Guaranteed to be positive by Jensen's
        # inequality and the convexity of the t function.
        self.gg_Pm = two_n*(self.gg_t_mean - self.gg_t_mu)

        # Loop over k values, computing Gm and Dm for each k value in gg_kvect
        for k in self.gg_kvect:
            # Create a dataset representing the compromise "action"
            if self.gg_bin_patterns:
                a = []
                for b in range(7):
                    next_a = (self.gg_binned_mu[b] + self.gg_binned_y[b])/(k + 1.0)
                    a.append(next_a) 
                t_a = self.calcBinnedT(a)
            else:
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

        self.phycas.output('Pm = %f' % self.gg_Pm)
        for i,k in enumerate(self.gg_kvect):
            self.phycas.output('k = %f:' % k)
            self.phycas.output('  Gm = %f' % self.gg_Gm[i])
            self.phycas.output('  Dm = %f' % self.gg_Dm[i])
        self.phycas.output()

        self.phycas.output('no. patterns in original dataset   = %d' % self.gg_y.getNUniquePatterns())
        self.phycas.output('no. patterns in mean dataset       = %d' % self.gg_mu.getNUniquePatterns())
        sum_npat = 0.0
        for npat in self.gg_npatterns:
            sum_npat += float(npat)
        self.phycas.output('mean no. patterns across datasets  = %f' % (sum_npat/float(len(self.gg_npatterns))))

        self.phycas.output('t for original dataset             = %f' % self.gg_t_y)
        self.phycas.output('t for mean dataset                 = %f' % self.gg_t_mu)
        self.phycas.output('mean of t across datasets          = %f' % self.gg_t_mean)
        for i,k in enumerate(self.gg_kvect):
            self.phycas.output('t of compromise action for k = %.1f = %f' % (k,self.gg_t_a[i]))
            
        ttotal = len(self.gg_t)
        assert ttotal == self.gg_total, 'mismatch between self.gg_total and len(self.gg_t)'
        tsumsq = 0.0
        for t in self.gg_t:
            tsumsq += t*t
        tvar = tsumsq - float(ttotal)*self.gg_t_mean*self.gg_t_mean
        self.phycas.output('std. dev. of t across datasets     = %f' % math.sqrt(tvar))
        for i,k in enumerate(self.gg_kvect):
            self.phycas.output('t of compromise action (k = %6f) = %f' % (k, self.gg_t_a[i]))

    def addPaupBlock(self, fn, tree, pheaders, pvalues):
        # this function not yet tested
        # add a paup block to allow easy estimation of parameters under maximum likelihood
        simf = file(fn, 'a')
        simf.write('\nbegin trees;\n')
        simf.write('  translate\n')
        for num,name in enumerate(self.taxon_labels):
            simf.write("  %d '%s'%s\n" % (num, name, num == self.ntax - 1 and ';' or ','))
        simf.write('  ;\n')
        simf.write('  utree one = %s;\n' % tree.makeNewick())
        simf.write('end;\n')
        simf.write('\nbegin paup;\n')
        simf.write('  log file=%s.log start replace;\n' % fn)
        simf.write('  set criterion=likelihood;\n')
        simf.write('  lset nst=6 basefreq=estimate tratio=estimate rmatrix=estimate pinvar=0.0 rates=gamma shape=estimate;\n')
        simf.write('  lscores 1;\n')
        simf.write('  describe 1 / brlens=sumonly;\n')
        simf.write('  log stop;\n')
        simf.write('  quit;\n')
        simf.write('end;\n')
        simf.write('[\n')
        simf.write('\n')
        simf.write('TL                     = %f\n' % tree.edgeLenSum())
        #if self.using_hyperprior:
        #    for p in self.chain_manager.getEdgeLenHyperparams():
        #        simf.write('edge length hyperparam = %f\n' % p.getCurrValue())
        simf.write('param headers          = %s\n' % pheaders)
        simf.write('param values           = %s\n' % pvalues)
        simf.write(']\n')
        simf.close()

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs simulations and computes the Gelfand-Ghosh measures Pm, Gm,
        and Dm = Pm + Gm for the parameters and trees specified when this
        object was created.
        
        """
        tree = Phylogeny.Tree()

        # Let gg_y contain the observed pattern counts            
        self.likelihood.addDataTo(self.gg_y)
        if self.gg_bin_patterns:
            self.gg_binned_y = self.gg_y.getBinnedCounts()

        if self.gg_bin_patterns and self.gg_bincount_filename:
            binf = open(self.gg_bincount_filename,'w')

        self.phycas.output('Performing posterior-predictive simulations:')
        prev_pct_done = 0.0        
        stopwatch = ProbDist.StopWatch()
        stopwatch.start()
        prev_secs = 0.0
        for i, t in enumerate(self.trees):
            if i >= self.gg_burnin:     # POL use slice to get rid of this waste
                tree_description = t.newick
                tree_name = t.name
                tree_rooted = t.rooted
                
                t.buildTree(tree) 

                #print '*** simulating from tree %s: %s' % (tree_name, tree_description)

                # Set up for simulation
                self.parameterizeModelForSimulation(i)
                    
                # Report progress so user doesn't give up
                pct_done = 100.0*float(i - self.gg_burnin + 1)/float(self.ntrees)
                if pct_done - prev_pct_done >= 10.0:
                    prev_pct_done = pct_done
                    secs = stopwatch.elapsedSeconds()
                    proportion_finished = pct_done/100.0
                    proportion_remaining = 1.0 - proportion_finished
                    eta = secs*proportion_remaining/proportion_finished
                    self.phycas.output('  %.0f%% done (%.1fs remaining)...' % (pct_done, eta))
                    prev_secs = secs
                    cum_caltbinned_secs = 0.0

                # TreeLikelihood object needs to be informed that model has changed
                self.likelihood.replaceModel(self.model)

                # Prepare the tree for simulation (i.e. equip nodes with transition matrices)
                self.likelihood.prepareForSimulation(tree)

                for j in range(self.gg_nreps):
                    self.gg_num_post_pred_reps += 1.0
                    
                    # Simulate from the posterior
                    # Use the function simulateFirst (rather than just simulate) in order
                    # to force recalculation of transition probabilities
                    self.gg_simdata.clear()
                    self.likelihood.simulateFirst(self.gg_simdata, tree, self.lot, self.nchar)

                    # Save the simulated data set if desired
                    if self.gg_save_postpreds:
                        fn = '%s_cycle%d_rep%d.nex' % (self.gg_postpred_prefix, i+1, j)
                        self.gg_simdata.saveToNexusFile(fn, self.taxon_labels, 'dna', ('a','c','g','t'))
                        self.addPaupBlock(fn, tree, self.params[1], self.params[i+2])

                    if self.gg_bin_patterns:
                        # Compute the t function for the simulated dataset                
                        curr_t = self.gg_simdata.calctBinned(4)
                        self.gg_binned_simdata = self.gg_simdata.getBinnedCounts()

                        if self.gg_bincount_filename:
                            bstr = ['%.1f' % x for x in self.gg_binned_simdata]
                            binf.write('%s\tposterior predictive replicate\n' % '\t'.join(bstr))
                    else:
                        # Compute the t function for the simulated dataset                
                        curr_t = self.gg_simdata.calct(4)

                    # Add this value of t to the list (later the mean t will be computed)                
                    self.gg_t.append(curr_t)

                    # Add the number of patterns in self.gg_simdata to the gg_npatterns list
                    self.gg_npatterns.append(self.gg_simdata.getNUniquePatterns())

                    # Update running mean vector gg_mu. A running mean is maintained because
                    # it is easy for the number of counts of constant patterns to overflow
                    # if you wait until the end of the MCMC run to divide by the total.
                    # Here is how the running mean is kept. Assume there will be four numbers
                    # (a, b, c, d) averaged. Thus, the desired quantity is (a+b+c+d)/4.
                    #
                    # gg_num_post_pred_reps   self.gg_mu
                    # ------------------------------------------------------------
                    #           1             a                      = (a)/1
                    #           2             (1/2)a + b/2           = (a+b)/2
                    #           3             (2/3)[(a+b)/2] + c/3   = (a+b+c)/3
                    #           4             (3/4)[(a+b+c)/3] + d/4 = (a+b+c+d)/4
                    # ------------------------------------------------------------
                    #
                    # Note that it is ok if gg_num_post_pred_reps = 1 (in which case
                    # gg_mu is multiplied by zero). Because gg_mu is empty, multBy is
                    # a no-op in this case
                    p = 1.0/self.gg_num_post_pred_reps
                    not_p = 1.0 - p
                    if self.gg_bin_patterns:
                        assert len(self.gg_binned_mu) == 7, 'expecting gg_binned_mu to have length 7, but instead it was %d' % len(self.gg_binned_mu)
                        assert len(self.gg_binned_simdata) == 7, 'expecting gg_binned_simdata to have length 7, but instead it was %d' % len(self.gg_binned_simdata)
                        for b in range(7):
                            self.gg_binned_mu[b] *= not_p
                            self.gg_binned_mu[b] += p*self.gg_binned_simdata[b]
                    else:
                        self.gg_mu.multBy(1.0 - p)
                        self.gg_simdata.multBy(p)
                        self.gg_simdata.addDataTo(self.gg_mu, 1.0)

                    # Increment count of the total number of simulated datasets created
                    # This value is used to later compute the mean t for all simulated datasets
                    # and the mean counts for all simulated data sets
                    self.gg_total += 1

        self.ggCalculate()
                
        # Close files
        self.pfile.close()
        self.tfile.close()

        if self.gg_bin_patterns and self.gg_bincount_filename:
            # write observed bin counts
            bstr = ['%.1f' % x for x in self.gg_binned_y]
            binf.write('%s\tobserved\n' % '\t'.join(bstr))

            # write mu bin counts
            bstr = ['%.1f' % x for x in self.gg_binned_mu]
            binf.write('%s\tmu\n' % '\t'.join(bstr))

            # write compromise action bin counts
            for i,k in enumerate(self.gg_kvect):
                bstr = ['%.1f' % x for x in self.gg_a[i]]
                binf.write('%s\ta for k=%.1f\n' % ('\t'.join(bstr),k))

            binf.close()

        return (self.gg_Pm, self.gg_Gm, self.gg_Dm);        
        
    
