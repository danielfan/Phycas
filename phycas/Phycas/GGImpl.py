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
        self.outfname               = phycas.gg_outfile
        self.datafname              = phycas.data_file_name
        self.paramfname             = phycas.gg_pfile
        self.treefname              = phycas.gg_tfile
        self.gg_nreps               = phycas.gg_nreps
        self.gg_kvect               = phycas.gg_kvect
        self.gg_burnin              = phycas.gg_burnin
        self.rnseed                 = phycas.random_seed
        self.num_rates              = phycas.num_rates
        self.gg_bin_patterns        = phycas.gg_bin_patterns
        self.gg_save_postpreds      = phycas.gg_save_postpreds
        self.gg_postpred_prefix     = phycas.gg_postpred_prefix

        # for debugging and research
        self.gg_save_spectra        = False    # adds all 256 counts for posterior predictive simulated data sets to a file named spectra.txt, with counts separated by tabs (only use for four-taxon problems)

        # initialize quantities used in Gelfand-Ghosh calculations
        self.gg_simdata             = Likelihood.SimData()     # temporary container used to hold nascent posterior predictive simulation results until they have been analyzed
        self.gg_y                   = Likelihood.SimData()     # observed dataset
        self.gg_mu                  = Likelihood.SimData()     # mean of all posterior predictive datasets
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

        # initialize quantities used for debugging/research code
        self.gg_spectrum           = Likelihood.SimData() # workspace used if gg_save_spectra is True
        self.gg_spectrum_points    = ''        # used for creating surface plot in Maple for spectrum
        self.gg_spectrum_row       = 0

        self.lot = ProbDist.Lot()
        if self.rnseed != 0:
            self.lot.setSeed(self.rnseed)
            self.rnseed = self.lot.getSeed()
            
        self.outfile = None
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
        self.is_flex_model = False
        self.ncat_index = 0
        self.first_flexrate_index = 0
        self.num_edgelen_hyperparams = 0
        self.taxon_labels = []

        self.ok = True
        if self.outfname:
            if os.path.exists(self.outfname):
                self.ok = False
                print 'Output file name (%s) specified for Gelfand-Ghosh analysis already exists' % self.outfname
                print 'Remove or rename it and try again'
                sys.exit()
            else:
                self.outfile = file(self.outfname, 'w')
                self.outputHeader()
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
        self.setupModel()
        self.getData()
        self.getTrees()

    def outputHeader(self):
        if self.outfile:
            self.outfile.write('Gelfand-Ghosh Analysis\n\n')
            self.outfile.write('  data file:          %s\n' % self.datafname)
            self.outfile.write('  parameter file:     %s\n' % self.paramfname)
            self.outfile.write('  tree file:          %s\n' % self.treefname)
            self.outfile.write('  sims/sample:        %d\n' % self.gg_nreps)
            self.outfile.write('  k values:\n')
            for kvalue in self.gg_kvect:
                self.outfile.write('    %f\n' % kvalue)
            self.outfile.write('  samples skipped:    %d\n' % self.gg_burnin)
            self.outfile.write('  random number seed: %d\n' % self.rnseed)
            self.outfile.write('\n')

    def outputDataInfo(self):
        if self.outfile:
            self.outfile.write('  Information about the data:\n\n')
            self.outfile.write('    number of taxa:         %d\n' % self.ntax)
            self.outfile.write('    number of characters:   %d\n' % self.nchar)
            self.outfile.write('    number of patterns:     %d\n' % self.npatterns)
            self.outfile.write('\n')
        
    def outputTreesInfo(self):
        if self.outfile:
            self.outfile.write('  Information about the sampled trees:\n\n')
            self.outfile.write('    number of trees skipped:  %d\n' % self.gg_burnin)
            self.outfile.write('    number of trees included: %d\n' % self.ntrees)
            self.outfile.write('    total trees in file:      %d\n' % (self.gg_burnin + self.ntrees))
            self.outfile.write('\n')

    def outputModelInfo(self):
        if self.outfile:
            pinvar_str = self.is_invariable_sites_model and '+I' or ''
            dgamma_str = self.is_discrete_gamma_model and '+G' or ''
            flex_str = self.is_flex_model and '+FLEX' or ''
            model_str = '%s%s%s%s' % (self.model_type, pinvar_str, dgamma_str, flex_str)
            self.outfile.write('  Information about the substitution model:\n\n')
            self.outfile.write('    model name:  %s\n' % model_str)
            self.outfile.write('\n')

    def getData(self):
        data_reader = ReadNexus.NexusReader()
        data_reader.readFile(self.datafname)
        self.data_matrix = ReadNexus.getDiscreteMatrix(data_reader, 0)
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
        self.ntrees = 0
        for i,t in enumerate(self.trees):
            if i >= self.gg_burnin:
                self.ntrees += 1
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
        if 'rates_probs' in headers:
            self.is_flex_model = True
            self.ncat_index = headers.index('ncat')
            self.first_flexrate_index = headers.index('rates_probs')
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
        # HKY+FLEX model: 
        #   Gen          LnL     TL    kappa    freqA    freqC    freqG    freqT  ncat  hyper(all)  rates_probs
        #   19980  -5259.917  0.570  3.97815  0.21272  0.30260  0.29701  0.18767     2     0.14564  0.14355  3.46689  0.74229  0.25771
        #   20000  -5260.538  0.573  5.39023  0.21662  0.29483  0.28883  0.19972     3     0.21985  0.07716  3.73220  3.75277  0.74855  0.06705  0.18440
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
        if self.is_flex_model:
            self.model.setFlexModel()
            self.num_rates = int(parameters[self.ncat_index])
            self.model.setNGammaRates(self.num_rates)
            for rindex in range(self.num_rates):
                rate = float(parameters[self.first_flexrate_index + rindex])
                self.model.setFlexRateUnnorm(rindex, rate)
                prob = float(parameters[self.first_flexrate_index + self.num_rates + rindex])
                self.model.setFlexProbUnnorm(rindex, prob)
                #print 'rate =', rate,', prob =', prob
            #raw_input('GGImpl.py stopped in parameterizeModelForSimulation')
        elif self.is_invariable_sites_model and self.is_discrete_gamma_model:
            self.model.setPinvar(float(parameters[8]))
            self.model.setShape(float(parameters[9]))
            self.model.setNGammaRates(self.num_rates)
        elif self.is_invariable_sites_model:
            self.model.setPinvar(float(parameters[8]))
        elif self.is_discrete_gamma_model:
            self.model.setShape(float(parameters[8]))
            self.model.setNGammaRates(self.num_rates)
        
    def ggCalculate(self):
        two_n = 2.0*float(self.nchar)

        # Compute the t function for the observed dataset
        if self.gg_bin_patterns:
            self.gg_t_y = self.gg_y.calctBinned(4)
        else:                
            self.gg_t_y = self.gg_y.calct(4)

        # temporary debugging code
        #print 't for observed data =',self.gg_t_y
        #print self.gg_y.patternTable('A C G T'.split())

        # gg_mu is the mean of all the posterior predictive datasets.
        # Compute the t function for the mean dataset
        if self.gg_bin_patterns:
            self.gg_t_mu = self.gg_mu.calctBinned(4)
        else:
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
            if self.gg_bin_patterns:
                t_a = a.calctBinned(4)
            else:
                t_a = a.calct(4)
            self.gg_t_a.append(t_a)
            self.gg_a.append(a)

            # Compute the goodness-of-fit term
            Gkm = (float(k) + 1.0)*two_n*((self.gg_t_mu + k*self.gg_t_y)/(k + 1.0) - t_a)
            self.gg_Gm.append(Gkm)

            # Compute the overall measure            
            Dkm = self.gg_Pm + Gkm
            self.gg_Dm.append(Dkm)

        if self.outfile:
            assert not self.outfile.closed, 'Error: could not save results because the output file is closed'
            self.outfile.write('# Pm = %f\n' % self.gg_Pm)
            for i,k in enumerate(self.gg_kvect):
                self.outfile.write('# k = %f:\n' % k)
                self.outfile.write('#   Gm = %f\n' % self.gg_Gm[i])
                self.outfile.write('#   Dm = %f\n' % self.gg_Dm[i])
            self.outfile.write('\n')

            self.outfile.write('# no. patterns in original dataset   = %d\n' % self.gg_y.getNUniquePatterns())
            self.outfile.write('# no. patterns in mean dataset       = %d\n' % self.gg_mu.getNUniquePatterns())
            sum_npat = 0.0
            for npat in self.gg_npatterns:
                sum_npat += float(npat)
            self.outfile.write('# mean no. patterns across datasets  = %f\n' % (sum_npat/float(len(self.gg_npatterns))))

            self.outfile.write('# t for original dataset             = %f\n' % self.gg_t_y)
            self.outfile.write('# t for mean dataset                 = %f\n' % self.gg_t_mu)
            self.outfile.write('# mean of t across datasets          = %f\n' % self.gg_t_mean)
            ttotal = len(self.gg_t)
            assert ttotal == self.gg_total, 'mismatch between self.gg_total and len(self.gg_t)'
            tsumsq = 0.0
            for t in self.gg_t:
                tsumsq += t*t
            tvar = tsumsq - float(ttotal)*self.gg_t_mean*self.gg_t_mean
            self.outfile.write('# std. dev. of t across datasets     = %f\n' % math.sqrt(tvar))
            for i,k in enumerate(self.gg_kvect):
                self.outfile.write('# t of compromise action (k = %6f) = %f\n' % (k, self.gg_t_a[i]))

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
        Performs simulations and computes the Gelfand-Ghosh measures Pm, Gm,
        and Dm = Pm + Gm for the parameters and trees specified when this
        object was created.
        
        """
        tree = Phylogeny.Tree()

        # Let gg_y contain the observed pattern counts            
        self.likelihood.addDataTo(self.gg_y)

        if self.gg_bin_patterns:
            # temporary
            binf = open('bins.txt','w')

        #    # If saving spectra, save the spectrum from the original data set
        #    if self.gg_save_spectra:
        #        assert self.ntax == 4, 'gg_save_spectra is designed for 4-taxon problems only (i.e. ntax = 4); ntax is %d in this case' % self.ntax
        #        self.fillSpectrum()
        #        self.gg_spectrum.zeroCounts()
        #        self.gg_y.addDataTo(self.gg_spectrum, 1.0)
        #        yobs_row = self.gg_spectrum.createMapleTuples(self.gg_spectrum_row, 100)
        #        self.gg_spectrum_points += yobs_row
        #        self.gg_spectrum_row += 1
        #        self.gg_spectrum.appendCountsToFile('spectra.txt', False)
        #        for z in range(9):
        #            self.gg_spectrum_points += ','
        #            self.gg_spectrum_points += yobs_row
        #            self.gg_spectrum_row += 1
        #            self.gg_spectrum.appendCountsToFile('spectra.txt', False)

        print 'Performing posterior-predictive simulations:'        

        prev_pct_done = 0.0        
        for i,t in enumerate(self.trees):
            if i >= self.gg_burnin:
                tree_description = t.newick
                # build the next tree
                tree_name = t.name
                tree_rooted = t.rooted
                
                # Second argument (zero-based tip node numbers) must be True
                # because we read these tree descriptions from a nexus file, and
                # the nexus reader always returns zero-based tree descriptions
                tree.buildFromString(tree_description, True) 

                #print '*** simulating from tree %s: %s' % (tree_name, tree_description)

                # set up for simulation
                self.parameterizeModelForSimulation(i)
                    
                # report progress so user doesn't give up
                pct_done = 100.0*float(i - self.gg_burnin + 1)/float(self.ntrees)
                if pct_done - prev_pct_done >= 10.0:
                    prev_pct_done = pct_done
                    print '  %.0f%% done...' % (pct_done)

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

                    # debugging/research code
                    if self.gg_save_spectra:
                        self.gg_spectrum.zeroCounts()
                        self.gg_simdata.addDataTo(self.gg_spectrum, 1.0)
                        self.gg_spectrum_points += ','
                        self.gg_spectrum_points += self.gg_spectrum.createMapleTuples(self.gg_spectrum_row, 100)
                        self.gg_spectrum_row += 1
                        self.gg_spectrum.appendCountsToFile('spectra.txt', False)

                    # Save the simulated data set if desired
                    if self.gg_save_postpreds:
                        fn = '%s_cycle%d_rep%d.nex' % (self.gg_postpred_prefix, i+1, j)
                        self.gg_simdata.saveToNexusFile(fn, self.taxon_labels, 'dna', ('a','c','g','t'))
                        self.addPaupBlock(fn, tree, self.params[1], self.params[i+2])
                        
                    if self.gg_bin_patterns:
                        # temporary
                        b = self.gg_simdata.getBinnedCounts()
                        bstr = ['%.1f' % x for x in b]
                        binf.write('%s\tpostpred rep.\n' % '\t'.join(bstr))
                        
                        # Compute the t function for the simulated dataset                
                        curr_t = self.gg_simdata.calctBinned(4)
                    else:
                        # Compute the t function for the simulated dataset                
                        curr_t = self.gg_simdata.calct(4)

                    # temporary debugging code
                    #print '\nt =',curr_t,', patterns =',self.gg_simdata.getNUniquePatterns()
                    #print self.gg_simdata.patternTable('A C G T'.split())

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
                    # gg_mu is multiplied by zero) because multBy is a no-op in this
                    # case since gg_mu is empty
                    p = 1.0/self.gg_num_post_pred_reps
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
        if self.outfile:
            self.outfile.close()

        if self.gg_bin_patterns:
            # temporary
            b = self.gg_y.getBinnedCounts()
            bstr = ['%.1f' % x for x in b]
            binf.write('%s\t(observed)\n' % '\t'.join(bstr))
            binf.close()

        return (self.gg_Pm, self.gg_Gm, self.gg_Dm);        
        
    