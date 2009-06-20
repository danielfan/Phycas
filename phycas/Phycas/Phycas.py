import os, sys, math, threading, types, copy
import MCMCManager  # poorly named, as MCMCManager is now only one of many classes within
from phycas.Conversions import *
from phycas.Likelihood import *
#from phycas.PDFGen import *
from phycas.Phylogeny import *
from phycas.ProbDist import *
from phycas.ReadNexus import *

# see http://mail.python.org/pipermail/python-list/2002-January/121376.html
import inspect

class Phycas(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Performs Bayesian phylogenetic MCMC analyses. The tree topology and
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
        self.quiet                  = False     # If True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well

        # Variables associated with the brownian command
        self.brownian_input_tree_file    = None           # Set to the name of the input tree file. This setting should not be None at the time the brownian method is called.

        # Variables associated with SAMC analyses
        self.doing_samc             = False     # If True, using Cheon and Liang "SSAMC" method
        self.samc_move_edgelen_mean = 1.0       # Specifies mean of exponential edge length generation distribution used by SamcMove when new edges are created
        self.samc_move_debug        = False     # If set to True, output will be saved to a file named cf.txt (if it doesn't exist) or, if cf.txt already exists, the output will be compared to cf.txt and the program will halt if a discrepency is found
        self.samc_t0                = 10000.0   # Samc_gain_factor = samc_t0/max(samc_t0, cycle)
        self.samc_move_weight       = 1         # Number of times per cycle that SAMC moves will be performed (currently unused because SAMC moves are not used in standard MCMC analyses)
        self.samc_temperature       = 0.6       # Temperature used in extrapolate move to smooth out differences in probabilities of different possible attachment points

        # Variables associated with Gelfand-Ghosh calculation
        self.gg_outfile             = 'gg.txt'  # File in which to save gg results (use None to not save results)
        self.gg_nreps               = 1         # The number of replicate simulations to do every MCMC sample
        self.gg_kvect               = [1.0]     # Vector of k values to use when computing Gm and Dm
        self.gg_save_postpreds      = False     # If True, all posterior predictive data sets will be saved
        self.gg_postpred_prefix     = 'pp'      # Prefix to use for posterior predictive dataset filenames (only used if gg_save_postpreds is True)
        self.gg_burnin              = 1         # Number of starting samples to skip when computing Gelfand-Ghosh measures
        self.gg_pfile               = None      # Name of parameter file to use for Gelfand-Ghosh calculations
        self.gg_tfile               = None      # Name of tree file to use for Gelfand-Ghosh calculations
        self.gg_bin_patterns        = False     # If True, patterns will be classified into 7 bins, corresponding to 'A only', 'C only', 'G only', 'T only', 'any 2 states', 'any 3 states' and 'any 4 states'. Gelfand-Ghosh statistics will be computed on this vector of counts instead of the complete vector of pattern counts. Can only be used for DNA/RNA data.
        self.gg_bincount_filename   = None      # If not None, and if gg_bin_patterns is True, the binned counts for the original dataset and all posterior predictive data sets will be saved to a file by this name

        # ***** IT IS BEST NOT TO CHANGE ANYTHING BELOW HERE *****
        self.debugging              = False      # If set to True expect lots of debug output (e.g. data pattern table)
        self.data_matrix            = None
        self.file_name_data_stored  = None
        self.file_name_trees_stored = None
        self.do_marginal_like       = False
        self.mcmc_manager           = MCMCManager.MCMCManager(self)
        self.heat_vector            = None      # Leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using self.nchains and self.heating_lambda
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
        self.tmp_simdata            = SimData()
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
        self.ss_delta_beta          = 0.0
        self.doing_steppingstone_sampling = False
        self.path_sample            = None
        self.psf                    = None
        #self.pdf_splits_to_plot     = None
        self.param_file_name        = None  
        self.tree_file_name         = None
        self.nsamples               = None
        self.ss_beta                = 1.0
        self.wangang_sampled_betas  = None
        self.wangang_sampled_likes  = None
        self.unimap_manager         = None
        self.nsamples               = 0

        # make a copy of the vector of keys from __dict__ so that we can detect (in check_settings)
        # whether the user has accidentally introduced a new variable (by misspelling, for example)
        self.dict_keys = copy.copy(self.__dict__.keys())        
        
    # see http://mail.python.org/pipermail/python-list/2002-January/121376.html
    def source_line():
        return inspect.getouterframes(inspect.currentframe())[1][2]

    def check_settings(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        At the end of the __init__ function, a variable dict_keys is
        created consisting of a copy of all keys in __dict__. Now,
        we expect len(__dict__) to be one greater than len(dict_keys)
        because __dict__ now includes an entry for dict_keys. If 
        len(__dict__) is longer than this, then probably the user
        misspelled one of the variable names, thus accidentally
        adding another entry to __dict__. This gives us an 
        opportunity to catch this kind of mistake.
        
        """
        if len(self.__dict__) > len(self.dict_keys) + 1:
            for k in self.__dict__.keys():
                if k not in self.dict_keys and not k == 'dict_keys':
                    print 'Error:',k,'is not a valid Phycas setting'

                    f = open('dist_keys.txt', 'w')
                    for kk in self.dict_keys:
                        f.write('%s\n' % kk)
                    f.close()
                    
                    f = open('__dict__.txt', 'w')
                    for kk in self.__dict__.keys():
                        f.write('%s\n' % kk)
                    f.close()

                    sys.exit(0)
        
    def setEdgelenPrior(self, dist):
        self.internal_edgelen_prior = self.external_edgelen_prior = dist
        
    def getEdgelenPrior(self):
        self.phycassert(self.internal_edgelen_prior is self.external_edgelen_prior, "There are separate distributions for internal and external edge lengths")
        return self.internal_edgelen_prior

    edgelen_prior = property(getEdgelenPrior, setEdgelenPrior)
    
    def runSAMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the SSAMC analysis.
        
        """
        if self.samc_move_debug:
            stop_at_cycle = -1  # set this to the cycle you want to examine in detail, or -1 to ignore

        nsamples = 0
        max_lnL = None
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

        self.stopwatch.start()
        self.mcmc_manager.resetNumLikelihoodEvals()

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

            if current_level == max_level:
                nsamples += 1
                if not max_lnL or lnL > max_lnL:
                    max_lnL = lnL
                self.mcmc_manager.recordSample(cycle)
                outstr = "nsamples = %d, cycle = %d, lnL = %f, best = %f" % (nsamples,cycle,lnL,max_lnL)
                print outstr
            else:
                assert current_level < max_level, 'max_level is not max. level'
                outstr = "cycle = %d, level = %d: lnL =%f %s" % (cycle, current_level, lnL, c)
            #outstr = "cycle = %d, level = %d: lnL =%*f %s" % (cycle, current_level, 15*current_level, lnL, c)
            #print outstr

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

        total_evals = self.mcmc_manager.getTotalEvals()
        total_secs = self.stopwatch.elapsedSeconds()
        self.output('%d likelihood evaluations in %.5f seconds' % (total_evals, total_secs))
        if (total_secs > 0.0):
            self.output('  = %.5f likelihood evaluations/sec' % (total_evals/total_secs))
    
    def readTreesFromFile(self):
        if not self.file_name_trees_stored or (self.tree_file_name != self.file_name_trees_stored):
            self.reader.readFile(self.tree_file_name)
            self.taxon_labels = self.reader.getTaxLabels()  # shouldn't overwrite taxon_labels stored previously
            self.stored_tree_defs = self.reader.getTrees()
            self.phycassert(len(self.stored_tree_defs) > 0, 'expecting a trees block defining at least one tree in the nexus data file %s' % self.tree_file_name)
            self.file_name_trees_stored = self.tree_file_name    # prevents rereading same tree file later

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
            #print "level", level, " self.addition_sequence =>", self.addition_sequence
            
    def getSamcStartingTree(self):
        #POL TEMP: not the best way to deal with this. See LikelihoodCore.__init__()
        r = Lot()
        r.setSeed(int(self.random_seed))
        self.starting_edgelen_dist.setLot(r)

        addseq = self.addition_sequence
        brlens = [self.starting_edgelen_dist.sample() for i in range(5)]
        self.tree_topology = Newick("((%d:%.5f,%d:%.5f):%.5f,%d:%.5f,%d:%.5f)" % (
                             addseq[0] + 1, brlens[0], addseq[1] + 1, brlens[1], brlens[4], addseq[2] + 1,
                             brlens[2], addseq[3] + 1, brlens[3]), Newick.ONE_BASED_TAXA_NUMBERS)
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
        
    #def setLogFile(self, filename):
    #    # Open a log file if requested
    #    if self.logf:
    #        self.logf.close()
    #        self.logf = None
    #    if not filename:
    #        self._logFileName = None
    #    else:
    #        # TODO check first to see if it exists before blindly overwriting
    #        self.logf = file(filename, 'w')
    #        self._logFileName = filename
    #
    #def getLogFile(self):
    #    return self._logFileName
    #    
    #log_file_name = property(getLogFile, setLogFile)

	#     def pathsampling(self):
	#         #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
	#         """
	#         Performs an MCMC analysis the purpose of which is to obtain an
	#         accurate estimate of the marginal likelihood of the model by path
	#         sampling. See Lartillot, N., and H. Philippe. 2006. Syst. Biol. 55(2):
	#         195-207.
	#         
	#         """
	#         self.check_settings()
	#         self.phycassert(self.ss_maxbeta > self.ss_minbeta, 'ss_maxbeta must be greater than ss_minbeta')
	#         self.nchains = 1
	#         self.ncycles = self.ss_burnin + (self.ss_Q*(self.ss_nbetavals))
	#         self.ss_delta_beta = (self.ss_maxbeta - self.ss_minbeta)/float(self.ss_nbetavals - 1)
	#         self.doing_steppingstone_sampling = True
	#         self.setupMCMC()
	#         self.runMCMC()

    def samc(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a sequential stochastic Markov chain analysis.
        
        """
        self.check_settings()
        self.doing_samc = True;
        self.setupSAMC()
        self.runSAMC()

    def likelihoods(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current model of all trees
        in the file whose name is stored in tree_file_name.
        
        """
        self.check_settings()
        self.phycassert(len(self.tree_file_name) > 0, 'specify tree_file_name before calling the likelihoods function')
        self.readTreesFromFile()
        self.phycassert(self.data_source == 'file', "set data_source to 'file' and specify data_file_name before calling the likelihoods function")
        self.readDataFromFile()
        for t, topology in enumerate(self.stored_tree_defs):
            self.starting_tree = topology
            core = MCMCManager.LikelihoodCore(self)
            core.setupCore(True)    # specify zero_based_tips = True because topology came from file
            core.prepareForLikelihood()
            print 'Setting all edge lengths to 0.1'
            core.tree.setAllEdgeLens(0.1)
            print 'length of tree %d = %.5f' % (t,core.tree.edgeLenSum())
            print 'log-likelihood of tree %d = %.5f' % (t, core.calcLnLikelihood())
        
    def gg(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes Gelfand-Ghosh on a pre-existing MCMC sample defined in the
        files self.gg_pfile and self.gg_tfile.
        
        """
        self.check_settings()
        self.phycassert(self.gg_pfile, 'gg_pfile cannot be None if gg function called')
        self.phycassert(self.gg_tfile, 'gg_pfile cannot be None if gg function called')
        import GGImpl
        gelfand_ghosh = GGImpl.GelfandGhosh(self)
        self.gg_Pm, self.gg_Gm, self.gg_Dm = gelfand_ghosh.run()
        return (self.gg_Pm, self.gg_Gm, self.gg_Dm)

    def obsolete_unimap(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a uniformized mapping MCMC analysis. 
        
        """
        self.check_settings()
        import UnimapImpl
        self.unimap_manager = UnimapImpl.UnimapManager(self)
        self.unimap_manager.run()

    def brownian(self):
        self.check_settings()
        import BrownianImpl
        brownian = BrownianImpl.Brownian(self)
        brownian.run()

    def calcDistances(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes a matrix of pairwise distances between sequences.
        
        """
        self.check_settings()
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
    #CPPCompiledInDebug = False

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
