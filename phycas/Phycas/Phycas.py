import os, sys, math, threading, types, copy
import os, sys, math, threading, types, copy
import MCMCManager  # poorly named, as MCMCManager is now only one of many classes within
from phycas.Conversions import *
from phycas.DataMatrix import *
from phycas.Likelihood import *
from phycas.PDFGen import *
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

        # Variables associated with PDF tree drawing (used in pdftree() function)
        # The 14 standard fonts guaranteed to be available in all PDF consumer applications:
        #   Times-Roman      Helvetica             Courier             Symbol
        #   Times-Bold       Helvetica-Bold        Courier-Bold        ZapfDingbats
        #   Times-Italic     Helvetica-Oblique     Courier-Oblique
        #   Times-BoldItalic Helvetica-BoldOblique Courier-BoldOblique
        self.pdf_filename              = 'trees.pdf'    # Set to desired name of pdf file to create
        self.pdf_edge_support_file     = None           # File containing PAUP* output with table of support values; if specified, the support values will be shown on trees plotted
        self.pdf_tip_label_font        = 'Times-Italic' # Font used for tip node names; should be one of the 14 standard fonts listed above
        self.pdf_tip_label_height      = 12             # Height in points of tip node name font
        self.pdf_plot_label_font       = 'Helvetica'    # Font used for plot axis labels; should be one of the 14 standard fonts listed above
        self.pdf_plot_label_height     = 12             # Height in points of plot axis label font
        self.pdf_title_font            = 'Helvetica'    # Font used for scalebar text; should be one of the 14 standard fonts listed above
        self.pdf_title_height          = 14             # Height in points of scalebar text font
        self.pdf_scalebar_position     = 'bottom'       # Valid values are 'top', 'bottom' or None
        self.pdf_scalebar_label_font   = 'Helvetica'    # Font used for scalebar text; should be one of the 14 standard fonts listed above
        self.pdf_scalebar_label_height = 10             # Height in points of scalebar text font
        self.pdf_support_label_font    = 'Times-Roman'  # Font used for edge support values; should be one of the 14 standard fonts listed above
        self.pdf_support_label_height  = 8              # Height in points of edge support font
        self.pdf_support_as_percent    = True           # If True, support values will be shown as percentages (e.g. 93.1) rather than proportions (e.g. 0.931)
        self.pdf_support_decimals      = 1              # The number of decimal places shown in support values (e.g. to get 93.7, specify 1; to round up to 94, specify 0)
        self.pdf_ladderize             = 'right'        # Valid values are 'right', 'left' or None
        self.pdf_page_width            = 8.5            # Page width in inches
        self.pdf_page_height           = 11.0           # Page length in inches
        self.pdf_line_width            = 1.0            # Width of lines representing edges in the tree
        self.pdf_left_margin           = 1.0            # Left margin in inches (1 inch = 72 points)
        self.pdf_right_margin          = 1.0            # Right margin in inches (1 inch = 72 points)
        self.pdf_top_margin            = 1.0            # Top margin in inches (1 inch = 72 points)
        self.pdf_bottom_margin         = 1.0            # Bottom margin in inches (1 inch = 72 points)
        self.pdf_treefile              = None           # Set to tree file name if you want to make one pdf file with each tree from tree file on a separate page
        self.pdf_newick                = None           # Set to the tree description to print if only want to save one tree to a pdf file
        self.pdf_outgroup_taxon        = None           # Set to taxon name of tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)
        

        # Variables associated with the brownian command
        self.brownian_input_tree_file    = None           # Set to the name of the input tree file. This setting should not be None at the time the brownian method is called.

        # Variables associated with SAMC analyses
        self.doing_samc             = False     # If True, using Cheon and Liang "SSAMC" method
        self.samc_move_edgelen_mean = 1.0       # Specifies mean of exponential edge length generation distribution used by SamcMove when new edges are created
        self.samc_move_debug        = False     # If set to True, output will be saved to a file named cf.txt (if it doesn't exist) or, if cf.txt already exists, the output will be compared to cf.txt and the program will halt if a discrepency is found
        self.samc_t0                = 10000.0   # Samc_gain_factor = samc_t0/max(samc_t0, cycle)
        self.samc_move_weight       = 1         # Number of times per cycle that SAMC moves will be performed (currently unused because SAMC moves are not used in standard MCMC analyses)
        self.samc_temperature       = 0.6       # Temperature used in extrapolate move to smooth out differences in probabilities of different possible attachment points

        # Variables associated with path sampling (i.e. thermodynamic integration)
        # If pathsampling() function called, ncycles will be ignored and instead the number of cycles
        # will be ps_burnin + (ps_Q*ps_nbetavals)
        self.ps_toward_posterior    = True      # If True, chain will start with beta = 0.0 (exploring prior) and end up exploring posterior; otherwise, chain will begin by exploring posterior and end exploring prior
        self.ps_burnin              = 1000      # Number of cycles used to equilibrate before increasing beta
        self.ps_Q                   = 100       # Number of cycles between changes in beta
        self.ps_sample_every        = 1         # Determines number of times likelihood will be sampled. Number of samples for each value of beta will be ps_Q/sample_every
        self.ps_nbetavals           = 10        # The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 0, 1/3, 2/3, 1
        self.ps_minbeta             = 0.0       # The first beta value that will be sampled if ps_toward_posterior = True
        self.ps_maxbeta             = 1.0       # The last beta value that will be sampled if ps_toward_posterior = True
        self.ps_filename            = None      # If defined, a file by this name will be created containing the intermediate results (average log-likelihood at each step of the path)

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

    def abort(self, msg):
        s = '\n***** Fatal error: %s' % msg
        #self.output(s)
        sys.exit(s)
        
    def warning(self, msg):
        s = '\n***** Warning: %s' % msg
        self.output(s)
        
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
        self.mcmc_manager.resetNEvals()

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

        total_evals = self.mcmc_manager.getTotalEvals() #self.likelihood.getNEvals()
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

    def pathsampling(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs an MCMC analysis the purpose of which is to obtain an
        accurate estimate of the marginal likelihood of the model by path
        sampling. See Lartillot, N., and H. Philippe. 2006. Syst. Biol. 55(2):
        195-207.
        
        """
        self.check_settings()
        self.phycassert(self.ps_maxbeta > self.ps_minbeta, 'ps_maxbeta must be greater than ps_minbeta')
        self.nchains = 1
        self.ncycles = self.ps_burnin + (self.ps_Q*(self.ps_nbetavals))
        self.ps_delta_beta = (self.ps_maxbeta - self.ps_minbeta)/float(self.ps_nbetavals - 1)
        self.doing_path_sampling = True
        self.setupMCMC()
        self.runMCMC()

    def samc(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a sequential stochastic Markov chain analysis.
        
        """
        self.check_settings()
        self.doing_samc = True;
        self.setupSAMC()
        self.runSAMC()

    def obsolete_simulateDNA(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simulates a DNA dataset and stores it in NEXUS format in the supplied
        filename.
        
        """
        self.check_settings()
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
        
    def obsolete_likelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current tree and current
        model.
        
        """
        self.check_settings()
        self.starting_tree = self.tree_topology
        self.phycassert(self.data_source == 'file', "set data_source to 'file' and specify data_file_name before calling the likelihood function")
        self.readDataFromFile()
        core = MCMCManager.LikelihoodCore(self)
        core.setupCore()
        core.prepareForLikelihood()
        return core.calcLnLikelihood()
        
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

    def unimap(self):
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

    def pdftree(self):
        #complex_outgroup = type(self.pdf_outgroup_taxon) in (types.ListType,types.TupleType)
        self.check_settings()
        simple_outgroup = type(self.pdf_outgroup_taxon) == types.StringType
        self.phycassert(simple_outgroup, 'Phycas cannot yet deal with pdf_outgroup_taxon containing more than one outgroup taxon')
        self.phycassert((self.pdf_treefile and not self.pdf_newick) or (self.pdf_newick and not self.pdf_treefile), 'set either pdf_newick or pdf_treefile, but not both')
        if self.pdf_edge_support_file and os.path.exists(self.pdf_edge_support_file):
            # Read splits file and store all splits found along with their frequencies
            contents_of_file = open(self.pdf_edge_support_file,'r').read()
            regex = re.compile('([*.]+)\s+([0-9.]+)', re.M)
            matches = regex.findall(contents_of_file)
            self.phycassert(matches, 'could not find any splits defined in the pdf_edge_support_file named %s' % self.pdf_edge_support_file)
            self.pdf_splits_to_plot = {}
            for p,f in matches:
                self.pdf_splits_to_plot[p] = float(f)
        if self.pdf_newick:        
            # Build tree the newick description of which is in self.newick
            tree = self.pdf_newick.buildTree()
            
            if self.pdf_outgroup_taxon:
                num = tree.findTipByName(self.pdf_outgroup_taxon)
                self.phycassert(num, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.pdf_outgroup_taxon)
                tree.rerootAtTip(num)
                
            if self.pdf_ladderize:
                if self.pdf_ladderize == 'right':
                    tree.ladderizeRight()
                else:
                    tree.ladderizeLeft()

            # Save tree in PDF  
            pdf = PDFGenerator(self.pdf_page_width, self.pdf_page_height)
            pdf.overwrite = True
            pdf.newPage()
            self.tree2pdf(pdf, tree)
            pdf.saveDocument(self.pdf_filename)
        else:
            # Open pdf_treefile and read trees therein
            self.tree_file_name = self.pdf_treefile
            self.readTreesFromFile()

            # Build each tree and determine its height
            tree = Tree()
            max_height = 0.0
            for tree_def in self.stored_tree_defs:
                tree_def.buildTree(tree)
                tree.rectifyNames(self.taxon_labels)
                if self.pdf_outgroup_taxon:
                    num = tree.findTipByName(self.pdf_outgroup_taxon)
                    self.phycassert(num, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.pdf_outgroup_taxon)
                    tree.rerootAtTip(num)
                h = tree.calcTotalHeight()
                if h > max_height:
                    max_height = h

            # Build each tree again and save in PDF file            
            pdf = PDFGenerator(self.pdf_page_width, self.pdf_page_height)
            pdf.overwrite = True
            for tree_def in self.stored_tree_defs:
                tree_def.buildTree(tree)
                tree.rectifyNames(self.taxon_labels)
                if self.pdf_outgroup_taxon:
                    num = tree.findTipByName(self.pdf_outgroup_taxon)
                    tree.rerootAtTip(num)
                if self.pdf_ladderize:
                    if self.pdf_ladderize == 'right':
                        tree.ladderizeRight()
                    else:
                        tree.ladderizeLeft()
                tree.rectifyNames(self.taxon_labels)
                pdf.newPage()
                self.tree2pdf(pdf, tree, None, max_height)
            pdf.saveDocument(self.pdf_filename)
            
        # Prevent unintentional spillover
        self.pdf_splits_to_plot = None
        
    def tree2pdf(self, pdf, tree, title = None, xscalemax = 0.0, show_support = False):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Prints tree on a pdf object (instance of class PDFGenerator). If title
        is specified, the supplied string will be centered at the top of the
        page. The optional argument xscalemax represents the maximum height
        of a group of trees being printed on separate pages in the same pdf
        document. If xscalemax is left unspecified, each tree will be scaled
        to fit the page and the scalebar will be adjusted accordingly. If
        xscalemax is specified, it will be used to determine the scalebar, and
        the scalebar will remain the same size for all trees printed with the
        same xcalemax value.
        
        """
        # TODO: max_label_points should be calculated outside this function and passed in as an argument
        inch = 72.0
        spacer = 5.0
        max_label_points = 0.0
        rooted_tree = tree.isRooted()
        nodes = []

        # Perform a preorder traversal:
        # 1) for each node, set x-value to height above root (in units of edge length)
        # 2) for each tip, set y-value to tip index, with root tip being 0, and other
        #    tips being numbered from left to right
        # 3) find the length of the longest taxon label as it will be rendered in the
        #    PDF file so that the margin calculations can be made
        # 4) for each internal, just set y-value to 0.0 for now; these internal y-values
        #    will be calculated on the subsequent postorder traversal

        if self.pdf_splits_to_plot:
            tree.recalcAllSplits(tree.getNObservables())
            
        # Record information about the tip serving as the root
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'first preorder node should be the root'
        if not rooted_tree:
            nodes.append(nd)
        subroot = nd.getLeftChild()
        height = subroot.getEdgeLen()
        nd.setX(height) 
        if self.pdf_ladderize and self.pdf_ladderize == 'left':
            last_tip_index = float(tree.getNObservables() - 1)
            nd.setY(last_tip_index) #--> Y is irrelevant if rooted
            ntips = 0.0
        else:
            nd.setY(0.0)
            if rooted_tree:
                ntips = 0.0
            else:
                ntips = 1.0
        max_height = height
        if self.pdf_tip_label_font and not rooted_tree:
            taxon_label = nd.getNodeName()
            label_width = float(self.pdf_tip_label_height)*pdf.calcStringWidth(self.pdf_tip_label_font, taxon_label)
            if label_width > max_label_points:
                max_label_points = label_width
        
        # Record information about the internal node serving as the subroot
        nd = nd.getNextPreorder()
        assert nd.getParent().isRoot(), 'second preorder node should be the subroot'
        nodes.append(nd)
        nd.setX(0.0)
        nd.setY(0.0)
        subroot = nd
        
        # Record information about the remaining nodes in the tree
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                ndpar = nd.getParent()
                nodes.append(nd)
                height = nd.getEdgeLen() + ndpar.getX()
                nd.setX(height)
                if height > max_height:
                    max_height = height
                if nd.isTip():
                    nd.setY(ntips)
                    ntips += 1.0
                    if self.pdf_tip_label_font:
                        taxon_label = nd.getNodeName()
                        label_width = float(self.pdf_tip_label_height)*pdf.calcStringWidth(self.pdf_tip_label_font, taxon_label)
                        if label_width > max_label_points:
                            max_label_points = label_width
                else:
                    nd.setY(0.0)

        # Compute length represented by scale bar. For example,
        #  xscalemax     = 0.00275
        #  log_xscalemax = -2.56
        #  ten_to_power  = 10^floor(-2.56)
        #                = 10^{-3}
        #                = 0.001
        #  scalebar      = 0.001*floor(0.00275/0.001)
        #                = 0.001*floor(2.75)
        #                = 0.002
        #  ndecimals     = -floor(-2.56)
        #                = 3.0
        if xscalemax == 0.0:
            xscalemax = max_height
        half_xscalemax = xscalemax/2.0
        log_xscalemax = math.log10(half_xscalemax)
        ten_to_power = 10**math.floor(log_xscalemax)
        scalebar = ten_to_power*math.floor(half_xscalemax/ten_to_power)
        ndecimals = -int(math.floor(log_xscalemax))
        if ndecimals < 0:
            ndecimals = 0
        format_str = '%%.%df' % (ndecimals)
        scalebar_str = format_str % scalebar
        scalebar_str_extent = float(self.pdf_scalebar_label_height)*pdf.calcStringWidth(self.pdf_scalebar_label_font, scalebar_str)
        scalebar_height = float(self.pdf_scalebar_label_height) + 2*spacer + self.pdf_line_width

        # Find xscaler (amount by which branch lengths must be multiplied to give x-coordinate)
        # and yscaler (amount by which the tip position must be multiplied to give y-coordinate).
        xheight = 0.0
        if self.pdf_tip_label_font:
            xheight = float(self.pdf_tip_label_height)*pdf.getXHeight(self.pdf_tip_label_font)
        half_xheight = xheight/2.0
        ntips = tree.getNObservables()
        label_width   = max_label_points + spacer
        right_margin  = self.pdf_right_margin*inch
        left_margin   = self.pdf_left_margin*inch
        top_margin    = self.pdf_top_margin*inch
        bottom_margin = self.pdf_bottom_margin*inch
        plot_right = self.pdf_page_width*inch
        plot_width = plot_right - left_margin - right_margin
        plot_top = self.pdf_page_height*inch
        plot_height = plot_top - top_margin - bottom_margin

        tree_width = plot_width - label_width
        tree_height = plot_height
        if self.pdf_scalebar_position:
            tree_height -= scalebar_height
        if title:
            tree_height -= 3.0*float(self.pdf_title_height)
        tree_x0 = left_margin
        tree_y0 = bottom_margin + scalebar_height

        xscaler = tree_width/xscalemax
        yscaler = tree_height/float(ntips - 1)

        #pdf.addRectangle(left_margin, bottom_margin, plot_width, plot_height, 1, 'dotted')

        if title and self.pdf_title_height > 0:
            # Draw title centered at top of page
            title_str_extent = float(self.pdf_title_height)*pdf.calcStringWidth(self.pdf_title_font, title)
            title_x = left_margin + (plot_width - title_str_extent)/2.0
            title_y = tree_y0 + tree_height + 2.0*float(self.pdf_title_height)
            pdf.addText(title_x, title_y, self.pdf_title_font, self.pdf_title_height, title)

        if self.pdf_scalebar_position:
            if self.pdf_scalebar_position == 'top':
                # Draw scalebar horizontally starting at top left corner
                scalebar_width = scalebar*xscaler
                scalebar_y = tree_x0 + tree_height - scalebar_height + spacer
                pdf.addLine(left_margin, scalebar_y, left_margin + scalebar_width, scalebar_y, self.pdf_line_width)

                # Draw scalebar text centered above the scalebar
                scalebar_x = left_margin + (scalebar_width - scalebar_str_extent)/2.0
                scalebar_y = tree_x0 + tree_height - float(self.pdf_scalebar_label_height)
                pdf.addText(scalebar_x, scalebar_y, self.pdf_scalebar_label_font, self.pdf_scalebar_label_height, scalebar_str)
            else:
                # Draw scalebar horizontally starting at bottom left corner
                scalebar_width = scalebar*xscaler
                pdf.addLine(left_margin, bottom_margin, left_margin + scalebar_width, bottom_margin, self.pdf_line_width)

                # Draw scalebar text centered above the scalebar
                scalebar_x = left_margin + (scalebar_width - scalebar_str_extent)/2.0
                scalebar_y = bottom_margin + spacer
                pdf.addText(scalebar_x, scalebar_y, self.pdf_scalebar_label_font, self.pdf_scalebar_label_height, scalebar_str)

        # add enough to left margin to center smaller trees horizontally
        left_margin += (xscaler*(xscalemax - max_height) + label_width*(1.0 - max_height/xscalemax))/2.0

        # add enough to the top margin to center smaller trees vertically
        top_margin += (tree_height*(1.0 - max_height/xscalemax))/2.0
        #top_margin += (plot_height*(1.0 - max_height/xscalemax))/2.0

        # adjust yscaler to keep vertical tree dimension proportional to its horizontal dimension
        yscaler *= max_height/xscalemax

        # adjust tip label height (in points) to make size of tip labels commensurate with size of tree
        tip_font_points = self.pdf_tip_label_height*max_height/xscalemax

        # Perform a postorder traversal:
        # 1) scale each x-value
        # 2) calculate y-value of each internal node as the average y-value of its children
        # 3) scale each y-value
        # 4) plot each edge
        # 5) plot names of tips
        # 6) for each internal node, draw shoulder from leftmost child to rightmost
        nodes.reverse()
        for nd in nodes:
            node_x = left_margin + nd.getX()*xscaler
            if nd.isTip():
                node_y = tree_y0 + tree_height - nd.getY()*yscaler
                if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                    node_y -= scalebar_height
                brlen = nd.isRoot() and xscaler*nd.getX() or xscaler*nd.getEdgeLen()
                # draw tip node name
                if self.pdf_tip_label_font:
                    pdf.addText(node_x + spacer, node_y - half_xheight, self.pdf_tip_label_font, tip_font_points, nd.getNodeName())
                # draw line representing edge leading to tip node
                pdf.addLine(node_x, node_y, node_x - brlen, node_y, self.pdf_line_width)
            else:
                nchildren = 1.0
                child = nd.getLeftChild()
                left_child = right_child = child
                childY = child.getY()
                while True:
                    child = child.getRightSib()
                    if child:
                        right_child = child
                        childY += child.getY()
                        nchildren += 1.0
                    else:
                        break
                if (not rooted_tree) and (nd is subroot):
                    if self.pdf_ladderize and self.pdf_ladderize == 'left':
                        right_child = nd.getParent()
                    else:
                        left_child = nd.getParent()
                else:
                    nd.setY(childY/nchildren)
                    node_y = tree_y0 + tree_height - childY*yscaler/nchildren
                    if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                        node_y -= scalebar_height
                    brlen = xscaler*nd.getEdgeLen()
                    # draw line representing edge leading to internal node
                    pdf.addLine(node_x, node_y, node_x - brlen, node_y, self.pdf_line_width)

                # draw line representing shoulders of internal node
                left_y = tree_y0 + tree_height - left_child.getY()*yscaler
                right_y = tree_y0 + tree_height - right_child.getY()*yscaler
                if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                    left_y -= scalebar_height
                    right_y -= scalebar_height
                pdf.addLine(node_x, left_y, node_x, right_y, self.pdf_line_width)

                # if specified, plot support value
                if show_support and self.pdf_splits_to_plot:
                    for p in self.pdf_splits_to_plot.keys():
                        s = Split()
                        s.setOnSymbol('*')
                        s.setOffSymbol('.')
                        s.createFromPattern(p)
                        if s.equals(nd.getSplit()):
                            support_x = node_x + spacer
                            support_y = (left_y + right_y)/2.0 - half_xheight
                            support_str = '%.1f' % self.pdf_splits_to_plot[p]
                            pdf.addText(support_x, support_y, self.pdf_support_label_font, self.pdf_support_label_height, support_str)
                            break
                elif show_support and nd is not subroot:
                    # Expecting each node's support data member to be set already
                    support_format = '%%.%df' % self.pdf_support_decimals
                    if self.pdf_support_as_percent:
                        support_str = support_format % (100.0*nd.getSupport(),)
                    else:
                        support_str = support_format % (nd.getSupport(),)
                    support_str_extent = float(self.pdf_support_label_height)*pdf.calcStringWidth(self.pdf_support_label_font, support_str)
                    support_x = node_x - (brlen + support_str_extent)/2.0
                    support_y = (left_y + right_y)/2.0 + half_xheight
                    pdf.addText(support_x, support_y, self.pdf_support_label_font, self.pdf_support_label_height, support_str)
                    
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
