import os,sys,math,random
from phycas import *
from MCMCManager import MCMCManager
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class InflatedDensityRatio(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Estimates marginal likelihood using the IDR method and saves result 
    in the log file. Method outlined in the following paper:
    Arima, Serena, and Luca Tardella. 2010. An alternative marginal 
    likelihood estimator for phylogenetic models. arXiv:1001:2136v2
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes InflatedDensityRatio object.
        
        """
        CommonFunctions.__init__(self, opts)
        self.sample_size = None         # number of tree definitions and parameter samples
        self.stored_trees = None        # list of trees built from tree definitions in the trees file
        self.param_file_lines = None    # list of lines from the params file (header lines excluded)
        self.starting_tree = None       # the tree to be processed
        self.data_matrix = None         # the data matrix itself
        self.taxon_labels = None        # list of taxon labels from the data matrix
        self.ntax = None                # number of taxa in the data matrix
        self.nchar = None               # number of sites in the data matrix
        self.param_names = None         # list of the names of the parameters
        self.tree_objects = None        # dictionary in which keys are unique tree identifiers, and values are tree objects
        self.parameters = None          # dictionary in which keys are unique tree identifiers, and values are lists
                                        #   e.g. self.parameters[tree_id][j] = {'rAG': -4.2432, 'freqC': -2.3243, ...}
                                        #   where tree_id is a list of split representations uniquely identifying a particular tree
                                        #   and j is the index of the (j+1)th sample pertaining to that tree. The keys in the 
                                        #   dictionary (e.g. 'rAG') are parameter names and the values are transformed parameter values
        self.edge_lengths = None        # dictionary in which keys are unique tree identifiers, and values are lists
                                        #   e.g. self.edge_lengths[tree_id][j] = {'-**-**--': -1.2345, '-*****--':-2.14225, ...}
                                        #   where tree_id is a list of split representations uniquely identifying a particular tree
                                        #   and j is the index of the (j+1)th sample from that tree. The keys in the dictionary (e.g. '-**-**--') 
                                        #   are string representations of splits and the values are log-transformed edge lengths
        self.models = None              # list of models used in the defined partition
        self.model_names = None         # list of names of models used in the defined partition

        # data members below are needed because must use an MCMCManager to compute posteriors
        # This is overkill and will be simplified later
        self.mcmc_manager = None        # manages MarkovChain used to compute posterior
        self.heat_vector = [1.0]        # specifies that just one chain will be created with power 1.0
        self.curr_treeid = None         # specifies the current tree id (tuple containing string representations of all splits)
        
    def storeTrees(self, input_trees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Open specified tree file and read trees therein, storing the list of
        tree objects in self.stored_trees.
        
        """
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_trees = list(input_trees)   # this triggers reading the tree file
        self.taxon_labels = input_trees.taxon_labels # this line must follow the coercion of the trees to a list
        num_stored_trees = len(self.stored_trees)
        self.stdout.phycassert(num_stored_trees > 0, 'Specified tree source (%s) contained no stored trees' %  str(input_trees))
        
    def storeParams(self, input_params):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Open specified params file and read parameters therein, storing the
        list of lines that represent sampled parameter vectors in 
        self.param_file_lines and the header (containing the names of the 
        parameters) in the list self.param_names.
        
        """
        burnin = self.opts.burnin
        self.param_file_lines = open(input_params, 'r').readlines()
        self.stdout.phycassert(len(self.param_file_lines) >= 3 + burnin, "File '%s' does not look like a parameter file (too few lines)")
        self.param_names = self.param_file_lines[1].split()
        self.stdout.phycassert(self.param_names[1] != 'beta', "File '%s' appears to be the result of a stepping-stone analysis. This method requires a sample from the posterior (not power posterior) distribution." % input_params)

    def fillParamDict(self, param_vect):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build and return a dictionary containing the parameter names (keys)
        and transformed parameter values (values) for the sample supplied in
        param_vect.
        
        """
        param_dict = {}
        
        param_list = param_vect.split()
        assert len(param_list) == len(self.param_names), 'Number of parameters in sample (%d) of parameter file differs from the number of parameters listed in the header line (%d)' % (len(param_list).len(self.param_names))
        
        log_freqA = None
        log_rAC = None
        for h,p in zip(self.param_names,param_list):
            if 'brlen' in h:
                # branch length parameters are gleaned from trees, so ignore the (redundant) edge lengths in the params file
                pass
            elif h in ['Gen', 'lnL', 'lnPrior', 'TL']:
                # quantities that should be saved but not transformed
                param_dict[h] = p
            elif 'freqA' in h:
                # first nucleotide frequency in a subset is not stored but is used to transform the other frequencies
                assert float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p))
                log_freqA = math.log(float(p))
            elif ('freqC' in h) or ('freqG' in h) or ('freqT' in h):
                # nucleotide frequencies for C, G and T are transformed and stored
                assert float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p))
                param_dict[h] = math.log(float(p)) - log_freqA
            elif 'rAC' in h:
                # first GTR exchangeability in a subset is not stored but is used to transform the other exchangeabilities
                assert float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p))
                log_rAC = math.log(float(p))
            elif ('rAG' in h) or ('rAT' in h) or ('rCG' in h) or ('rCT' in h) or ('rGT' in h):
                # GTR exchangeabilities other than rAC are transformed and stored
                assert float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p))
                param_dict[h] = math.log(float(p)) - log_rAC
            else:
                # everything else is log-transformed and stored
                assert float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p))
                param_dict[h] = math.log(float(p))
        return param_dict

    def fillEdgeLenDict(self, tree):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build and return a tuple containing: (1) a dictionary containing the 
        string representations (keys) and log-transformed edge lengths 
        (values) for the tree definition supplied in tree_def; and (2) the
        tree id, which is a tuple of split representations for this tree that
        can be used to uniquely identify the tree topology.
        
        """
        edgelen_dict = {}
        
        #tree = Phylogeny.Tree()
        #tree_def.buildTree(tree)
        tree.rectifyNames(self.taxon_labels)
        ntips = tree.getNObservables()
        tree.recalcAllSplits(ntips)
        
        # Traverse the tree
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'The first preorder node should be the root'
        treelen = 0.0
        assert tree.hasEdgeLens(), 'Sampled tree has no edge lengths'
        split_list = []    # this will be a list of split representations (strings) that uniquely identifies the tree
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                # Determine whether this split represents an internal or tip node
                is_tip_node = nd.isTip() or nd.getParent().isRoot()

                # Grab the edge length
                edge_len = nd.getEdgeLen()
                treelen += edge_len
                
                # Grab the split and invert it if necessary to attain a standard polarity
                s = nd.getSplit()
                if s.isBitSet(0):
                    s.invertSplit()
                    
                # Create a string representation of the split
                ss = s.createPatternRepresentation()
                split_list.append(ss)

                # Add string represention of the split to the tree_key list, which
                # will be used to uniquely identify the tree topology
                #if not is_tip_node:                        
                #    tree_key.append(ss)

                # Create a dictionary entry with ss as key and the log-transformed edge_len as value
                edgelen_dict[ss] = math.log(edge_len)
                
        #print 'tree ==>',tree.newick
        #for ss in edgelen_dict.keys():
        #    edgelen = math.exp(edgelen_dict[ss])
        #    print ss,edgelen
        #raw_input('check edge lengths...')
        return edgelen_dict, tuple(split_list)
        
    def harvestTreesAndParams(self, input_trees, input_params):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build dictionaries self.parameters and self.edge_lengths from 
        supplied input_trees and input_params.
        
        """
        self.storeTrees(input_trees)
        self.storeParams(input_params)
        
        # Initialize the data member that will hold the newick-format tree definitions for each tree sampled
        self.tree_objects = {}
        
        # Initialize the data member that will hold all information about edge lengths
        self.edge_lengths = {}
        
        # Initialize the data member that will hold all information about parameter samples
        self.parameters = {}
        
        self.sample_size = 0
        post_burnin_paramvects = self.param_file_lines[self.opts.burnin+3:] # eliminate burnin samples as well as 2 header lines and starting state (generation 0)
        post_burnin_trees = self.stored_trees[self.opts.burnin+1:]   # eliminate burnin samples as well as starting state (generation 0)
        for tree, param_vect in zip(post_burnin_trees,post_burnin_paramvects): 
            self.sample_size += 1

            # This dictionary will be filled with parameter values: e.g. sample_dict['freqA'] = 0.24335
            param_dict = self.fillParamDict(param_vect)
            
            # This dictionary will be filled with log-transformed edge lengths: e.g. edgelen_dict[s] = -1.10293
            # where the split object s associated with the edge is used as the key
            edgelen_dict, tree_id = self.fillEdgeLenDict(tree)
            
            # add tree to the self.tree_objects dictionary if not already present
            if not tree_id in self.tree_objects.keys():
                self.tree_objects[tree_id] = tree
                
            # add edgelen_dict to the appropriate list in the self.edge_lengths dictionary
            if tree_id in self.edge_lengths.keys():
                # this tree has already been seen, so add to list already established for this tree
                self.edge_lengths[tree_id].append(edgelen_dict)
            else:
                # this tree has not yet been seen, so start a new list
                self.edge_lengths[tree_id] = [edgelen_dict]
                
            # add param_dict to the appropriate list in the self.parameters dictionary
            if tree_id in self.parameters.keys():
                # this tree has already been seen, so add to list already established for this tree
                self.parameters[tree_id].append(param_dict)
            else:
                # this tree has not yet been seen, so start a new list
                self.parameters[tree_id] = [param_dict]
    
    def loadData(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Loads data from data_source and sets data members taxon_labels, ntax,
        nchar, and data_matrix.
        
        """
        ds = self.opts.data_source
        self.phycassert(ds is not None, 'Data source is not allowed to be None for an IDR analysis')
        self.data_matrix = ds.getMatrix()
        self.phycassert(self.data_matrix is not None, 'Data matrix could not be input')
        self.taxon_labels = self.data_matrix.getTaxLabels()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar()
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")
        self.phycassert(self.ntax > 0, 'Number of taxa in data matrix was 0')

    def checkModel(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Freezes model, creating a default partition if partition not defined
        by user, then performs various sanity checks to make sure model is
        internally consistent.
        
        """
        # Create a default partition if a partition was not defined by the user
        partition.validate(self.nchar)
        
        # Ask for a partition report, passing self as the reporter (object that has an output function)
        partition.partitionReport(self)

        self.models         = [m for (n,s,m) in partition.subset]
        self.model_names    = [n for (n,s,m) in partition.subset]
        
        # Perform sanity checks on models
        for m,n in zip(self.models, self.model_names):
            #print '==> checking model %s' % n
            bad_priors = m.checkPriorSupport()
            self.phycassert(len(bad_priors) == 0, 'In model %s, prior support is incorrect for these parameters:\n%s' % (n,'  \n'.join([p for p in bad_priors])))
            
            if m.edgelen_prior is not None:
                # set both internal and external edge length priors to edgelen_prior
                m.internal_edgelen_prior = m.edgelen_prior
                m.external_edgelen_prior = m.edgelen_prior
            else:
                # Ensure that user has specified both internal and external edge length priors
                self.phycassert(m.internal_edgelen_prior is not None, 'In model %s, internal_edgelen_prior cannot be None if edgelen_prior is None' % n)
                self.phycassert(m.external_edgelen_prior is not None, 'In model %s, external_edgelen_prior cannot be None if edgelen_prior is None' % n)
                
            if m.edgelen_hyperprior is not None:
                # Ensure that both internal and external edgelen priors are Exponential
                if m.internal_edgelen_prior.getDistName() != 'Exponential':
                    m.internal_edgelen_prior = Exponential(1.0)
                    self.warning('In model %s, internal_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)
                if m.external_edgelen_prior.getDistName() != 'Exponential':
                    m.external_edgelen_prior = Exponential(1.0)
                    self.warning('In model %s, external_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)
        
    def getStartingTree(self):
        self.starting_tree = self.tree_objects[self.curr_treeid]
        return self.starting_tree
        
    def checkPosterior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculates log-posterior for sampled parameters and edge_lengths and
        compares result to the lnL and lnPrior values from the original MCMC
        analysis. This function is intended as a debugging tool.
        
        """
        self.checkModel()
        
        for tid in self.parameters.keys():
            self.curr_treeid = tid
            
            # create chain manager and (single) chain
            self.mcmc_manager = MCMCManager(self)
            self.mcmc_manager.createChains()
            c = self.mcmc_manager.getColdChain()
            
            #updaters = {}
            #for m in c.chain_manager.getAllUpdaters():
            #    nm = m.getName()
            #    updaters[nm] = m

            # updater names:
            #   state_freqs
            #   relrates
            #   larget_simon_local
            #   master_edgelen
            #   gamma_shape
            
            for e,p in zip(self.edge_lengths[tid],self.parameters[tid]):
                # set parameters
                if 'rAG' in p.keys():
                    rAG_on_rAC = math.exp(p['rAG'])
                    rAT_on_rAC = math.exp(p['rAT'])
                    rCG_on_rAC = math.exp(p['rCG'])
                    rCT_on_rAC = math.exp(p['rCT'])
                    rGT_on_rAC = math.exp(p['rGT'])
                    rAC = 1.0/(1.0 + rAG_on_rAC + rAT_on_rAC + rCG_on_rAC + rCT_on_rAC + rGT_on_rAC)
                    rAG = rAC*rAG_on_rAC
                    rAT = rAC*rAT_on_rAC
                    rCG = rAC*rCG_on_rAC
                    rCT = rAC*rCT_on_rAC
                    rGT = rAC*rGT_on_rAC
                    print 'GTR rates:',rAC,rAG,rAT,rCG,rCT,rGT
                    #u = updaters['relrates']
                    m = c.partition_model.getModel(0)
                    m.setRelRates([rAC,rAG,rAT,rCG,rCT,rGT])
                    
                if 'freqC' in p.keys():
                    freqC_on_freqA = math.exp(p['freqC'])
                    freqG_on_freqA = math.exp(p['freqG'])
                    freqT_on_freqA = math.exp(p['freqT'])
                    freqA = 1.0/(1.0 + freqC_on_freqA + freqG_on_freqA + freqT_on_freqA)
                    freqC = freqA*freqC_on_freqA
                    freqG = freqA*freqG_on_freqA
                    freqT = freqA*freqT_on_freqA
                    print 'base freqs:',freqA,freqC,freqG,freqT
                    #u = updaters['state_freqs']
                    m = c.partition_model.getModel(0)
                    m.setStateFreqUnnorm(0, freqA)
                    m.setStateFreqUnnorm(1, freqC)
                    m.setStateFreqUnnorm(2, freqG)
                    m.setStateFreqUnnorm(3, freqT)

                if 'gamma_shape' in p.keys():
                    shape = math.exp(p['gamma_shape'])
                    print 'previous gamma shape =',m.getShape()
                    
                    print 'gamma_shape:',shape
                    #u = updaters['gamma_shape']
                    m.setShape(shape)
                    
                # temporary!
                print 'm.getModelName() =',m.getModelName()
                #for d in dir(c.partition_model):
                #    print d
                nsubsets = c.partition_model.getNumSubsets()
                print 'c.partition_model.getNumSubsets() =',nsubsets
                for s in range(nsubsets):
                    print 'c.partition_model.getSubsetRelRate(%d) =' % s,c.partition_model.getSubsetRelRate(s)
                    m = c.partition_model.getModel(s)
                    #for d in dir(m):
                    #    print d
                    print 'm.getStateFreqs()  =',m.getStateFreqs()
                    print 'm.getRelRates()    =',m.getRelRates()
                    print 'm.getShape()       =',m.getShape()
                    print 'm.getNGammaRates() =',m.getNGammaRates()
                raw_input('doof')

                print 'edge lengths:'
                for ss in e.keys():
                    edgelen = math.exp(e[ss])
                    print ss,edgelen
                                        
                # set edge lengths
                tree = c.getTree()
                #print 'tree:',tree.newick
                ntips = tree.getNObservables()
                #print 'ntips=',ntips
                tree.recalcAllSplits(ntips)
                
                nd = tree.getFirstPreorder()
                assert nd.isRoot(), 'The first preorder node should be the root'
                treelen = 0.0
                assert tree.hasEdgeLens(), 'tree has no edge lengths'
                while True:
                    nd = nd.getNextPreorder()
                    if not nd:
                        break
                    else:
                        # Grab the split and invert it if necessary to attain a standard polarity
                        was_inverted = False
                        s = nd.getSplit()
                        if s.isBitSet(0):
                            was_inverted = True
                            s.invertSplit()
                            
                        # Create a string representation of the split
                        ss = s.createPatternRepresentation()
                        
                        #print '------------------------------------------------------------'
                        #print 'nd.getNodeName()   =',nd.getNodeName()
                        #print 'nd.getNodeNumber() =',nd.getNodeNumber()
                        #print 'nd.getEdgeLen()    =',nd.getEdgeLen()
                        #print 'nd.getParent()     =',(nd.getParent() is None) and 'None' or nd.getParent().getNodeNumber()
                        #print 'nd.getLeftChild()  =',(nd.getLeftChild() is None) and 'None' or nd.getLeftChild().getNodeNumber()
                        #print 'nd.getRightSib()   =',(nd.getRightSib() is None) and 'None' or nd.getRightSib().getNodeNumber()
                        #print 'nd.isRoot()     =',nd.isRoot()
                        #print 'nd.isTip()      =',nd.isTip()
                        #print 'nd.isInternal() =',nd.isInternal()
                        #print 'was_inverted    =',was_inverted
                        #print 'ss              =',ss
                        
                        log_edge_len = e[ss]
                        edge_len = math.exp(log_edge_len)
                        print '%s edge_len = %g' % (ss,edge_len)
                        nd.setEdgeLen(edge_len)  
                print 'tree =',tree.newick
                                    
                c.likelihood.replaceModel(c.partition_model)
                c.chain_manager.refreshLastLnLike()
                self.output('log-likelihood = %s' % c.chain_manager.getLastLnLike())
                c.chain_manager.refreshLastLnPrior()
                self.output('log-prior = %s' % c.chain_manager.getLastLnPrior())
                raw_input('Press return to continue...')
            
    def summarize(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Summarizes what was found - mostly for debugging purposes.
        
        """
        self.output('\nSummary of edge length samples:')
        self.output('  Number of distinct tree topologies: %d' % len(self.edge_lengths.keys()))
        self.output('  Number of samples within each distinct tree topology:')
        self.output('  %12s\t%6s' % ('tree', 'samples'))
        for i,treeid in enumerate(self.edge_lengths.keys()):
            self.output('  %12d\t%6d' % (i+1,len(self.edge_lengths[treeid])))
            
        self.output('\nSummary of parameter samples:')
        self.output('  Number of distinct tree topologies: %d' % len(self.parameters.keys()))
        self.output('  Number of samples within each distinct tree topology:')
        self.output('  %12s\t%6s' % ('tree', 'samples'))
        for i,treeid in enumerate(self.parameters.keys()):
            self.output('  %12d\t%6d' % (i+1,len(self.parameters[treeid])))
            
    def marginal_likelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This is the main member function of the class. It reads a sample of 
        trees and edge lengths from the specified tree file and a 
        corresponding sample of parameters from the specified params file.
        From these, it estimates the marginal likelihood using the IDR method.
        
        """
        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        self.stdout.phycassert(input_trees is not None, 'trees cannot be None when the idr method is called')
        self.stdout.phycassert(input_trees.__class__.__name__ == 'TreeCollection', 'trees must be a TreeCollection object')
        
        # Check to make sure user specified an input params file
        input_params = self.opts.params
        self.stdout.phycassert(input_params is not None and len(input_params) > 0, 'params cannot be None or empty when the idr method is called')
        
        # Store transformed edge lengths and parameter values in dictionaries self.edge_lengths and self.parameters, respectively
        self.harvestTreesAndParams(input_trees, input_params)
        
        # Store transformed parameter values in dictionary self.parameters
        #self.harvestParameters(input_params)
                        
        # Store log-transformed edge lengths in dictionary self.edge_lengths
        #self.harvestEdgeLengths(input_trees)
        
        # Read in the data and store in self.data_matrix
        self.loadData()

        # These are debugging functions whose reports tell us if everything is working as expected
        self.checkPosterior()
        self.summarize()
        
