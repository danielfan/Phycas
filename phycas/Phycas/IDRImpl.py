import os,sys,math,random
from phycas import *
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
        self.stored_tree_defs = None    # list of tree definitions from the trees file
        self.param_file_lines = None    # list of lines from the params file (header lines excluded)
        self.starting_tree = None       # the tree to be processed
        self.data_matrix = None         # the data matrix itself
        self.taxon_labels = None        # list of taxon labels from the data matrix
        self.ntax = None                # number of taxa in the data matrix
        self.nchar = None               # number of sites in the data matrix
        self.param_names = None         # list of the names of the parameters
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
        
    def storeTrees(self, input_trees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Open specified tree file and read trees therein, storing the list of
        tree definitions (strings) in self.stored_tree_defs.
        
        """
        burnin = self.opts.burnin
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_tree_defs = list(input_trees)
        self.taxon_labels = input_trees.taxon_labels # this line must follow the coercion of the trees to a list (in case that is what triggers the readinf of the file with the taxon labels)
        num_stored_trees = len(self.stored_tree_defs)
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

    def fillEdgeLenDict(self, tree_def):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build and return a tuple containing: (1) a dictionary containing the 
        string representations (keys) and log-transformed edge lengths 
        (values) for the tree definition supplied in tree_def; and (2) the
        tree id, which is a tuple of split representations for this tree that
        can be used to uniquely identify the tree topology.
        
        """
        edgelen_dict = {}
        
        # Build the tree
        tree = Phylogeny.Tree()
        tree_def.buildTree(tree)
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

                # Create a dictionary entry with ss as key and edge_len as value
                edgelen_dict[ss] = edge_len
        return edgelen_dict, tuple(split_list)
        
    def harvestTreesAndParams(self, input_trees, input_params):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build dictionaries self.parameters and self.edge_lengths from 
        supplied input_trees and input_params.
        
        """
        self.storeTrees(input_trees)
        self.storeParams(input_params)
        
        # Initialize the data member that will hold all information about edge lengths
        self.edge_lengths = {}
        
        # Initialize the data member that will hold all information about parameter samples
        self.parameters = {}
        
        self.sample_size = 0
        post_burnin_paramvects = self.param_file_lines[self.opts.burnin+3:] # eliminate burnin samples as well as 2 header lines and starting state (generation 0)
        post_burnin_treedefs = self.stored_tree_defs[self.opts.burnin+1:]   # eliminate burnin samples as well as starting state (generation 0)
        for tree_def, param_vect in zip(post_burnin_treedefs,post_burnin_paramvects): 
            self.sample_size += 1

            # This dictionary will be filled with parameter values: e.g. sample_dict['freqA'] = 0.24335
            param_dict = self.fillParamDict(param_vect)
            
            # This dictionary will be filled with log-transformed edge lengths: e.g. edgelen_dict[s] = -1.10293
            # where the split object s associated with the edge is used as the key
            edgelen_dict, tree_id = self.fillEdgeLenDict(tree_def)
            
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

    #def checkLikelihoods(self):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Calculates log-likelihoods from the supplied parameters and 
    #    edge_lengths lists and compares with lnL stored in parameters.
    #    
    #    """
    #    pass
    #    #assert len(parameters) == len(edge_lengths), 'parameters (%d) and edge_lengths (%d) lists not the same length' % (len(parameters), len(edge_lengths))
    #    #assert len(parameters) == len(self.stored_tree_defs[self.opts.burnin:]), 'parameters and stored_tree_defs lists not the same length'
    #    #for p,e,d in zip(parameters, edge_lengths, self.stored_tree_defs[self.opts.burnin:]):
    #    #    self.starting_tree = Phylogeny.Tree()
    #    #    d.buildTree(self.starting_tree)
    #    #    ntips = self.starting_tree.getNObservables()
    #    #    self.starting_tree.recalcAllSplits(ntips)
    #    #    core = LikelihoodCore(self)
    #    #    core.setupCore()
    #    #    core.prepareForLikelihood()
    #    #    lnL = core.calcLnLikelihood()
    #    #    self.stdout.info('lnL = %g (%g)' % (lnL,p['lnL']))
            
    def summarize(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Summarizes what was found - mostly for debugging purposes.
        
        """
        self.stdout.info('\nSummary of edge length samples:')
        self.stdout.info('  Number of distinct tree topologies: %d' % len(self.edge_lengths.keys()))
        self.stdout.info('  Number of samples within each distinct tree topology:')
        self.stdout.info('  %12s\t%6s' % ('tree', 'samples'))
        for i,treeid in enumerate(self.edge_lengths.keys()):
            self.stdout.info('  %12d\t%6d' % (i+1,len(self.edge_lengths[treeid])))
            
        self.stdout.info('\nSummary of parameter samples:')
        self.stdout.info('  Number of distinct tree topologies: %d' % len(self.parameters.keys()))
        self.stdout.info('  Number of samples within each distinct tree topology:')
        self.stdout.info('  %12s\t%6s' % ('tree', 'samples'))
        for i,treeid in enumerate(self.parameters.keys()):
            self.stdout.info('  %12d\t%6d' % (i+1,len(self.parameters[treeid])))

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
        
        self.summarize()
        
