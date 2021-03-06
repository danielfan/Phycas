import os,sys,math
from phycas import *
import phycas.Phylogeny as Phylogeny
import phycas.ProbDist as ProbDist
import phycas.Likelihood as Likelihood
from LikelihoodCore import LikelihoodCore

#from threading import Thread
class UnimapSpreadingWrapper(object):
    def __init__(self, mcmc):
        self.unimap_ls_move_list = [Likelihood.UnimapLSMove(mcmc.likelihood) for i in range(mcmc.parent.opts.unimap_thread_count)]
        self.unimap_spreader_move = Likelihood.UnimapTopoMoveSpreader()
        self.unimap_spreader_move.setName("unimap_topo_move_spreader")
        self.unimap_spreader_move.setWeight(mcmc.parent.opts.unimap_ls_move_weight)
        self.unimap_spreader_move.setTree(mcmc.tree)
        self.unimap_spreader_move.setLot(mcmc.r)
        for n, m in enumerate(self.unimap_ls_move_list):
            self.unimap_spreader_move.addTopoMoveToSpreader(m)
            m.setName("unimap_LS_move %d" % n)
            m.setWeight(0)
            m.setTree(mcmc.tree)
            m.setModel(mcmc.partition_model.getModel(0))
            m.setLot(mcmc.r)
            mcmc.chain_manager.addMove(m)

    def setSaveDebugInfo(self, f):
        self.unimap_spreader_move.setSaveDebugInfo(f)
        for m in self.unimap_ls_move_list:
            m.setSaveDebugInfo(f)
    def getName(self):
        return self.unimap_spreader_move.getName()
    def getDebugInfo(self):
        return self.unimap_spreader_move.getDebugInfo()
    def update(self):
        #import pdb; pdb.set_trace()
        self.unimap_spreader_move.update()
        return True
        # below here is old code that failed to be truly parallel due to the Python GIL
        #threads = [Thread(target=i.update) for i in self.unimap_ls_move_list]
        #for t in threads:
        #   t.start()
        #for t in threads:
        #   t.join()
        #print "all threads done!"
        #return True
    
class MarkovChain(LikelihoodCore):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The MarkovChain class encapsulates the notion of a Markov chain used
    in Bayesian analyses. The MarkovChain class has the ability to
    orchestrate a Markov chain Monte Carlo (MCMC) analysis. In 
    Metropolis-coupled MCMC, the cold chain and each of the heated chains
    are MarkovChain objects. The primary addition MarkovChain adds to the
    base class LikelihoodCore are prior distributions and the 
    self.heating_power data member.
    
    """
    def __init__(self, parent, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The MarkovChain constructor makes a copy of the supplied parent
        object, clones all of the prior ProbabilityDistribution objects in the
        supplied parent object, and sets self.heating_power to the supplied
        power.
        
        """
        LikelihoodCore.__init__(self, parent)
        
        self.parent                     = parent    # Note: self.parent is the MCMCImpl object
        self.boldness                   = 0.0
        self.heating_power              = power
        self.chain_manager              = None
        self.tree_scaler_move           = None      
        self.subset_relrates_move       = None
        self.edge_move                  = None
        self.unimap_fast_nni_move       = None
        self.unimap_sample_ambig_move   = None
        self.unimap_nni_move            = None
        self.unimap_node_slide_move     = None
        self.unimapping_move            = None
        self.unimap_edge_move           = None
        self.larget_simon_move          = None
        # FLEXCAT_MOVE
        #self.ncat_move                 = None
        self.bush_move                  = None
        self.topo_prior_calculator      = None
        
        self.state_freq_moves           = []
        self.rel_rate_moves             = []
        self.all_updaters_list = None

        self.setupChain()

    def resetNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the resetNumLikelihoodEvals function of self.likelihood. This
        resets the number of likelihood evaluations performed to 0.

        """
        # Note: likelihood data member inherited from LikelihoodCore
        return self.likelihood.resetNumLikelihoodEvals()
    
    def getNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the getNumLikelihoodEvals function of self.likelihood. This 
        returns the number of likelihood evaluations performed since the last
        call of resetNumLikelihoodEvals.

        """
        # Note: likelihood data member inherited from LikelihoodCore
        return self.likelihood.getNumLikelihoodEvals()
        
    def paramFileHeader(self, paramf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the parameter file. The
        parameter paramf should be an (already open for writing) file object.
        The parameter file corresponds to the *.p file produced by MrBayes. It
        begins with a line containing the initial random number seed, followed
        by a line of column titles for the (tab-delimited) sampled parameter
        information on the subsequent lines.

        """
        # Note: r data member inherited from LikelihoodCore
        paramf.write('[ID: %d]\n' % self.r.getInitSeed())
        if self.parent.opts.doing_steppingstone_sampling:
            paramf.write('Gen\tbeta\tlnL\tlnPrior')
        else:
            paramf.write('Gen\tlnL\tlnPrior')

        if self.parent.opts.doing_steppingstone_sampling and not self.parent.opts.ssobj.ti:
            paramf.write('\tlnRefDens')

        # If the user has defined a reference tree, add a column for the Robinson-Foulds
        # distance between the sampled tree and the reference tree
        if self.parent.ref_tree is not None:
            paramf.write('\tdRF')
                        
        # If using a model that allows polytomies, include a column indicating the 
        # resolution class of the tree
        if self.parent.opts.allow_polytomies:
            paramf.write('\tResClass')
                        
        paramf.write('\tTL')
        if self.parent.opts.fix_topology:
            nbrlens = self.tree.getNNodes() - 1
            for i in range(nbrlens):
                paramf.write('\tbrlen%d' % (i+1))
            self.parent.output('\nKey to the edges (preorder traversal):\n%s' % self.tree.keyToEdges())
        
        nmodels = self.partition_model.getNumSubsets()
        if nmodels > 1:
            for i in range(nmodels):
                paramf.write('\tm_%d' % (i+1,))
        for i in range(nmodels):
            m = self.partition_model.getModel(i)
            paramf.write(m.paramHeader())
            
    def treeFileHeader(self, treef):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the tree file. The
        parameter treef should be an (already open for writing) file object.
        The tree file corresponds to the *.t file produced by MrBayes, and is
        in NEXUS format. This function writes the opening "#NEXUS" line, a 
        comment containing the initial random number seed, and the beginning
        of a NEXUS TREES block, including the TRANSLATE command.

        """
        treef.write('#NEXUS\n')
        treef.write('[ID: %d]\n' % self.r.getInitSeed())
        treef.write('begin trees;\n')
        if self.parent.ntax > 0:
            treef.write('\ttranslate\n')
            for i in range(self.parent.ntax):
                if self.parent.taxon_labels[i].find(' ') < 0:
                    # no spaces found in name
                    treef.write('\t\t%d %s%s\n' % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))
                else:
                    # at least one space in taxon name, so enclose name in quotes
                    treef.write("\t\t%d '%s'%s\n" % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))

    def setPower(self, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the heating_power data member and calls the setPower method for
        every updater so that all updaters are immediately informed that the
        power for this chain has changed.
        
        """
        self.heating_power = power
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(power)

    def setBoldness(self, boldness):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the boldness data member and calls the setBoldness method for
        every updater so that all updaters are immediately informed that the
        boldness for this chain has changed. The boldness is a value from 0 to
        100 that specifies the boldness of Metropolis-Hastings moves. Setting
        the boldness has no effect on slice sampling based updaters. Each 
        move class defines what is meant by boldness. The boldness is changed
        during an MCMC run in some circumstances, such as during a path 
        sampling analysis where the target distribution changes during the 
        run.
        
        """
        self.boldness = boldness
        for updater in self.chain_manager.getAllUpdaters():
            updater.setBoldness(boldness)

    def setupChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The setupChain method prepares a MarkovChain object before an MCMC
        analysis is started. This method is called at the end of the 
        MarkovChain constructor. The duties performed include:
        1) calling the base class setupCore method; 
        2) preparing the starting tree by adding data structures needed for
        computing likelihood (transition matrices and conditional likelihood
        arrays); 
        3) setting up prior distributions for model parameters; 
        4) creating an MCMCManager object and adding all relevant updaters to
        it; 
        and 5) makes sure each updater knows the type of heating and the 
        power.
        
        """
        LikelihoodCore.setupCore(self)

        from phycas import partition,model
        if self.parent.opts.partition.noData():
            self.likelihood.setNoData()
            
        LikelihoodCore.prepareForLikelihood(self)

        self.python_only_moves = []

        # add priors to models already added (by LikelihoodCore.setupCore) to partition_model
        modelspecs = partition.getModels()  # partition here refers to the global object (associated with the Phycas partition command)
        print 'modelspecs has length %d:' % len(modelspecs)
        nmodels = self.partition_model.getNumSubsets()
        if nmodels == 1:
            print 'partition_model contains 1 model (i.e. unpartitioned)'
        else:
            print 'partition_model contains %d models' % nmodels
        self.chain_manager = Likelihood.MCMCChainManager()
        for i in range(nmodels):
            # get the Model (as defined in likelihood_models.cpp) associated with partition subset i
            m = self.partition_model.getModel(i)

            # get the model specification (defined in Model.py) stored in (Python) partition object
            mspec = modelspecs[i]

            #implemented = not (self.parent.opts.fix_topology and mspec.edgelen_hyperprior is not None)
            #self.parent.phycassert(implemented, 'Cannot currently specify an edge length hyperprior and fix the topology at the same time')
            
            # Copy priors related to edge lengths
            # Note: while these priors are copied for every subset, only those for the
            # first subset are actually used (at this writing, 31 Jan 2010) because both
            # tree topology and edge lengths are always (at this writing) linked across subsets
            
            # POLPY_NEWWAY  // no touch
            separate_edge_len_dists = mspec.separate_edgelen_hyper
            # else
            #separate_edge_len_dists = mspec.internal_edgelen_prior is not mspec.external_edgelen_prior
            
            #raw_input('separate_edge_len_dists = %s' % (separate_edge_len_dists and 'yes' or 'no'))
            m.separateInternalExternalEdgeLenPriors(separate_edge_len_dists)
            if mspec.external_edgelen_prior is not None:
                m.setExternalEdgeLenPrior(mspec.external_edgelen_prior.cloneAndSetLot(self.r))
            if mspec.internal_edgelen_prior is not None:
                m.setInternalEdgeLenPrior(mspec.internal_edgelen_prior.cloneAndSetLot(self.r))
            if mspec.edgelen_hyperprior is not None:
                m.setEdgeLenHyperPrior(mspec.edgelen_hyperprior.cloneAndSetLot(self.r))
                if mspec.fix_edgelen_hyperparam:
                    m.fixEdgeLenHyperprior()    #@POL should be named fixEdgeLenHyperparam
            else:
                m.setEdgeLenHyperPrior(None)
                
            # Copy priors for this model's parameters
            if mspec.type == 'codon':
                m.setKappaPrior(mspec.kappa_prior.cloneAndSetLot(self.r))
                m.setOmegaPrior(mspec.omega_prior.cloneAndSetLot(self.r))
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
            elif mspec.type == 'gtr':
                if mspec.update_relrates_separately:
                    m.setRelRateParamPrior(mspec.relrate_param_prior.cloneAndSetLot(self.r))
                else:
                    self.parent.phycassert(mspec.relrate_prior.getDistName() == 'Dirichlet', 'mspec.relrate_prior must be of type Dirichlet')
                    m.setRelRatePrior(mspec.relrate_prior)
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    self.parent.phycassert(mspec.state_freq_prior.getDistName() == 'Dirichlet', 'mspec.state_freq_prior must be of type Dirichlet')
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
            elif mspec.type == 'hky':
                m.setKappaPrior(mspec.kappa_prior.cloneAndSetLot(self.r))
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))

            # Copy priors related to among-site rate heterogeneity
            if mspec.num_rates > 1:
                m.setDiscreteGammaShapePrior(mspec.gamma_shape_prior.cloneAndSetLot(self.r))
            if mspec.pinvar_model:
                m.setPinvarPrior(mspec.pinvar_prior.cloneAndSetLot(self.r))
            
            # If user specifies fixed tree topology, make each edge length a separate parameter; otherwise, use edge length master parameters
            if self.parent.opts.fix_topology:
                m.setEdgeSpecificParams(True)
            else:
                m.setEdgeSpecificParams(False)
                
            # Add all necessary updaters to the MCMCManager
            if nmodels == 1:
                subset_pos = -1
            else:
                subset_pos = i
                
            self.chain_manager.addMCMCUpdaters(
                m,                                  # substitution model
                self.tree,                          # tree
                self.likelihood,                    # likelihood calculation machinery
                self.r,                             # pseudorandom number generator
                #False,                             # separate_edgelen_params (separate_edgelen_params is now True only if topology is fixed)
                self.parent.opts.slice_max_units,   # maximum number of slice units allowed
                self.parent.opts.slice_weight,      # weight for each parameter added
                subset_pos)                         # i is the subset (needed so that add edge length params will only be added for for first subset)
                
            if self.parent.opts.doing_steppingstone_sampling:
                self.chain_manager.setMinSSWPSampleSize(self.parent.opts.ssobj.minsample)
                
            # Add subset-model-specific moves
            if not mspec.type == 'jc' and not mspec.update_freqs_separately:
                # Create a StateFreqMove to update entire state frequency vector
                sfm = Likelihood.StateFreqMove()
                if mspec.type == 'codon': 
                    sfm.setDimension(61)
                else:
                    sfm.setDimension(4)
                if nmodels > 1:
                    sfm.setName("state_freqs_%d" % (i+1,))
                else:
                    sfm.setName("state_freqs")
                sfm.setWeight(self.parent.opts.state_freq_weight)
                sfm.setPosteriorTuningParam(self.parent.opts.state_freq_psi)
                sfm.setPriorTuningParam(self.parent.opts.state_freq_psi0)
                sfm.setTree(self.tree)
                sfm.setModel(m)
                sfm.setTreeLikelihood(self.likelihood)
                sfm.setLot(self.r)
                if m.stateFreqsFixed():
                    sfm.fixParameter()
                sfm.setMultivarPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
                self.chain_manager.addMove(sfm)
                self.state_freq_moves.append(sfm)
            
            if mspec.type == 'gtr' and not mspec.update_relrates_separately:
                # Create a RelRateMove to update entire relative rates vector
                rrm = Likelihood.RelRatesMove()
                if nmodels > 1:
                    rrm.setName("relrates_%d" % (i+1,))
                else:
                    rrm.setName("relrates")
                rrm.setWeight(self.parent.opts.rel_rate_weight)
                rrm.setPosteriorTuningParam(self.parent.opts.rel_rate_psi)
                rrm.setPriorTuningParam(self.parent.opts.rel_rate_psi0)
                rrm.setTree(self.tree)
                rrm.setModel(m)
                rrm.setTreeLikelihood(self.likelihood)
                rrm.setLot(self.r)
                #if self.model.relRatesFixed():
                #    rrm.fixParameter()
                rrm.setMultivarPrior(mspec.relrate_prior.cloneAndSetLot(self.r))
                self.chain_manager.addMove(rrm)
                self.rel_rate_moves.append(rrm)

        self.likelihood.replaceModel(self.partition_model)          
                
        if self.parent.opts.data_source is None:
            self.likelihood.setNoData() # user apparently wants to run MCMC with no data
            
        model0 = self.partition_model.getModel(0)
        
        # Create a TreeScalerMove object to handle scaling the entire tree to allow faster
        # convergence in edge lengths. This move is unusual in using slice sampling rather
        # than Metropolis-Hastings updates: most "moves" in parent are Metropolis-Hastings.
        if self.parent.opts.tree_scaler_weight > 0:
            self.tree_scaler_move = Likelihood.TreeScalerMove()
            self.tree_scaler_move.setName("tree_scaler")
            self.tree_scaler_move.setWeight(self.parent.opts.tree_scaler_weight)
            self.tree_scaler_move.setPosteriorTuningParam(self.parent.opts.tree_scaler_lambda)
            self.tree_scaler_move.setPriorTuningParam(self.parent.opts.tree_scaler_lambda0)
            self.tree_scaler_move.setTree(self.tree)
            self.tree_scaler_move.setModel(model0)
            self.tree_scaler_move.setTreeLikelihood(self.likelihood)
            self.tree_scaler_move.setLot(self.r)
            if model0.edgeLengthsFixed():
                self.tree_scaler_move.fixParameter()
            self.chain_manager.addMove(self.tree_scaler_move)

        # If more than one partition subset, add a SubsetRelRate move to modify the 
        # vector of relative substitution rates for each subset
        if (nmodels > 1):
            self.subset_relrates_move = Likelihood.SubsetRelRatesMove()
            self.subset_relrates_move.setDimension(nmodels)
            self.subset_relrates_move.setName("subset_relrates")
            self.subset_relrates_move.setWeight(self.parent.opts.subset_relrates_weight)
            self.subset_relrates_move.setPosteriorTuningParam(self.parent.opts.subset_relrates_psi)
            self.subset_relrates_move.setPriorTuningParam(self.parent.opts.subset_relrates_psi0)
            self.subset_relrates_move.setTree(self.tree)
            self.subset_relrates_move.setModel(None)    # the model data member is ignored in this case; instead, the partition model stores the parameters
            self.subset_relrates_move.setPartitionModel(self.partition_model)
            self.subset_relrates_move.setTreeLikelihood(self.likelihood)
            self.subset_relrates_move.setLot(self.r)
            subset_proportions = partition.getSubsetProportions()
            self.subset_relrates_move.setSubsetProportions(subset_proportions)
            if partition.fix_subset_relrates:
                self.subset_relrates_move.fixParameter()
            else:
                self.subset_relrates_move.freeParameter()
                # only assign a prior distribution if subset relative rates are not fixed
                if partition.subset_relrates_prior is None:
                    param_list = tuple([1.0]*nmodels)
                    d = ProbDist.RelativeRateDistribution(param_list)
                    d.setCoefficients(subset_proportions)
                    self.subset_relrates_move.setMultivarPrior(d.cloneAndSetLot(self.r))
                    self.partition_model.setSubsetRelRatePrior(d.cloneAndSetLot(self.r))
                else:
                    self.parent.phycassert(partition.subset_relrates_prior.getDistName() == 'RelativeRateDistribution', 'partition.subset_relrates_prior must be of type RelativeRateDistribution')
                    self.parent.phycassert(partition.subset_relrates_prior.getNParams() == nmodels, 'partition.subset_relrates_prior has dimension %d, but there are %d subsets in the partition. Try setting partion.subset_relrates_prior = None to get default flat Dirichlet prior of the appropriate dimension' % (partition.subset_relrates_prior.getNParams(), nmodels))
                    partition.subset_relrates_prior.setCoefficients(subset_proportions)
                    self.subset_relrates_move.setMultivarPrior(partition.subset_relrates_prior.cloneAndSetLot(self.r))
                    self.partition_model.setSubsetRelRatePrior(partition.subset_relrates_prior.cloneAndSetLot(self.r))
            self.chain_manager.addMove(self.subset_relrates_move)
        
        #OLDWAY
        #if self.parent.opts.fix_topology:
        #   # Create an EdgeMove object to handle Metropolis-Hastings
        #   # updates to the edge lengths only (does not change the topology)
        #   self.edge_move = Likelihood.EdgeMove()
        #   self.edge_move.setName("edge_move")
        #   self.edge_move.setWeight(self.parent.opts.edge_move_weight)
        #   self.edge_move.setPosteriorTuningParam(self.parent.opts.edge_move_lambda)
        #   self.edge_move.setPriorTuningParam(self.parent.opts.edge_move_lambda0)
        #   self.edge_move.setTree(self.tree)
        #   self.edge_move.setModel(model0)
        #   self.edge_move.setTreeLikelihood(self.likelihood)
        #   self.edge_move.setLot(self.r)
        #   self.edge_move.setLambda(self.parent.opts.edge_move_lambda)
        #   if model0.edgeLengthsFixed():
        #       self.edge_move.fixParameter()
        #   self.chain_manager.addMove(self.edge_move)

        if self.parent.opts.use_unimap:
            if False:
                # Create a UnimapFastNNIMove (replaces LargetSimonMove for unimap analyses)
                self.unimap_fast_nni_move = Likelihood.UnimapFastNNIMove()
                self.unimap_fast_nni_move.setName("unimap_fastNNI_move")
                self.unimap_fast_nni_move.setWeight(self.parent.opts.unimap_fast_nni_move_weight)
                self.unimap_fast_nni_move.setTree(self.tree)
                self.unimap_fast_nni_move.setModel(model0)
                self.unimap_fast_nni_move.setTreeLikelihood(self.likelihood)
                self.unimap_fast_nni_move.setLot(self.r)
                self.chain_manager.addMove(self.unimap_fast_nni_move)

            # Create a UnimapSampleAmbigMove
            wt = self.parent.opts.unimap_sample_ambig_move_weight
            self.unimap_sample_ambig_move = Likelihood.UnimapSampleAmbigMove(self.likelihood, self.tree, wt)
            self.unimap_sample_ambig_move.setName("unimap_sample_ambig_move")
            num_ambig = self.unimap_sample_ambig_move.getNumAmbigNodes()
            if num_ambig > 0:
                if wt > 0.0:
                    self.unimap_sample_ambig_move.setLot(self.r)
                    self.chain_manager.addMove(self.unimap_sample_ambig_move)
                else:
                    self.parent.phycassert(False, "unimap_sample_ambig_move_weight was set to 0, but %d ambiguous leaves were found" % num_ambig)

            # Create a UnimapLSMove (replaces LargetSimonMove for unimap analyses)
            
            if self.parent.opts.unimap_thread_count > 1 and (self.parent.opts.unimap_ls_move_weight > 0):
                self.unimap_spreader_move = UnimapSpreadingWrapper(self)
                self.python_only_moves.append((self.unimap_spreader_move, self.parent.opts.unimap_ls_move_weight))
                #self.chain_manager.addMove(self.unimap_spreader_move)
            else:
                # Create a UnimapLSMove (replaces LargetSimonMove for unimap analyses)
                self.unimap_ls_move = Likelihood.UnimapLSMove(self.likelihood)
                self.unimap_ls_move.setName("unimap_LS_move")
                self.unimap_ls_move.setWeight(self.parent.opts.unimap_ls_move_weight)
                self.unimap_ls_move.setTree(self.tree)
                self.unimap_ls_move.setModel(model0)
                self.unimap_ls_move.setLot(self.r)
                self.chain_manager.addMove(self.unimap_ls_move)

            # Create a UnimapNNIMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_nni_move = Likelihood.UnimapNNIMove(self.likelihood)
            self.unimap_nni_move.setName("unimap_NNI_move")
            self.unimap_nni_move.setWeight(self.parent.opts.unimap_nni_move_weight)
            self.unimap_nni_move.setTree(self.tree)
            self.unimap_nni_move.setModel(model0)
            self.unimap_nni_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_nni_move)


            # Create a UnimapNodeSlideMove (replaces LargetSimonMove for unimap analyses)
            self.unimap_node_slide_move = Likelihood.UnimapNodeSlideMove(self.likelihood)
            self.unimap_node_slide_move.setName("unimap_nodeslide_move")
            self.unimap_node_slide_move.setWeight(self.parent.opts.unimap_node_slide_move_weight)
            self.unimap_node_slide_move.setTree(self.tree)
            self.unimap_node_slide_move.setModel(model0)
            self.unimap_node_slide_move.setLot(self.r)
            self.chain_manager.addMove(self.unimap_node_slide_move)

            # Create a MappingMove (to refresh the mapping for all sites)
            self.unimapping_move = Likelihood.MappingMove()
            self.unimapping_move.setName("univent_mapping_move")
            self.unimapping_move.setWeight(self.parent.opts.mapping_move_weight)
            self.unimapping_move.setTree(self.tree)
            self.unimapping_move.setModel(model0)
            self.unimapping_move.setTreeLikelihood(self.likelihood)
            self.unimapping_move.setLot(self.r)
            self.chain_manager.addMove(self.unimapping_move)
            
            # Create a UnimapEdgeMove object to handle Metropolis-Hastings
            # updates to the edge lengths only (does not change the topology)
            self.unimap_edge_move = Likelihood.UnimapEdgeMove()
            self.unimap_edge_move.setName("unimap_edge_length_move")
            self.unimap_edge_move.setWeight(self.parent.opts.unimap_edge_move_weight)
            self.unimap_edge_move.setPosteriorTuningParam(self.parent.opts.unimap_edge_move_lambda)
            self.unimap_edge_move.setPriorTuningParam(self.parent.opts.unimap_edge_move_lambda0)
            self.unimap_edge_move.setTree(self.tree)
            self.unimap_edge_move.setModel(model0)
            self.unimap_edge_move.setTreeLikelihood(self.likelihood)
            self.unimap_edge_move.setLot(self.r)
            self.unimap_edge_move.setLambda(self.parent.opts.unimap_edge_move_lambda)
            if model0.edgeLengthsFixed():
                self.unimap_edge_move.fixParameter()
            self.chain_manager.addMove(self.unimap_edge_move)

        elif not self.parent.opts.fix_topology:
            # Create a LargetSimonMove object to handle Metropolis-Hastings
            # updates to the tree topology and edge lengths
            self.larget_simon_move = Likelihood.LargetSimonMove()
            self.larget_simon_move.setName("larget_simon_local")
            self.larget_simon_move.setWeight(self.parent.opts.ls_move_weight)
            self.larget_simon_move.setPosteriorTuningParam(self.parent.opts.ls_move_lambda)
            self.larget_simon_move.setPriorTuningParam(self.parent.opts.ls_move_lambda0)
            self.larget_simon_move.setTree(self.tree)
            self.larget_simon_move.setModel(model0)
            self.larget_simon_move.setTreeLikelihood(self.likelihood)
            self.larget_simon_move.setLot(self.r)
            self.larget_simon_move.setLambda(self.parent.opts.ls_move_lambda)
            if model0.edgeLengthsFixed():
                self.larget_simon_move.fixParameter()
            self.chain_manager.addMove(self.larget_simon_move)

        # If requested, create a BushMove object to allow polytomous trees
        if self.parent.opts.allow_polytomies:
            # Create a BushMove object
            self.bush_move = Likelihood.BushMove()

            # Set up the topology prior
            self.topo_prior_calculator = self.bush_move.getPolytomyTopoPriorCalculator()
            self.topo_prior_calculator.chooseUnrooted()
            self.topo_prior_calculator.setC(self.parent.opts.topo_prior_C)
            if self.parent.opts.polytomy_prior:
                self.topo_prior_calculator.choosePolytomyPrior()
            else:
                self.topo_prior_calculator.chooseResolutionClassPrior()
                
            # Continue setting up BushMove object
            self.bush_move.setName("bush_move")
            self.bush_move.setWeight(self.parent.opts.bush_move_weight)
            self.bush_move.setTree(self.tree)
            self.bush_move.setModel(model0)
            self.bush_move.setTreeLikelihood(self.likelihood)
            self.bush_move.setLot(self.r)
            self.bush_move.setEdgeLenDistMean(self.parent.opts.bush_move_edgelen_mean)
            #self.bush_move.viewProposedMove(self.parent.bush_move_debug)
            if model0.edgeLengthsFixed():
                self.bush_move.fixParameter()
            self.bush_move.finalize()
            
            self.chain_manager.addMove(self.bush_move)

        self.chain_manager.finalize()

        if self.parent.opts.use_unimap:
            if num_ambig > 0:
                self.unimap_sample_ambig_move.sampleTipsAsDisconnected()
            self.likelihood.fullRemapping(self.tree, self.r, True) 

        # Calculate relative rates if among-site rate heterogeneity is part of the model
        # Currently, the unnormalized rates are stored in the model; the function call
        # below normalizes the rate means and probabilities so that they will be accurately
        # recorded in the first line of the parameter file
        self.likelihood.recalcRelativeRates()

        # Make sure each updater knows the heating power and heating type
        # and set boldness to 0 for each updater (0% boldness is the default,
        # with more boldness used only in steppingstone sampling where bolder
        # moves are needed as the influence of the likelihood is progressively
        # diminished)
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(self.heating_power)
            updater.setBoldness(0.0)
            if self.parent.opts.ss_heating_likelihood:
                updater.setLikelihoodHeating()
            else:
                updater.setStandardHeating()

