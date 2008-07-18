from phycas.Phycas.PhycasCommand import *
from phycas import model
from phycas.Phycas.MCMCImpl import MCMCImpl
from phycas.ProbDist import BetaDist, ExponentialDist, InverseGammaDist
class MCMC(PhycasCommand):
    def __init__(self, p):
        args = (
                ("random_seed",            0,                               "Determines the random number seed used; specify 0 to generate seed automatically from system clock", IntArgValidate(min=0)),
                ("ncycles",                10000,                           "The number of update cycles (a cycle is analogous to, but different than, a 'generation' in MrBayes; Phycas does in one cycle what MrBayes does in about 100 generations for a simple model such as JC)", IntArgValidate(min=0)),
                ("sample_every",           100,                             "The current tree topology and model parameter values will be sampled after this many cycles have elapsed since the last sample was taken", IntArgValidate(min=0)),
                ("report_every",           100,                             "A progress report will be displayed after this many cycles have elapsed since the last progress report", IntArgValidate(min=0)),
                ("verbose",                True,                            "You will get more output if True, less output if False", BoolArgValidate),
                ("quiet",                  False,                           "If True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well", BoolArgValidate),
                ("outfile_prefix",         None,                            "If None, parameter and tree files created will have a name beginning with the name of the data file; if provided, this prefix will form the first part of the parameter (e.g. <outfile_prefix>.p) and tree file (e.g. <outfile_prefix>.t) names"),
                ("model",                  model,                           "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("data_source",            'file',                          "not yet documented"),
                ("data_file_name",         '',                              "Used to specify the nexus data file name to be used for subsequent analyses"),
                ("starting_tree_source",   'random',                        "Source of the starting tree topology: can be either 'random' or 'usertree'. Note that this setting does not determine the edge lengths in the starting tree, only the topology. Starting edge lengths are determined by the probability distribution specified in starting_edgelen_dist"),
                ("tree_topology",          None,                            "Unused unless starting_tree_source is 'usertree', in which case this should be a standard newick string representation of the tree topology; e.g. '(A:0.01,B:0.071,(C:0.013,D:0.021):0.037)'"),
                ("fix_topology",           False,                           "If True, an EdgeMove move will be substituted for the LargetSimonMove, so edge lengths will be updated by slice sampling but the topology will remain unchanged during an MCMC analysis", BoolArgValidate),
                ("ls_move_lambda",         0.2,                             "The value of the tuning parameter for the Larget-Simon move", FloatArgValidate(min=0.01)),
                ("ls_move_weight",         100,                             "Larget-Simon moves will be performed this many times per cycle", IntArgValidate(min=0)),
                ("ls_move_debug",          False,                           "If set to true, TreeViewer will popup on each Larget-Simon move update showing edges affected by the proposed move", BoolArgValidate),
                ("edge_move_lambda",       0.2,                             "The value of the tuning parameter for the EdgeMove", FloatArgValidate(min=0.01)),
                ("edge_move_weight",       0,                               "Only used if fix_topology is True. Makes sense to set this to some multiple of the number of edges since each EdgeMove affects a single randomly-chosen edge ", IntArgValidate(min=0)),
                ("mapping_move_weight",     1,                              "Univent mapping will be performed this many times per cycle", IntArgValidate(min=0)),
                ("unimap_nni_move_weight", 100,                             "Unimap NNI moves will be performed this many times per cycle", IntArgValidate(min=0)),
                ("tree_scaler_weight",     0,                               "Whole-tree scaling will be performed this many times per cycle", IntArgValidate(min=0)),
                ("use_unimap",              False,                          "if True, MCMC analyses will use the uniformized mapping approach", BoolArgValidate),
                ("allow_polytomies",        False,                          "If True, do Bush moves in addition to Larget-Simon moves; if False, do Larget-Simon moves only", BoolArgValidate),
                ("polytomy_prior",          True,                           "If True, use polytomy prior; if False, use resolution class prior", BoolArgValidate),
                ("topo_prior_C",            2.0,                            "Specifies the strength of the prior (C = 1 is flat prior; C > 1 favors less resolved topologies)", FloatArgValidate(min=0.01)),
                ("bush_move_edgelen_mean",  1.0,                            "Specifies mean of exponential edge length generation distribution used by BushMove when new edges are created", FloatArgValidate(min=0.01)),
                ("bush_move_weight",        100,                            "Bush moves will be performed this many times per cycle if", IntArgValidate(min=0)),
                ("bush_move_debug",         False,                          "If set to true, TreeViewer will pop up on each Bush move update showing edges affected by the proposed move", BoolArgValidate),
                ("slice_weight",            1,                              "Slice sampled parameters will be updated this many times per cycle", IntArgValidate(min=0)),
                ("slice_max_units",         1000,                           "Max. number of units used in slice sampling", IntArgValidate(min=0)),
                ("adapt_first",             100,                            "Adaptation of slice samplers is performed the first time at cycle adapt_first. Subsequent adaptations wait twice the number of cycles as the previous adaptation. Thus, adaptation n occurs at cycle adapt_first*(2**(n - 1)). The total number of adaptations that will occur during an MCMC run is [ln(adapt_first + ncycles) - ln(adapt_first)]/ln(2)", IntArgValidate(min=0)),
                ("adapt_simple_param",      0.5,                            "Slice sampler adaptation parameter", FloatArgValidate(min=0.01)),
                ("using_hyperprior",        True,                           "not yet documented", BoolArgValidate),
                ("edgelen_hyperprior",      InverseGammaDist(2.1,1.0/1.1),  "not yet documented"),
                ("fix_edgelen_hyperparam",  False,                          "not yet documented", BoolArgValidate),
                ("starting_edgelen_hyperparam", 0.05,                       "not yet documented", FloatArgValidate(min=0.01)),
                ("internal_edgelen_dist",   None,                           "Can be used to set a prior distribution for internal edges that differs from that applied to external edges. If this is set to something besides None, you should also set external_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("external_edgelen_dist",   None,                           "Can be used to set a prior distribution for external edges that differs from that applied to internal edges. If this is set to something besides None, you should also set internal_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("edgelen_dist",            ExponentialDist(2.0),           "Sets both internal_edgelen_dist and external_edgelen_dist to the supplied value. Use this setting if you want all edges in the tree to have the same prior distribution. Using this setting will overwrite any values previously supplied for internal_edgelen_dist and external_edgelen_dist"),
                ("fix_edgelens",            False,                          "not yet documented", BoolArgValidate),
                ("starting_edgelen_dist",  ExponentialDist(10.0),           "Used to select the starting edge lengths when starting_tree_source is 'random'"),
                ("heating_lambda",         0.2,                             "not yet documented", FloatArgValidate(min=0.01)),
                ("nchains",                1,                               "The number of Markov chains to run simultaneously. One chain serves as the cold chain from which samples are drawn, the other chains are heated to varying degrees and serve to enhance mixing in the cold chain.", IntArgValidate(min=1,max=1)), # only allowing 1 chain now because multiple chains not yet fully implemented
                ("is_standard_heating",    True,                            "not yet documented", BoolArgValidate),
                ("use_flex_model",         False,                           "not yet documented", BoolArgValidate),
                ("flex_ncat_move_weight",  1,                               "Number of times each cycle to attempt an ncat move", IntArgValidate(min=0)),
                ("flex_num_spacers",       1,                               "Number of fake rates between each adjacent pair of real rates", IntArgValidate(min=1)),
                ("flex_phi",               0.25,                            "Proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)", ProbArgValidate()),
                ("flex_L",                 1.0,                             "Upper bound of interval used for unnormalized relative rate parameter values", FloatArgValidate(min=0.01)),
                ("flex_lambda",            1.0,                             "Parameter of Poisson prior on the number of extra categories", FloatArgValidate(min=0.01)),
                ("flex_prob_param_prior",  ExponentialDist(1.0),            "not yet documented"),
                ("uf_num_edges",           50,                              "Number of edges to traverse before taking action to prevent underflow", IntArgValidate(min=1)),
                ("ntax",                   0,                               "To explore the prior, set to some positive value. Also set data_source to None", IntArgValidate(min=0)),
                )
        PhycasCommand.__init__(self, p, args, "mcmc", "The mcmc command is used to conduct a Bayesian Markov chain Monte Carlo analysis.")

        # The data members added below should be hidden from the user because they are for use by phycas developers.
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["debugging"] = False
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        mcmc_impl = MCMCImpl(self.phycas, self)
        mcmc_impl.run()
        
#                ("default_model",          'hky',                           "Can be 'jc', 'hky' or 'gtr'"),
#                ("relrate_prior",          ExponentialDist(1.0),            "The prior distribution for individual GTR relative rate parameters"),
#                ("relrates",               [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] , "The starting values for GTR relative rates"),
#                ("fix_relrates",           False,                           "If True, GTR relative rates will not be modified during the course of an MCMC analysis", BoolArgValidate),
#                ("kappa_prior",            ExponentialDist(1.0),            "The prior distribution for the kappa parameter in an HKY model"),
#                ("kappa",                  4.0,                             "The starting value for the kappa parameter in an HKY model", FloatArgValidate(min=0.01)),
#                ("fix_kappa",              False,                           "If True, the HKY kappa parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
#                ("num_rates",              1,                               "The number of relative rates used for the discrete gamma rate heterogeneity submodel; default is rate homogeneity (i.e. 1 rate)", IntArgValidate(min=1)),
#                ("gamma_shape_prior",      ExponentialDist(1.0),            "The prior distribution for the shape parameter of the gamma among-site rate distribution"),
#                ("gamma_shape",            0.5,                             "The starting value for the gamma shape parameter", FloatArgValidate(min=0.01)),
#                ("fix_shape",              False,                           "If True, the gamma shape parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
#                ("use_inverse_shape",      False,                           "If True, gamma_shape_prior is applied to 1/shape rather than shape", BoolArgValidate),
#                ("pinvar_model",           False,                           "If True, an invariable sites submodel will be applied and the parameter representing the proportion of invariable sites will be estimated", BoolArgValidate),
#                ("pinvar_prior",           BetaDist(1.0, 1.0),              "The prior distribution for pinvar, the proportion of invariable sites parameter"),
#                ("pinvar",                 0.2,                             "The starting value of pinvar, the proportion of invariable sites parameter", ProbArgValidate()),
#                ("fix_pinvar",             False,                           "If True, the proportion of invariable sites parameter (pinvar) will not be modified during the course of an MCMC analysis", BoolArgValidate),
#                ("base_freq_param_prior",  ExponentialDist(1.0),            "The prior distribution for the individual base frequency parameters; these parameters, when normalized to sum to 1, represent the equilibrium proportions of the nucleotide states"),
#                ("base_freqs",             [1.0, 1.0, 1.0, 1.0],            "The starting values for the four base frequency parameters"),
#                ("fix_freqs",              False,                           "If True, the base frequencies will not be modified during the course of an MCMC analysis", BoolArgValidate),
