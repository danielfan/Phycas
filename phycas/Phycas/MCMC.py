from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import model, randomtree, P
from phycas.Phycas.MCMCImpl import MCMCImpl
#from phycas.ProbDist import Beta, Exponential, InverseGamma
import copy

class MCMC(PhycasCommand):
    def __init__(self):
        args = tuple(
                PhycasCommand._getRNGOptions() + 
                [("burnin",                 0,                               "The number of update cycles to ignore before sampling begins. If burnin=1000 and ncycles=10000, the total number of cycles will be 11000.", IntArgValidate(min=0)),
                ("ncycles",                10000,                           "The number of update cycles (a cycle is analogous to, but different than, a 'generation' in MrBayes; Phycas does in one cycle what MrBayes does in about 100 generations for a simple model such as JC)", IntArgValidate(min=0)),
                ("sample_every",           100,                             "The current tree topology and model parameter values will be sampled after this many cycles have elapsed since the last sample was taken", IntArgValidate(min=0)),
                ("report_every",           100,                             "A progress report will be displayed after this many cycles have elapsed since the last progress report", IntArgValidate(min=0)),
                ("verbose",                True,                            "You will get more output if True, less output if False", BoolArgValidate),
                ("quiet",                  False,                           "If True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well", BoolArgValidate),
                ("model",                  model,                           "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("data_source",            P.characters,                    "The DataSource that provides the data, if any, to be used in the MCMC analysis. Should be a DataSource object", DataSourceValidate),
                ("starting_tree_source",   None,                            "A TreeCollection that will serve as the source of the starting tree topology. If a string is passed in, it is interpretted as a the path to a file with trees.", TreeSourceValidate),
                ("tree_topology",          None,                            "Unused unless starting_tree_source is 'usertree', in which case this should be a standard newick string representation of the tree topology; e.g. '(A:0.01,B:0.071,(C:0.013,D:0.021):0.037)'"),
                ("fix_topology",           False,                           "If True, an EdgeMove move will be substituted for the LargetSimonMove, so edge lengths will be updated by slice sampling but the topology will remain unchanged during an MCMC analysis", BoolArgValidate),
                ("ls_move_lambda",         0.2,                             "The value of the tuning parameter for the Larget-Simon move", FloatArgValidate(min=0.01)),
                ("ls_move_weight",         100,                             "Larget-Simon moves will be performed this many times per cycle", IntArgValidate(min=0)),
                ("ls_move_debug",          False,                           "If set to true, TreeViewer will popup on each Larget-Simon move update showing edges affected by the proposed move", BoolArgValidate),
                ("edge_move_lambda",       0.2,                             "The value of the tuning parameter for the EdgeMove", FloatArgValidate(min=0.01)),
                ("edge_move_weight",       0,                               "Only used if fix_topology is True. Makes sense to set this to some multiple of the number of edges since each EdgeMove affects a single randomly-chosen edge ", IntArgValidate(min=0)),
                #("mapping_move_weight",     1,                              "Univent mapping will be performed this many times per cycle", IntArgValidate(min=0)),
                #("unimap_nni_move_weight", 100,                             "Unimap NNI moves will be performed this many times per cycle", IntArgValidate(min=0)),
                ("rel_rate_weight",        1,                               "Updates of GTR relative rates will occur this many times per cycle if relative rates are being updated jointly", IntArgValidate(min=0)),
                ("state_freq_weight",      1,                               "Updates of state frequencies will occur this many times per cycle if state frequencies are being updated jointly", IntArgValidate(min=0)),
                ("tree_scaler_weight",     0,                               "Whole-tree scaling will be performed this many times per cycle", IntArgValidate(min=0)),
                #("use_unimap",              False,                          "if True, MCMC analyses will use the uniformized mapping approach", BoolArgValidate),
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
                ("heating_lambda",         0.2,                             "not yet documented", FloatArgValidate(min=0.01)),
                ("uf_num_edges",           50,                              "Number of edges to traverse before taking action to prevent underflow", IntArgValidate(min=1)),
                ("ntax",                   0,                               "To explore the prior, set to some positive value. Also set data_source to None", IntArgValidate(min=0)),
                ("ndecimals",              8,                               "Number of decimal places used for sampled parameter values", IntArgValidate(min=1)),
                ])

        # Specify output options
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["log", "trees", "params"]
        logf_spec = TextOutputSpec(prefix='mcmcoutput', help_str="The file specified by this setting saves the console output generated by mcmc(). If set to None, console output will not be saved to a file.")
        o.__dict__["log"] = logf_spec
        t = TextOutputSpec(prefix='trees', suffix=".t", help_str="The nexus tree file in which all sampled tree topologies are saved. This file is equivalent to the MrBayes *.t file.")
        o.__dict__["trees"] = t
        p = TextOutputSpec(prefix='params', suffix=".p", help_str="The text file in which all sampled parameter values are saved. This file is equivalent to the MrBayes *.p file.")
        o.__dict__["params"] = p
        PhycasCommand.__init__(self, args, "mcmc", "The mcmc command is used to conduct a Bayesian Markov chain Monte Carlo analysis.", o)

        # The data members added below should be hidden from the user because they are for use by phycas developers.
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["debugging"] = False
        self.__dict__["use_unimap"] = False                 # if True, MCMC analyses will use the uniformized mapping approach.
        self.__dict__["mapping_move_weight"] = 1            # Univent mapping will be performed this many times per cycle
        self.__dict__["unimap_nni_move_weight"] = 100       # Unimap NNI moves will be performed this many times per cycle
        self.__dict__["draw_directly_from_prior"] = True    # If True, MCMCImpl.explorePrior function is used to draw samples directly from the prior during path sampling, which dramatically improves mixing compared to using MCMC proposals to explore the prior
        self.__dict__["nchains"] = 1                        # The number of Markov chains to run simultaneously. One chain serves as the cold chain from which samples are drawn, the other chains are heated to varying degrees and serve to enhance mixing in the cold chain.

        # The data members added below are hidden from the user because they are set by the ps command
        self.__dict__["doing_path_sampling"] = False
        self.__dict__["ps_nbetavals"] = 101
        self.__dict__["ps_maxbeta"] = 1.0
        self.__dict__["ps_minbeta"] = 0.0
        self.__dict__["ps_shape1"] = 1.0
        self.__dict__["ps_shape2"] = 1.0
        self.__dict__["ps_heating_likelihood"] = False
        
        # The data members added below are hidden from the user because they are set by the cpo command
        self.__dict__["saving_sitelikes"] = False
        self.__dict__["sitelikef"] = None
        
        # The data members added below are hidden from the user because they are set when the mcmc command runs
        self.__dict__["ps_sampled_likes"] = None
        self.__dict__["ps_sampled_betas"] = None
        
    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of self.run().
        """
        cf = CommonFunctions(self)
        cf.phycassert(self.ps_nbetavals > 0, 'ps_nbetavals cannot be less than 1')
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        #self.checkSanity()
        c = copy.deepcopy(self)
        mcmc_impl = MCMCImpl(c)
        mcmc_impl.setSiteLikeFile(self.sitelikef)
        mcmc_impl.run()
        self.ps_sampled_betas = mcmc_impl.ps_sampled_betas
        self.ps_sampled_likes = mcmc_impl.ps_sampled_likes
        mcmc_impl.unsetSiteLikeFile()