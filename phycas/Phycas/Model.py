import copy
from phycas.Utilities.PhycasCommand import *
from phycas.ProbDist import Beta, Exponential, InverseGamma
class Model(PhycasCommand):
    def __init__(self):
        args = ( 
                ("type",                   'hky',                           "Can be 'jc', 'hky' or 'gtr'", EnumArgValidate(['jc', 'hky', 'gtr'])),
                ("relrate_prior",          Exponential(1.0),                "The prior distribution for individual GTR relative rate parameters"),
                ("relrates",               [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] , "The current values for GTR relative rates. These should be specified in this order: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T."),
                ("fix_relrates",           False,                           "If True, GTR relative rates will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("kappa_prior",            Exponential(1.0),                "The prior distribution for the kappa parameter in an HKY model"),
                ("kappa",                  4.0,                             "The current value for the kappa parameter in an HKY model", FloatArgValidate(greaterthan=0.0)),
                ("fix_kappa",              False,                           "If True, the HKY kappa parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("num_rates",              1,                               "The number of relative rates used for the discrete gamma rate heterogeneity submodel; default is rate homogeneity (i.e. 1 rate)", IntArgValidate(min=1)),
                ("gamma_shape_prior",      Exponential(1.0),                "The prior distribution for the shape parameter of the gamma among-site rate distribution"),
                ("gamma_shape",            0.5,                             "The current value for the gamma shape parameter", FloatArgValidate(greaterthan=0.0)),
                ("fix_shape",              False,                           "If True, the gamma shape parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("use_inverse_shape",      False,                           "If True, gamma_shape_prior is applied to 1/shape rather than shape", BoolArgValidate),
                ("pinvar_model",           False,                           "If True, an invariable sites submodel will be applied and the parameter representing the proportion of invariable sites will be estimated", BoolArgValidate),
                ("pinvar_prior",           Beta(1.0, 1.0),                  "The prior distribution for pinvar, the proportion of invariable sites parameter"),
                ("pinvar",                 0.2,                             "The current value of pinvar, the proportion of invariable sites parameter", ProbArgValidate()),
                ("fix_pinvar",             False,                           "If True, the proportion of invariable sites parameter (pinvar) will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("base_freq_param_prior",  Exponential(1.0),                "The prior distribution for the individual base frequency parameters; these parameters, when normalized to sum to 1, represent the equilibrium proportions of the nucleotide states"),
                ("base_freqs",             [1.0, 1.0, 1.0, 1.0],            "The current values for the four base frequency parameters"),
                ("fix_freqs",              False,                           "If True, the base frequencies will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("edgelen_hyperprior",     InverseGamma(2.1,1.0/1.1),       "The prior distribution for the hyperparameter that serves as the mean of an Exponential edge length prior. If set to None, a non-hierarchical model will be used with respect to edge lengths."),
                ("fix_edgelen_hyperparam", False,                           "If True, the hyperparameter that governs the mean of the Exponential edge length prior will be fixed at the value edgelen_hyperparam.", BoolArgValidate),
                ("edgelen_hyperparam",     0.05,                            "The current value of the edge length hyperparameter - setting this currently has no effect", FloatArgValidate(greaterthan=0.0)),
                ("internal_edgelen_dist",  None,                            "Can be used to set a prior distribution for internal edges that differs from that applied to external edges. If this is set to something besides None, you should also set external_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("external_edgelen_dist",  None,                            "Can be used to set a prior distribution for external edges that differs from that applied to internal edges. If this is set to something besides None, you should also set internal_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("edgelen_dist",           Exponential(2.0),                "Sets both internal_edgelen_dist and external_edgelen_dist to the supplied value. Use this setting if you want all edges in the tree to have the same prior distribution. Using this setting will overwrite any values previously supplied for internal_edgelen_dist and external_edgelen_dist"),
                ("fix_edgelens",           False,                           "not yet documented", BoolArgValidate),
                ("starting_edgelen_dist",  Exponential(10.0),               "Used to select the starting edge lengths when starting_tree_source is 'random'"),
                ("use_flex_model",         False,                           "not yet documented", BoolArgValidate),
                ("flex_ncat_move_weight",  1,                               "Number of times each cycle to attempt an ncat move", IntArgValidate(min=0)),
                ("flex_num_spacers",       1,                               "Number of fake rates between each adjacent pair of real rates", IntArgValidate(min=1)),
                ("flex_phi",               0.25,                            "Proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)", ProbArgValidate()),
                ("flex_L",                 1.0,                             "Upper bound of interval used for unnormalized relative rate parameter values", FloatArgValidate(greaterthan=0.0)),
                ("flex_lambda",            1.0,                             "Parameter of Poisson prior on the number of extra categories", FloatArgValidate(greaterthan=0.0)),
                ("flex_prob_param_prior",  Exponential(1.0),                "not yet documented"),
                )
        PhycasCommand.__init__(self, args, "model", "Defines a substitution model.")

    def saveas(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memo):
        c = memo.get(self)
        if c:
            return c
        new_model = Model()
        new_model.type                  = copy.deepcopy(self.type, memo)

        new_model.num_rates             = copy.deepcopy(self.num_rates, memo)
        new_model.use_inverse_shape     = copy.deepcopy(self.use_inverse_shape, memo)
        new_model.pinvar_model          = copy.deepcopy(self.pinvar_model, memo)
        
        new_model.relrate_prior         = copy.deepcopy(self.relrate_prior, memo)
        new_model.kappa_prior           = copy.deepcopy(self.kappa_prior, memo)
        new_model.gamma_shape_prior     = copy.deepcopy(self.gamma_shape_prior, memo)
        new_model.pinvar_prior          = copy.deepcopy(self.pinvar_prior, memo)
        new_model.base_freq_param_prior = copy.deepcopy(self.base_freq_param_prior, memo)

        new_model.relrates              = copy.deepcopy(self.relrates, memo)
        new_model.kappa                 = copy.deepcopy(self.kappa, memo)
        new_model.gamma_shape           = copy.deepcopy(self.gamma_shape, memo)
        new_model.pinvar                = copy.deepcopy(self.pinvar, memo)
        new_model.base_freqs            = copy.deepcopy(self.base_freqs, memo)

        new_model.fix_relrates          = copy.deepcopy(self.fix_relrates, memo)
        new_model.fix_kappa             = copy.deepcopy(self.fix_kappa, memo)
        new_model.fix_shape             = copy.deepcopy(self.fix_shape, memo)
        new_model.fix_pinvar            = copy.deepcopy(self.fix_pinvar, memo)
        new_model.fix_freqs             = copy.deepcopy(self.fix_freqs, memo)
        
        memo[self] = new_model
        return new_model

    def __call__(self, **kwargs):
        self.set(**kwargs)
        return self.saveas()
        








