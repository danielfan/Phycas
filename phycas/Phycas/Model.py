from copy import copy
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
                ("kappa",                  4.0,                             "The current value for the kappa parameter in an HKY model", FloatArgValidate(min=0.01)),
                ("fix_kappa",              False,                           "If True, the HKY kappa parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("num_rates",              1,                               "The number of relative rates used for the discrete gamma rate heterogeneity submodel; default is rate homogeneity (i.e. 1 rate)", IntArgValidate(min=1)),
                ("gamma_shape_prior",      Exponential(1.0),                "The prior distribution for the shape parameter of the gamma among-site rate distribution"),
                ("gamma_shape",            0.5,                             "The current value for the gamma shape parameter", FloatArgValidate(min=0.01)),
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
                ("edgelen_hyperparam",     0.05,                            "The current value of the edge length hyperparameter - setting this currently has no effect", FloatArgValidate(min=0.01)),
                ("internal_edgelen_dist",  None,                            "Can be used to set a prior distribution for internal edges that differs from that applied to external edges. If this is set to something besides None, you should also set external_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("external_edgelen_dist",  None,                            "Can be used to set a prior distribution for external edges that differs from that applied to internal edges. If this is set to something besides None, you should also set internal_edgelen_dist appropriately. Setting the edgelen_dist property sets both external_edgelen_dist and internal_edgelen_dist to the same value"),
                ("edgelen_dist",           Exponential(2.0),                "Sets both internal_edgelen_dist and external_edgelen_dist to the supplied value. Use this setting if you want all edges in the tree to have the same prior distribution. Using this setting will overwrite any values previously supplied for internal_edgelen_dist and external_edgelen_dist"),
                ("fix_edgelens",           False,                           "not yet documented", BoolArgValidate),
                ("starting_edgelen_dist",  Exponential(10.0),               "Used to select the starting edge lengths when starting_tree_source is 'random'"),
                ("use_flex_model",         False,                           "not yet documented", BoolArgValidate),
                ("flex_ncat_move_weight",  1,                               "Number of times each cycle to attempt an ncat move", IntArgValidate(min=0)),
                ("flex_num_spacers",       1,                               "Number of fake rates between each adjacent pair of real rates", IntArgValidate(min=1)),
                ("flex_phi",               0.25,                            "Proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)", ProbArgValidate()),
                ("flex_L",                 1.0,                             "Upper bound of interval used for unnormalized relative rate parameter values", FloatArgValidate(min=0.01)),
                ("flex_lambda",            1.0,                             "Parameter of Poisson prior on the number of extra categories", FloatArgValidate(min=0.01)),
                ("flex_prob_param_prior",  Exponential(1.0),                "not yet documented"),
                )
        PhycasCommand.__init__(self, args, "model", "Defines a substitution model.")

    def saveas(self):
        new_model = Model(self.phycas)
        new_model.type                  = copy(self.type)

        new_model.num_rates             = copy(self.num_rates)
        new_model.use_inverse_shape     = copy(self.use_inverse_shape)
        new_model.pinvar_model          = copy(self.pinvar_model)

        new_model.relrate_prior         = self.relrate_prior.clone()
        new_model.kappa_prior           = self.kappa_prior.clone()
        new_model.gamma_shape_prior     = self.gamma_shape_prior.clone()
        new_model.pinvar_prior          = self.pinvar_prior.clone()
        new_model.base_freq_param_prior = self.base_freq_param_prior.clone()

        new_model.relrates              = copy(self.relrates)
        new_model.kappa                 = copy(self.kappa)
        new_model.gamma_shape           = copy(self.gamma_shape)
        new_model.pinvar                = copy(self.pinvar)
        new_model.base_freqs            = copy(self.base_freqs)

        new_model.fix_relrates          = copy(self.fix_relrates)
        new_model.fix_kappa             = copy(self.fix_kappa)
        new_model.fix_shape             = copy(self.fix_shape)
        new_model.fix_pinvar            = copy(self.fix_pinvar)
        new_model.fix_freqs             = copy(self.fix_freqs)
        
        return new_model
    
    def __call__(self, **kwargs):
        self.set(**kwargs)
        return self.saveas()
        








