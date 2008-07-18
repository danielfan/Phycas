from copy import copy
from phycas.Utilities.PhycasCommand import *
from phycas.ProbDist import BetaDist, ExponentialDist, InverseGammaDist
class Model(PhycasCommand):
    def __init__(self, p):
        args = ( 
                ("type",                   'hky',                           "Can be 'jc', 'hky' or 'gtr'", EnumArgValidate(['jc', 'hky', 'gtr'])),
                ("relrate_prior",          ExponentialDist(1.0),            "The prior distribution for individual GTR relative rate parameters"),
                ("relrates",               [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] , "The current values for GTR relative rates. These should be specified in this order: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T."),
                ("fix_relrates",           False,                           "If True, GTR relative rates will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("kappa_prior",            ExponentialDist(1.0),            "The prior distribution for the kappa parameter in an HKY model"),
                ("kappa",                  4.0,                             "The current value for the kappa parameter in an HKY model", FloatArgValidate(min=0.01)),
                ("fix_kappa",              False,                           "If True, the HKY kappa parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("num_rates",              1,                               "The number of relative rates used for the discrete gamma rate heterogeneity submodel; default is rate homogeneity (i.e. 1 rate)", IntArgValidate(min=1)),
                ("gamma_shape_prior",      ExponentialDist(1.0),            "The prior distribution for the shape parameter of the gamma among-site rate distribution"),
                ("gamma_shape",            0.5,                             "The current value for the gamma shape parameter", FloatArgValidate(min=0.01)),
                ("fix_shape",              False,                           "If True, the gamma shape parameter will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("use_inverse_shape",      False,                           "If True, gamma_shape_prior is applied to 1/shape rather than shape", BoolArgValidate),
                ("pinvar_model",           False,                           "If True, an invariable sites submodel will be applied and the parameter representing the proportion of invariable sites will be estimated", BoolArgValidate),
                ("pinvar_prior",           BetaDist(1.0, 1.0),              "The prior distribution for pinvar, the proportion of invariable sites parameter"),
                ("pinvar",                 0.2,                             "The current value of pinvar, the proportion of invariable sites parameter", ProbArgValidate()),
                ("fix_pinvar",             False,                           "If True, the proportion of invariable sites parameter (pinvar) will not be modified during the course of an MCMC analysis", BoolArgValidate),
                ("base_freq_param_prior",  ExponentialDist(1.0),            "The prior distribution for the individual base frequency parameters; these parameters, when normalized to sum to 1, represent the equilibrium proportions of the nucleotide states"),
                ("base_freqs",             [1.0, 1.0, 1.0, 1.0],            "The current values for the four base frequency parameters"),
                ("fix_freqs",              False,                           "If True, the base frequencies will not be modified during the course of an MCMC analysis", BoolArgValidate),
                )
        PhycasCommand.__init__(self, p, args, "model", "Defines a substitution model.")

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
        








