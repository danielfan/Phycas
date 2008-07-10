from phycas.Phycas.PhycasCommand import *
from phycas.ProbDist import ExponentialDist
from phycas.Phycas.SimImpl import SimImpl
class Sim(PhycasCommand):
    def __init__(self, p):
        args = (   
                ("file_name",       'simulated.nex',                          "Name of file in which to save simulated data'"),
                ("taxon_labels",    ['taxon1', 'taxon2', 'taxon3', 'taxon4'], "Names to use for taxa in simulated data set (number of labels defined determines the number of taxa in the simulated dataset)"),
                ("nchar",           1000,                                     "Number of characters to generate", IntArgValidate(min=1)),
                ("tree_source",     'usertree',                               "Currently, only 'usertree' is allowed, with the tree description supplied via tree_topology"),
                ("tree_topology",   '(1:0.02,2:0.02,(3:0.01,4:0.01):0.01)',   "The tree topology (with branch lengths) of the model tree to be used for simulation. Used only if tree_source is 'usertree'"),
                ("random_seed",     0,                                        "Determines the random number seed used; specify 0 to generate seed automatically from system clock", IntArgValidate(min=0)),
                ("default_model",   'hky',                                    "Can be 'jc', 'hky' or 'gtr'"),
                ("relrates",        [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] ,          "The GTR relative rates to use for simulation. Only used if default_model is 'gtr'"),
                ("kappa",           4.0,                                      "The kappa (transition/transversion rate ratio) parameter. Only used if default_model is 'hky'", FloatArgValidate(min=0.01)),
                ("num_rates",       1,                                        "The number of relative rates used for the discrete gamma rate heterogeneity submodel; equal rates are used if num_rates = 1", IntArgValidate(min=1)),
                ("gamma_shape",     0.5,                                      "The gamma shape parameter. Only used if num_rates > 1", FloatArgValidate(min=0.01)),
                ("pinvar_model",    False,                                    "If True, an invariable sites submodel will be applied", BoolArgValidate),
                ("pinvar",          0.2,                                      "The proportion of invariable sites. Only used if use_pinvar is True", ProbArgValidate()),
                ("base_freqs",      [1.0, 1.0, 1.0, 1.0],                     "The four base frequency parameters. These will be normalized before they are used, and thus do not need to sum to 1.0"),
                ("starting_edgelen_dist",  ExponentialDist(10.0),           "Used to select the starting edge lengths when tree_source is 'random'"),
                )                                                  
        PhycasCommand.__init__(self, p, args, "simulate", "Simulates DNA sequences according to the specified model.")
        
        # The data members added below should be hidden from the user because they are irrelevant to simulating data
        # They must be present, however, because they are referenced in the LikelihoodCore class, which is also used
        # for likelihood and Bayesian analyses. The roundabout way of introducing these data members is necessary because
        # PhycasCommand.__setattr__ tried to prevent users from adding new data members (to prevent accidental
        # misspellings from causing problems)
        self.__dict__["use_flex_model"] = False
        self.__dict__["fix_kappa"]      = False
        self.__dict__["fix_relrates"]   = False
        self.__dict__["fix_freqs"]      = False
        self.__dict__["fix_edgelens"]   = False
        self.__dict__["uf_num_edges"]   = False
        self.__dict__["use_unimap"]     = False
        self.__dict__["data_source"]    = None

    def __call__(self, **kwargs):
        self.set(**kwargs)
        simulate = SimImpl(self.phycas, self)
        simulate.run()
        
