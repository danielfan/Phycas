from phycas.Phycas.PhycasCommand import *
from phycas import model
from phycas.ProbDist import ExponentialDist
from phycas.Phycas.SimImpl import SimImpl
class Sim(PhycasCommand):
    def __init__(self, p):
        args = (   
                ("file_name",       'simulated.nex',                          "Name of file in which to save simulated data'"),
                ("model",           model,                                    "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("taxon_labels",    ['taxon1', 'taxon2', 'taxon3', 'taxon4'], "Names to use for taxa in simulated data set (number of labels defined determines the number of taxa in the simulated dataset)"),
                ("nchar",           1000,                                     "Number of characters to generate", IntArgValidate(min=1)),
                ("tree_source",     'usertree',                               "If 'usertree', the tree description should be supplied in 'tree_topology'; if 'random', an edge length distribution should be supplied in 'starting_edgelen_dist'", EnumArgValidate(['random','usertree'])),
                ("tree_topology",   '(1:0.02,2:0.02,(3:0.01,4:0.01):0.01)',   "The tree topology (with branch lengths) of the model tree to be used for simulation. Used only if tree_source is 'usertree'"),
                ("random_seed",     0,                                        "Determines the random number seed used; specify 0 to generate seed automatically from system clock", IntArgValidate(min=0)),
                ("starting_edgelen_dist",  ExponentialDist(10.0),             "Used to select the starting edge lengths when tree_source is 'random'"),
                )                                                  
        PhycasCommand.__init__(self, p, args, "sim", "Simulates DNA sequences according to the specified model.")
        
        # The data members added below should be hidden from the user because they are irrelevant to simulating data
        # They must be present, however, because they are referenced in the LikelihoodCore class, which is also used
        # for likelihood and Bayesian analyses. 
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["use_flex_model"] = False
        self.__dict__["fix_edgelens"]   = False
        self.__dict__["uf_num_edges"]   = 50
        self.__dict__["use_unimap"]     = False
        self.__dict__["data_source"]    = None

    def __call__(self, **kwargs):
        self.set(**kwargs)
        simulate = SimImpl(self.phycas, self)
        simulate.run()
        
