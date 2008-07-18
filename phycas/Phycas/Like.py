from phycas.Phycas.PhycasCommand import *
from phycas import model
from phycas.ProbDist import ExponentialDist
from phycas.Phycas.LikeImpl import LikeImpl
class Like(PhycasCommand):
    def __init__(self, p):
        args = ( 
                ("data_file_name",  None,                                     "Name of file in which to save simulated data'"),
                ("model",           model,                                    "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("tree_source",     'usertree',                               "If 'usertree', the tree description should be supplied in 'tree_topology'; if 'random', an edge length distribution should be supplied in 'starting_edgelen_dist'", EnumArgValidate(['random','usertree'])),
                ("tree_topology",   '(1:0.02,2:0.02,(3:0.01,4:0.01):0.01)',   "The tree topology (with branch lengths) of the model tree to be used for simulation. Used only if tree_source is 'usertree'"),
                ("starting_edgelen_dist",  ExponentialDist(10.0),             "Used to select the starting edge lengths when tree_source is 'random'"),
                )
        PhycasCommand.__init__(self, p, args, "like", "Calculates the log-likelihood under the current model.")

        # The data members added below should be hidden from the user because they are irrelevant to 
        # computing the likelihood. They must be present, however, because they are referenced in the 
        # LikelihoodCore class, which is also used for Bayesian analyses. 
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["random_seed"] = 0
        self.__dict__["fix_edgelens"]   = False
        self.__dict__["use_flex_model"]  = False
        self.__dict__["uf_num_edges"]   = 50
        self.__dict__["use_unimap"]     = False
        self.__dict__["data_source"]    = 'file'

    def __call__(self, **kwargs):
        self.set(**kwargs)
        calclike = LikeImpl(self.phycas, self)
        return calclike.run()
        
