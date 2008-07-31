from phycas.Utilities.PhycasCommand import *
from phycas import model, randomtree
from phycas.ProbDist import Exponential
from phycas.Phycas.LikeImpl import LikeImpl

class Like(PhycasCommand):
    def __init__(self):
        args = ( 
                 PhycasCommand._getRNGOptions() + 
                [("data_file_name",  None,                                     "Name of file in which to save simulated data'"),
                ("model",           model,                                    "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("tree_source",   randomtree(),                                  "TreeCollection that will provide the tree.", TreeSourceValidate),
                ("starting_edgelen_dist",  Exponential(10.0),                 "Used to select the starting edge lengths when tree_source is 'random'"),
                ]
                )
        PhycasCommand.__init__(self, args, "like", "Calculates the log-likelihood under the current model.")

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
        calclike = LikeImpl(self)
        return calclike.run()
        
