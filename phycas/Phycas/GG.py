from phycas.Utilities.PhycasCommand import *
from phycas import model, randomtree
from phycas.ProbDist import Exponential
from phycas.Phycas.GGImpl import GelfandGhosh

class GG(PhycasCommand):
    def __init__(self):
        args = tuple([   
                ("paramfile",       None,           "Name of file containing parameter samples'"),              #was gg_pfile
                ("treefile",        None,           "Name of file containing tree samples'"),                   #was gg_tfile
                ("nreps",           1,              "The number of replicate simulations to do every MCMC sample", IntArgValidate(min=0)),  #was gg_burnin
                ("burnin",          1,              "Number of starting samples to skip when computing Gelfand-Ghosh measures", IntArgValidate(min=1)),
                ("kvalues",         [1.0],          "Vector of overall measures (one for each k in kvalues)"),  # was gg_kvect
                ("postpred_prefix", 'pp',           "Prefix to use for posterior predictive dataset filenames (no datasets will be saved if set to None)"),  #was gg_postpred_prefix,gg_save_postpreds
                ] + PhycasCommand._getRNGOptions()
                )                                                  
        PhycasCommand.__init__(self, args, "sim", "Simulates DNA sequences according to the specified model.")
        
        # The data members added below should be hidden from the user because they are irrelevant to simulating data
        # They must be present, however, because they are referenced in the LikelihoodCore class, which is also used
        # for likelihood and Bayesian analyses. 
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["gg_Pm"]                = 0.0       # Penalty component (same for all k)
        self.__dict__["gg_Gm"]                = []        # Vector of goodness-of-fit components (one for each k in gg_kvect)
        self.__dict__["gg_Dm"]                = []        # Vector of overall measures (one for each k in gg_kvect)
        self.__dict__["gg_outfile"]           = 'gg.txt'  # File in which to save gg results (use None to not save results)
        self.__dict__["gg_bin_patterns"]      = False     # If True, patterns will be classified into 7 bins, corresponding to 'A only', 'C only', 'G only', 'T only', 'any 2 states', 'any 3 states' and 'any 4 states'. Gelfand-Ghosh statistics will be computed on this vector of counts instead of the complete vector of pattern counts. Can only be used for DNA/RNA data.
        self.__dict__["gg_bincount_filename"] = None      # If not None, and if gg_bin_patterns is True, the binned counts for the original dataset and all posterior predictive data sets will be saved to a file by this name
        
    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep Sim's name from being displayed 
        in the list of classes displayed when users type help. Change the return value to 
        False when it is ready to be advertised.
        """
        return True
        
    hidden = staticmethod(hidden)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        gg = GelfandGhosh(c)
        gg.run()
        return gg
