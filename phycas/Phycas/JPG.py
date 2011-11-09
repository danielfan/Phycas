import sys
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Phycas.JPGImpl import JPGImpl
from phycas import P
from phycas.ProbDist import Exponential,BetaPrime

class JPG(PhycasCommand):
    def __init__(self):
        args = tuple(PhycasCommand._getRNGOptions() + [
            ("data_source", P.characters, "The DataSource that provides the data to be analyzed in the JPG analysis. Should be a DataSource object", DataSourceValidate),
            ("tree_source", None, "A TreeCollection that will serve as the source of trees to analyze. If a string is passed in, it is interpreted as a the path to a tree file.", TreeSourceValidate),
            ("scaling_factor_prior", Exponential(1.0), "The prior distribution for the scaling factor parameter that adjusts edge lengths for the rate of the focal character"),
            ("rate_ratio_prior", BetaPrime(1.0,1.0), "The prior distribution for the forward/reverse rate ratio of the focal character (for a model that allows forward and reverse substitution)"),
            ("nreps", 100, "Number of replicates per tree", IntArgValidate(min=1)),
            ("fromtree", 1, "First tree in tree_source to evaluate (1 <= fromtree <= totree <= no. trees in tree_source)", IntArgValidate(min=1)),
            ("totree", sys.maxint, "Last tree in tree_source to evaluate (1 <= fromtree <= totree <= no. trees in tree_source)", IntArgValidate(min=1)),
            ])
            
        # Specify output options
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["sitelike"]
        p = TextOutputSpec(prefix='sitelike', suffix=".txt", help_str="The text file in which all sampled site log-likelihood values are saved.")
        o.__dict__["sitelike"] = p
        PhycasCommand.__init__(self, args, "jpg", "Tests fit of presence/absence character(s) on supplied trees under current model.", o)

        # The data members added below should be hidden from the user because they are for use by phycas developers.
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["uf_num_edges"] = 50      # necessary because LikelihoodCore looks for this variable
        self.__dict__["use_unimap"] = False     # necessary because LikelihoodCore looks for this variable
        
        #self.__dict__["sitelikef"] = None

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep Gogarten's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False
        
    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a CPO object.
        """
        cf = CommonFunctions(self)
        cf.phycassert(self.totree >= self.fromtree, 'totree cannot be less than fromtree')
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        c = copy.deepcopy(self)
        jpg_impl = JPGImpl(c)
        jpg_impl.run()