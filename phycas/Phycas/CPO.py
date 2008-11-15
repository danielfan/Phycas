from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Phycas.CPOImpl import CPOImpl
from phycas import mcmc

class CPO(PhycasCommand):
    def __init__(self):
        # **** IMPORTANT **** The arguments below are ignored for now - left over from ancestral PS command - replace with something meaningful
        args = (   ("nbetavals", 101, "The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 1, 2/3, 1/3, 0", IntArgValidate(min=1)),
                   ("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("shape1", 0.3, "The first shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                   ("shape2", 1.0, "The second shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                )
        # Specify output options
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["sitelike"]
        p = TextOutputSpec(prefix='sitelike', suffix=".txt", help_str="The text file in which all sampled site log-likelihood values are saved.")
        o.__dict__["sitelike"] = p
        PhycasCommand.__init__(self, args, "cpo", "Performs a Conditional Predictive Ordinate (CPO) analysis to determine the relative fit of the model to individual sites/characters.", o)

        # The data members added below are hidden from the user because they are set when the mcmc command runs
        #self.__dict__["sampled_likes"] = None
        #self.__dict__["sampled_betas"] = None
        
        self.__dict__["sitelikef"] = None

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep CPO's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return True
    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a CPO object.
        """
        #cf = CommonFunctions(self)
        #cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for path sampling')
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        c = copy.deepcopy(self)
        cpo_impl = CPOImpl(c)
        cpo_impl.siteLikeFileOpen()
        if cpo_impl.sitelikef is not None:
            mcmc.sitelikef = cpo_impl.sitelikef
            mcmc.saving_sitelikes = True
            mcmc()
            mcmc.saving_sitelikes = False
            cpo_impl.siteLikeFileClose()
        else:
            print 'Could not run the cpo command because the sitelike file could not be opened'
        