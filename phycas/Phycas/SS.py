from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class SS(PhycasCommand):
    def __init__(self):
        args = (   ("nbetavals", 101, "The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 1, 2/3, 1/3, 0", IntArgValidate(min=1)),
                   #("simpsons_method", True, "If True (recommended), Simpson's method will be used to approximate the integral; otherwise, the trapezoid method will be used.", BoolArgValidate),
                   ("aux", False, "NOT YET IMPLEMENTED, DO NOT USE! If True, steppingstone sampling will be performed using a multivariate normal auxiliary distribution fit to the posterior; if False, steppingstone sampling will be performed without using an auxiliary distribution.", BoolArgValidate),
                   ("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("shape1", 0.3, "The first shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                   ("shape2", 1.0, "The second shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                )
        # Specify output options
        #self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
        PhycasCommand.__init__(self, args, "ss", "Performs steppingstone sampling (and also path sampling/thermodynamic integration) for purposes of estimating the marginal likelihood of the current model.")

        # The data members added below are hidden from the user because they are set when the mcmc command runs
        self.__dict__["sampled_likes"] = None
        self.__dict__["sampled_betas"] = None

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep SS's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False
        
    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a SS object.
        """
        cf = CommonFunctions(self)
        cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for steppingstone sampling')
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        mcmc.doing_steppingstone_sampling = True
        #mcmc.ss_simpsons_method = self.simpsons_method
        mcmc.ss_nbetavals = self.nbetavals
        mcmc.ss_maxbeta = self.maxbeta
        mcmc.ss_minbeta = self.minbeta
        mcmc.ss_shape1 = self.shape1
        mcmc.ss_shape2 = self.shape2
        if self.aux:
            mcmc.ss_heating_likelihood = False
            mcmc.ss_aux = True
        else:
            mcmc.ss_heating_likelihood = True
            mcmc.ss_aux = False
        mcmc()
        mcmc.ss_heating_likelihood = False
        mcmc.doing_steppingstone_sampling = False
        mcmc.ss_aux = False
        self.sampled_betas = mcmc.ss_sampled_betas
        self.sampled_likes = mcmc.ss_sampled_likes
