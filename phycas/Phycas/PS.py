from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class PS(PhycasCommand):
    def __init__(self):
        args = (   ("nbetavals", 101, "The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 1, 2/3, 1/3, 0", IntArgValidate(min=1)),
                   ("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("shape1", 1.0, "The first shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                   ("shape2", 1.0, "The second shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from minbeta to maxbeta.", FloatArgValidate(greaterthan=0.0)),
                )
        # Specify output options
        self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
        PhycasCommand.__init__(self, args, "ps", "Performs path sampling (thermodynamic integration) for purposes of estimating the marginal likelihood of the current model.")

        # The data members added below are hidden from the user because they are set when the mcmc command runs
        self.__dict__["sampled_likes"] = None
        self.__dict__["sampled_betas"] = None

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep PS's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return True
    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a PS object.
        """
        cf = CommonFunctions(self)
        cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for path sampling')
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        mcmc.doing_path_sampling = True
        mcmc.ps_nbetavals = self.nbetavals
        mcmc.ps_maxbeta = self.maxbeta
        mcmc.ps_minbeta = self.minbeta
        mcmc.ps_shape1 = self.shape1
        mcmc.ps_shape2 = self.shape2
        mcmc.ps_heating_likelihood = True
        mcmc()
        mcmc.ps_heating_likelihood = False
        mcmc.doing_path_sampling = False
        self.sampled_betas = mcmc.ps_sampled_betas
        self.sampled_likes = mcmc.ps_sampled_likes
