from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class PS(PhycasCommand):
    def __init__(self):
        args = (   ("nbetavals", 101, "The number of values beta will take on during the run; for example, if this value is 4, then beta will take on these values: 1, 2/3, 1/3, 0", IntArgValidate(min=1)),
                   ("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   ("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                )
        # Specify output options
        PhycasCommand.__init__(self, args, "ps", "Performs path sampling (thermodynamic integration) for purposes of estimating the marginal likelihood of the current model.")

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
        mcmc()
        mcmc.doing_path_sampling = False
    
    # probably not good to give users this choice
    #("toward_posterior", True, "If True, chain will start with beta = 0.0 (exploring prior) and end up exploring posterior; otherwise, chain will begin by exploring posterior and end exploring prior", BoolArgValidate),

    # can just use mcmc's normal burnin, sample_every and ncycles for these:
    #("burnin", 1000, "Number of cycles used to equilibrate before increasing beta", IntArgValidate(min=0)),
    #("sample_every", 1, "Determines number of times likelihood will be sampled. Number of samples for each value of beta will be ps_Q/sample_every", IntArgValidate(min=0)),
    #("nbetacycles", 100, "Number of cycles between changes in beta", IntArgValidate(min=1)),