from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc

class Unimap(PhycasCommand):
    def __init__(self):
        args = (
                ("mapping_move_weight",       1,    "Univent mapping will be performed this many times per cycle", IntArgValidate(min=0)),
                ("unimap_nni_move_weight",  100,    "Unimap NNI moves will be performed this many times per cycle", IntArgValidate(min=0)),    
                )
        # Specify output options
        #self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
        PhycasCommand.__init__(self, args, "unimap", "Performs MCMC with uniformization enabled.")


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
        mcmc.use_unimap = True
        mcmc.mapping_move_weight = self.mapping_move_weight
        mcmc.unimap_nni_move_weight = self.unimap_nni_move_weight
        mcmc()
