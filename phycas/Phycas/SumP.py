from phycas import randomtree
from phycas.Utilities.PhycasCommand import *
from phycas.Phycas.SumPImpl import ParamSummarizer
import copy

class SumP(PhycasCommand):
    def __init__(self):
        args = (   
                   ("burnin",   1, "Number of lines from the input list of trees to skip (first line stores starting parameter values, so this value should normally be at least 1)", IntArgValidate(min=0)),
                   ("file",    "", "Name of file containing sampled parameter values", FileExistsValidate),
                )
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["log"]
        logf_spec = TextOutputSpec(prefix='sumpoutput', help_str="The file specified by this setting saves the console output generated by sumt(). If set to None, console output will not be saved to a file.")
        o.__dict__["log"] = logf_spec
        PhycasCommand.__init__(self, args, "sump", "The sump command is used to summarize a sample of parameter values produced by an MCMC simulation.", o)
        
    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep SumP's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return True
        
    hidden = staticmethod(hidden)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        param_summarizer = ParamSummarizer(c)
        param_summarizer.run()        
