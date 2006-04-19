#from _DataMatrixBase import *
from _LikelihoodBase import *

class MCMCChainManager(MCMCChainManagerBase):    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|    """    Sorry, no documentation yet.
        """
    def __init__(self):
        """
        Sorry, no documentation yet.
        
        """
        MCMCChainManagerBase.__init__(self)
        
    def finalize(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.finalize(self)
    def recalcEdgeLenPriors(self, mu, var):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.recalcEdgeLenPriors(self, mu, var)
    def calcJointLnPrior(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.calcJointLnPrior(self)
        def addMove(self, move):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addMove(self, move)        def addModelParam(self, param):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addModelParam(self, param)        def addEdgeLenParam(self, param):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addEdgeLenParam(self, param)        def addEdgeLenHyperparam(self, param):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addEdgeLenHyperparam(self, param)        def getLastLnLike(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getLastLnLike(self)
    def getMoves(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getMoves(self)
    def getModelParams(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getModelParams(self)
    def getEdgeLenParams(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getEdgeLenParams(self)
    def getEdgeLenHyperparam(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getEdgeLenHyperparam(self)
    def addMCMCUpdaters(self, model, tree, like, rng, separate_edgelens, max_units, weight):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addMCMCUpdaters(self, model, tree, like, rng, separate_edgelens, max_units, weight)
    def clear(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.clear(self)
    def getAllUpdaters(self):        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getAllUpdaters(self)
