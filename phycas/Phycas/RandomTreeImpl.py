import copy, os, sys, math, random, copy
from phycas.Utilities.PhycasCommand import _value_for_user
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Utilities.io import TreeCollection

class TreeSimulator(CommonFunctions, TreeCollection):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Saves consensus tree and QQ plots in pdf files.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object by assigning supplied phycas object
        to a data member variable.
        
        """
        CommonFunctions.__init__(self, opts)
        self.r = self._getLot()
        self.opts.edgelen_dist.setLot(self.r)

    def __str__(self):
        return "Collection of trees simulated by the Yule process with edgelen_dist = %s" % _value_for_user(self.opts.edgelen_dist)

    def writeTree(self, tree, name="", rooted=None):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")

    def finish(self):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")
    def __iter__(self):
        return self
    def next(self):
        t = Phylogeny.Tree()
        Phylogeny.TreeManip(t).randomTree(
                self.ntax,           # number of tips
                self.r,                     # pseudorandom number generator
                self.opts.edgelen_dist, # distribution from which to draw starting edge lengths
                False)


