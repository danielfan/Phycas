from phycas import P
from phycas.Utilities.PhycasCommand import *
from phycas.Phycas.RandomTreeImpl import TreeSimulator
from phycas.ProbDist import Exponential
import copy
class RandomTree(PhycasCommand):
    def __init__(self):
        args = tuple(
                    PhycasCommand._getRNGOptions() + 
                   [("taxon_labels", P.taxon_labels , "The names of the taxa to simulate"),
                    ("edgelen_dist",  Exponential(10.0), "Used to generate edge lengths"),
                    ("n_trees", 1, "The number of trees to generate (if 0 is specified, a bottomless collection of trees is generated)", IntArgValidate(min=0)),
                    ("n_taxa", 0, "The number of taxa to generate (only used if the taxon_labels attribute is not specified", IntArgValidate(min=0)),
                   ]
                   )
        o = PhycasCommandOutputOptions()
        PhycasCommand.__init__(self, args, "randomtree", "Produces a TreeCollection by simulating trees using the Yule Process", o)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        return TreeSimulator(c)

