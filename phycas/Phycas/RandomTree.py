from phycas import taxa
from phycas.Utilities.PhycasCommand import *
from phycas.Phycas.RandomTreeImpl import TreeSimulator
from phycas.ProbDist import Exponential
import copy
class RandomTree(PhycasCommand):
    def __init__(self):
        args = tuple(
                    PhycasCommand._getRNGOptions() + 
                   [("taxa", taxa, "The names of the taxa to simulate"),
                    ("edgelen_dist",  Exponential(10.0), "Used to generate edge lengths"),                   
                   ]
                   )
        o = PhycasCommandOutputOptions()
        PhycasCommand.__init__(self, args, "randomtree", "Produces a TreeCollection by simulating trees using the Yule Process", o)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        return TreeSimulator(c)

