'''Root of the phycas package.'''
#__all__ = [
#    'Conversions',
#    'DataMatrix', 
#    'Likelihood', 
#    'PDFGen',
#    'Phylogeny',
#    'ProbDist', 
#    'ReadNexus', 
#    'PhycasA',
#    'Examples',
#    ]

import Conversions
import DataMatrix
import Likelihood
import PDFGen
import wxPhycas
import Phylogeny
import ProbDist
import ReadNexus
import sys
from Phycas import Phycas
from Phycas.SumT import SumT
_user_ini_checked = False

intercept_python_exceptions = True
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)

phycas = Phycas()
def error_msg(msg):
    sys.stderr.write("Error: %s\n" % msg)

def phycas_except_hook(t, v, tb):
    error_msg(v)

if intercept_python_exceptions:
    sys.excepthook = phycas_except_hook

class Newick(object):
    """A class that holds a newick string to define a tree along with an
    a field `naming` that indicates whether the taxa are numbered from 0, 1 or 
    whether they newick string has taxon labels in it."""
    ZERO_BASED_TAXA_NUMBERS = 0
    ONE_BASED_TAXA_NUMBERS = 1
    TAXA_NAMES = 2
    def __init__(self, newick, naming=1):
        self.newick = newick
        self.naming = naming
    def buildTree(self, tree = None):
        """Calls buildFromString or other appropriate method to construct a tree 
        in the the variable `tree`
        Returns the `tree` instance."""
        if tree is None:
            tree = Phylogeny.Tree()
        tree.buildFromString(self.newick, self.naming == Newick.ZERO_BASED_TAXA_NUMBERS) # tree descriptions from NCL are 1-based not 0-based
        return tree
    def __str__(self):
        return self.newick



sumt = SumT(phycas)

#print 'importing phycas...'
