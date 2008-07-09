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

class OutFilter:
    DEBUGGING, VERBOSE, NORMAL, WARNINGS, ERRORS, SILENT = range(6)
    _names = ["DEBUGGING", "VERBOSE", "NORMAL", "WARNINGS", "ERRORS", "SILENT"]
    def to_str(v):
        if v < 0 or v > len(OutFilter._names):
            raise VauleError("Invalid OutFilter code (%s) specified" % str(v))       
        return OutFilter._names[v]
    to_str = staticmethod(to_str)

class OutputFilter(object):
    def __init__(self, level, stream):
        self.level = level
        self.stream = stream
    def _filter_output(self, msg, level):
        if self.level <= level:
            if level < OutFilter.WARNINGS:
                self.stream(msg)
            elif level == OutFilter.WARNINGS:
                self.stream("\n***** Warning: " + msg)
            else:
                self.stream("\n***** Error: " + msg)
    def info(self, msg):
        self._filter_output(msg, OutFilter.NORMAL)
    def warning(self, msg):
        self._filter_output(msg, OutFilter.WARNINGS)
    def error(self, msg):
        self._filter_output(msg, OutFilter.ERRORS)
    def verbose_info(self, msg):
        self._filter_output(msg, OutFilter.VERBOSE)
    def debugging(self, msg):
        self._filter_output(msg, OutFilter.DEBUGGING)

def getDefaultOutFilter():
    global default_verbosity_level
    return default_verbosity_level


# These globals are set here, so that reading in the startup.py gives
#   experienced users the chance to override the default behavior.
intercept_python_exceptions = True
default_verbosity_level = OutFilter.NORMAL
_user_ini_checked = False
_use_wx_phycas = True

from Phycas import Phycas
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)

def useWxPhycas():
    return _use_wx_phycas

import Conversions
import DataMatrix
import Likelihood
import PDFGen
import Phylogeny
import ProbDist
import ReadNexus
import sys
from Phycas.SumT import SumT
from Phycas.MCMC import MCMC
from Phycas.PhycasCommand import FileFormats, REPLACE, APPEND, ADD_NUMBER, phycas_help
# keep the wx import after the reading of the startup so that it can be optional
if useWxPhycas():
    import wxPhycas


phycas = Phycas()
python_help = help
help = phycas_help
def error_msg(msg):
    sys.stderr.write("Error: %s\n" % msg)

def phycas_except_hook(t, v, tb):
    #print '***** error message passed to phycas_except_hook is %s *****' % v
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
mcmc = MCMC(phycas)

#print 'importing phycas...'
