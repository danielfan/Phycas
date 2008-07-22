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
        self.other = []
    def add_mirror(self, o):
        if not o in self.other:
            self.other.append(o)
    def remove_mirror(self, o):
        try:
            n = self.other.index(o)
            self.other.pop(n)
        except:
            return

    def _filter_output(self, msg, level):
        if self.level <= level:
            m = msg
            if level == OutFilter.WARNINGS:
                m = "\n***** Warning: " + msg
            elif level > OutFilter.WARNINGS:
                m = "\n***** Error: " + msg
            self.stream(m)
            for o in self.other:
                o.write(m)
    def info(self, msg):
        self._filter_output(msg, OutFilter.NORMAL)
    def warning(self, msg):
        self._filter_output(msg, OutFilter.WARNINGS)
    def error(self, msg):
        self._filter_output(msg, OutFilter.ERRORS)
    def abort(self, msg):
        self._filter_output(msg, OutFilter.ERRORS)
        sys.exit('\n***** Fatal error: %s' % msg)
    def verbose_info(self, msg):
        self._filter_output(msg, OutFilter.VERBOSE)
    def debugging(self, msg):
        self._filter_output(msg, OutFilter.DEBUGGING)
    def phycassert(self, assumption, msg):
        if not assumption:
            if phycassertRaisesException:
                raise AssertionError(msg)
            sys.exit('Error: ' + msg)

def getDefaultOutFilter():
    global default_verbosity_level
    return default_verbosity_level


# These globals are set here, so that reading in the startup.py gives
#   experienced users the chance to override the default behavior.
help_double_space = True
current_double_space = True
current_follows_help = True
phycassertRaisesException = False
cppCompiledInDebug = False
intercept_python_exceptions = True
default_verbosity_level = OutFilter.NORMAL
_user_ini_checked = False
_use_wx_phycas = True
_check_for_updates = True
_phycas_update_url = "129.237.138.231" # change this to phycas.org url
_phycas_branch = "$HeadURL$"
_phycas_revision = "$Revision: 733 $"

from Phycas import Phycas
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)

def useWxPhycas():
    return _use_wx_phycas

from phycas.Utilities.io import getPhycasDir, getPhycasTestData, _runRegressionTests
import Conversions
from DataMatrix._DataMatrix import DataMatrix
import Likelihood
import PDFGen
import Phylogeny
import ProbDist
from ProbDist import BernoulliDist, BetaDist, BinomialDist, DirichletDist, ExponentialDist, GammaDist, ImproperUniformDist, InverseGammaDist, NormalDist, UniformDist
import ReadNexus
import sys, os
from Utilities.PhycasCommand import FileFormats, REPLACE, APPEND, ADD_NUMBER, phycas_help
# keep the wx import after the reading of the startup so that it can be optional
if useWxPhycas():
    import wxPhycas

#phycas = Phycas()
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

from Phycas.Model import Model
model = Model()

from Phycas.SumT import SumT
sumt = SumT()

from Phycas.MCMC import MCMC
mcmc = MCMC()

from Phycas.Sim import Sim
sim  = Sim()

from Phycas.Like import Like
like = Like()

from Phycas.PS import PS
ps = PS()

def simpleOutputter(msg):
    print msg

def runTests():
    #output_stream = OutputFilter(getDefaultOutFilter(), phycas.output)
    output_stream = OutputFilter(getDefaultOutFilter(), simpleOutputter)
    _runRegressionTests(output_stream)

if _check_for_updates and sys.argv and not sys.argv[0]:
    import phycas.Utilities.PhycasUpdateCheck as PhycasUpdateCheck
    #output_stream = OutputFilter(getDefaultOutFilter(), phycas.output)
    output_stream = OutputFilter(getDefaultOutFilter(), simpleOutputter)
    PhycasUpdateCheck.runPhycasUpdateChecker(output_stream, _phycas_update_url, _phycas_branch, _phycas_revision)

def touch(f):
    "updates the time stamp of the path `f` or creates a file with that uri."
    if os.path.exists(f):
        os.utime(f, None)
    else:
        fo = file(f, "w")
        fo.close()
