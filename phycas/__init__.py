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


def name_of(_obj):
    k=None
    v=None
    for k,v in locals().iteritems():
        if v is _obj and k != "_obj":
            return k
    for k,v in globals().iteritems():
        if v is _obj and k != "_obj":
            return k
    return None

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
            m = msg + '\n'
            if level == OutFilter.WARNINGS:
                m = "\n***** Warning: %s\n" % msg
            elif level > OutFilter.WARNINGS:
                m = "\n***** Error: %s\n" % msg
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

# These globals are set here, so that reading in the startup.py gives
#   experienced users the chance to override the default behavior.
# Mark: it would be better I think to put these in a file that can be
# rewritten by a Settings object. The Settings object could be defined
# like MCMC, SumT, etc., with options such as "intercept_python_exceptions"
# If a user wanted to change, say, help_double_space, he/she could type
# settings.help_double_space=False; settings(). The settings() call would
# overwrite the file containing these settings with the current values
# and immediately reload the module so that the new options would 
# immediately take effect. Can PhycasCommand be modified to allow classes
# that behave like Settings?
release_version = True
help_double_space = True
current_double_space = False
current_follows_help = True
phycassertRaisesException = False
cppCompiledInDebug = False
intercept_python_exceptions = True
default_verbosity_level = OutFilter.NORMAL
_user_ini_checked = False
_use_wx_phycas = False
_check_for_updates = True
_phycas_update_url = "129.237.138.231" # change this to phycas.org url
_phycas_branch = "$HeadURL$"

# Look for a file named svnver.txt in the same directory as this __init__.py file
# The svnver.txt file should contain only the svn version number, and should be
# created automatically when bjam is run (using a make rule in the Jamroot file)
_phycas_revision = '?'
try:
    import os
    svnverfile = 'svnver.txt'
    svnverdir = os.path.dirname(__file__)
    svnverpath = os.path.join(svnverdir,svnverfile)
    _phycas_revision = open(svnverpath, 'r').read()
except:
    print 'Warning: could not find the %s file in this directory:\n  %s' % (svnverfile,svnverdir)
#_phycas_revision = "PHYCAS_SVN_REVISION_NUMBER_HERE"

#from Phycas import Phycas
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.phycas/startup.py")
    if os.path.exists(p):
        execfile(p)
    del p

def useWxPhycas():
    return _use_wx_phycas

from phycas.Utilities.io import getPhycasDir, getPhycasTestData, _runRegressionTests, FileFormats, TreeCollection, DataSource
import Conversions
from DataMatrix._DataMatrix import DataMatrix
import Likelihood
import PDFGen
import Phylogeny
import ProbDist
from ProbDist import Bernoulli, Beta, Binomial, Dirichlet, Exponential, Gamma, ImproperUniform, InverseGamma, Normal, Uniform, Lot
import ReadNexus
import sys, os
from Utilities.PhycasCommand import REPLACE, APPEND, ADD_NUMBER, phycas_help, public
from phycas.Utilities.DefaultData import DefaultData

P = DefaultData.getInstance()

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

try:
    if intercept_python_exceptions:
        sys.excepthook = phycas_except_hook
    import IPython
    _ipython_api = IPython.ipapi.get()
    _ip = _ipython_api.IP
    _ip_phycas_except_hook  = lambda x, a, b, c: phycas_except_hook(a,b,c)
    _ip.set_custom_exc((Exception,), _ip_phycas_except_hook)
    _ip.magic_colors("NoColor")
except:
    pass

from Phycas.RandomTree import RandomTree
randomtree = RandomTree()

from Phycas.Model import Model
model = Model()

from Phycas.SumT import SumT
sumt = SumT()

from Phycas.SumP import SumP
sump = SumP()

from Phycas.MCMC import MCMC
mcmc = MCMC()

from Phycas.Sim import Sim
sim  = Sim()

from Phycas.Like import Like
like = Like()

from Phycas.PS import PS
ps = PS()

from Phycas.CPO import CPO
cpo = CPO()

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

def at(addr):
    """Look up an object by its id."""
    import gc
    for o in gc.get_objects():
        if id(o) == addr:
            return o
    return None
from phycas.Utilities.GlobalState import readFile

RNG = Lot
    
print """
  /////////////////////////////
 ///// Welcome to Phycas /////
/////////////////////////////
Version %s

Phycas is written by Paul O. Lewis, Mark Holder and David Swofford

Phycas is distributed under the GNU Public License (see License file for more
information).
""" % (_phycas_revision,)


