import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _Phylogeny import *
from _Tree import *
from _TreeManip import *

#print 'importing Phylogeny...'

def testExamples():
    import doctest
    doctest.testfile('_Tree.py')
    doctest.testfile('_TreeManip.py')
