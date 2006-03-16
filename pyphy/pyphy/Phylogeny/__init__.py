#import warnings
#warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from Conversions._Conversions import *
from _Phylogeny import *
from _Tree import *
from _TreeManip import *

def testExamples():
    import doctest
    doctest.testfile('_Tree.py')
    doctest.testfile('_TreeManip.py')
    doctest.testfile('doctest_tree.txt')
