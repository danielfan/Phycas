import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _PhylogenyExt import *
from _Tree import *
from _TreeManip import *
from _Split import *

#print 'importing Phylogeny...'

def testExamples():
    import doctest
    a = [0,0]
    r = doctest.testfile('_Tree.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_TreeManip.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_Split.py')
    a[0] += r[0] ; a[1] += r[1]
    return tuple(a)
    
