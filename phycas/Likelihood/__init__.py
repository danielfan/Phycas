import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _LikelihoodBase import *
from _TreeLikelihood import *
from _Model import *
from _MCMCChainManager import *
from _SimData import *
from _TopoPriorCalculator import *
from _QMatrix import *
#from _SamcMove import *

#print 'importing Likelihood...'

def testExamples():
    import doctest
    a = [0,0]
    r = doctest.testfile('_MCMCChainManager.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_Model.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_QMatrix.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_SimData.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_TopoPriorCalculator.py')
    a[0] += r[0] ; a[1] += r[1]
    r = doctest.testfile('_TreeLikelihood.py')
    a[0] += r[0] ; a[1] += r[1]
    return tuple(a)
