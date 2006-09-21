import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phypy.Conversions import *
from _LikelihoodBase import *
from _TreeLikelihood import *
from _Model import *
from _MCMCChainManager import *
from _SimData import *
from _TopoPriorCalculator import *
from _QMatrix import *

#print 'importing Likelihood...'

def testExamples():
    import doctest
    doctest.testfile('_MCMCChainManager.py')
    doctest.testfile('_Model.py')
    doctest.testfile('_QMatrix.py')
    doctest.testfile('_SimData.py')
    doctest.testfile('_TopoPriorCalculator.py')
    doctest.testfile('_TreeLikelihood.py')