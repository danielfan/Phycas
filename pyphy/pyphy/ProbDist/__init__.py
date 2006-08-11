#import warnings
#warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from pyphy.Conversions import *
from _ProbDist import *
from _Lot import *
from _SliceSampler import *
from _GammaDist import *
from _ExponentialDist import *
from _InverseGammaDist import *
from _UniformDist import *
from _BinomialDist import *
from _BernoulliDist import *
from _BetaDist import *
from _DirichletDist import *

#print 'importing ProbDist...'

def testExamples():
    import doctest
    doctest.testfile('_BernoulliDist.py')
    doctest.testfile('_BetaDist.py')
    doctest.testfile('_BinomialDist.py')
    doctest.testfile('_DirichletDist.py')
    doctest.testfile('_ExponentialDist.py')
    doctest.testfile('_GammaDist.py')
    doctest.testfile('_InverseGammaDist.py')
    doctest.testfile('_Lot.py')
    doctest.testfile('_SliceSampler.py')
    doctest.testfile('_UniformDist.py')
