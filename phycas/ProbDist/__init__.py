import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *
from _Lot import *
from _StopWatch import *
from _SliceSampler import *
from _NormalDist import *
from _GammaDist import *
from _ExponentialDist import *
from _InverseGammaDist import *
from _UniformDist import *
from _ImproperUniformDist import *
from _BinomialDist import *
from _BernoulliDist import *
from _BetaDist import *
from _BetaPrimeDist import *
from _DirichletDist import *
from _MVNormalDist import *
from _RelRateDist import *

#print 'importing ProbDist...'

def testExamples():
    import doctest
    a = [0,0]
    r = doctest.testfile('_BernoulliDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_BetaDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_BetaPrimeDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_BinomialDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_DirichletDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_RelRateDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_MVNormalDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_ExponentialDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_GammaDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_InverseGammaDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_Lot.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_StopWatch.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_SliceSampler.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_UniformDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    r = doctest.testfile('_ImproperUniformDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    return tuple(a)
    


