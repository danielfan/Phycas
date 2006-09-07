import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

#from pyphy.Conversions import *
from _DataMatrixBase import *
from _DataMatrix import *

#print 'importing DataMatrix...'

def testExamples():
    import doctest
    doctest.testfile('_DataMatrix.py')
