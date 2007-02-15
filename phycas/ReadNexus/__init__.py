import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _ReadNexus import *
from _NexusReader import *

#print 'importing ReadNexus...'

def testExamples():
    import doctest
    doctest.testfile('_NexusReader.py')
