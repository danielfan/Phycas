import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.Conversions import *
from _PDFGen import *
from _PDFGenerator import *

def testExamples():
    import doctest
    doctest.testfile('_PDFGenerator.py')
