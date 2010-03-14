import sys
import doctest
import phycas
from phycas.DataMatrix import testExamples as dataMatrixTestExamples
import os

# Set verbose to True if you want to track down which file was being tested
# when an assert was tripped
verbose = False

a = [0, 0]

if verbose: print '...testing examples in ProbDist module...'
r = phycas.ProbDist.testExamples(verbose)
a[0] += r[0] ; a[1] += r[1] 

if verbose: print '...testing examples in PDFGen module...'
r = phycas.PDFGen.testExamples(verbose)
a[0] += r[0] ; a[1] += r[1] 

if verbose: print '...testing examples in Phylogeny module...'
r = phycas.Phylogeny.testExamples(verbose)
a[0] += r[0] ; a[1] += r[1] 

if verbose: print '...testing examples in ReadNexus module...'
r = phycas.ReadNexus.testExamples(verbose)
r = doctest.testfile(os.path.join("..", "DataMatrix", '_DataMatrix.py'))
a[0] += r[0] ; a[1] += r[1] 

if verbose: print '...testing examples in Likelihood module...'
r = phycas.Likelihood.testExamples(verbose)
a[0] += r[0] ; a[1] += r[1] 

if verbose: print '...testing examples in Phycas.py...'
r = doctest.testfile('../Phycas/Phycas.py')
a[0] += r[0] ; a[1] += r[1] 

print "Finished testing examples"
print "%d failures out of %d tests." % (a[0], a[1])
if sys.platform == "win32":
    raw_input("Press any key to quit")
sys.exit(a[0])
