import sys
import doctest
import phycas
from phycas.DataMatrix import testExamples as dataMatrixTestExamples
import os
a = [0, 0]
r = phycas.ProbDist.testExamples()
a[0] += r[0] ; a[1] += r[1] 
r = phycas.PDFGen.testExamples()
a[0] += r[0] ; a[1] += r[1] 
r = phycas.Phylogeny.testExamples()
a[0] += r[0] ; a[1] += r[1] 
r = phycas.ReadNexus.testExamples()
r = doctest.testfile(os.path.join("..", "DataMatrix", '_DataMatrix.py'))
a[0] += r[0] ; a[1] += r[1] 
r = phycas.Likelihood.testExamples()
a[0] += r[0] ; a[1] += r[1] 
r = doctest.testfile('../Phycas/Phycas.py')
a[0] += r[0] ; a[1] += r[1] 

print "Finished testing examples"
print "%d failures out of %d tests." % (a[0], a[1])
if sys.platform == "win32":
    raw_input("Press any key to quit")
sys.exit(a[0])
