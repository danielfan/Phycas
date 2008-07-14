import sys
import doctest
import phycas
from phycas.DataMatrix import testExamples as dataMatrixTestExamples
import os
phycas.ProbDist.testExamples()
phycas.PDFGen.testExamples()
phycas.Phylogeny.testExamples()
phycas.ReadNexus.testExamples()
doctest.testfile(os.path.join("..", "DataMatrix", '_DataMatrix.py'))
phycas.Likelihood.testExamples()
doctest.testfile('../Phycas/Phycas.py')

print "Finished testing examples"
if sys.platform == "win32":
    raw_input("Press any key to quit")

