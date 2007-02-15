import sys
import doctest
import phycas

phycas.ProbDist.testExamples()
phycas.Phylogeny.testExamples()
phycas.ReadNexus.testExamples()
phycas.DataMatrix.testExamples()
phycas.Likelihood.testExamples()
doctest.testfile('../Phycas/Phycas.py')

print "Finished testing examples"
if sys.platform == "win32":
    raw_input("Press any key to quit")

