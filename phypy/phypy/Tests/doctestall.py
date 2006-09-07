import doctest
import phypy

phypy.ProbDist.testExamples()
phypy.Phylogeny.testExamples()
phypy.ReadNexus.testExamples()
phypy.DataMatrix.testExamples()
phypy.Likelihood.testExamples()
doctest.testfile('../Phycas/Phycas.py')

print 'Finished testing examples'
