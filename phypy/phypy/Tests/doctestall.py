import doctest
import pyphy

pyphy.ProbDist.testExamples()
pyphy.Phylogeny.testExamples()
pyphy.ReadNexus.testExamples()
pyphy.DataMatrix.testExamples()
pyphy.Likelihood.testExamples()
doctest.testfile('../Phycas/Phycas.py')

print 'Finished testing examples'
