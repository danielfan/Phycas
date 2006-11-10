import doctest
import phypy

phypy.ProbDist.testExamples()
phypy.Phylogeny.testExamples()
phypy.ReadNexus.testExamples()
phypy.DataMatrix.testExamples()
phypy.Likelihood.testExamples()
doctest.testfile('../Phycas/Phycas.py')

raw_input('Finished testing examples\nPress any key to quit')
