import doctest

#raw_input('About to test ProbDist')
import ProbDist
ProbDist.testExamples()

#raw_input('About to test Phylogeny')
import Phylogeny
Phylogeny.testExamples()

#raw_input('About to test ReadNexus')
import ReadNexus
ReadNexus.testExamples()

#raw_input('About to test DataMatrix')
import DataMatrix
DataMatrix.testExamples()

#raw_input('About to test Likelihood')
import Likelihood
Likelihood.testExamples()

#raw_input('About to test Phycas.py')
doctest.testfile('Phycas/Phycas.py')

raw_input('Press return to quit')
