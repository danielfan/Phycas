#import Conversions

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

# Test examples within the example apps
import doctest
#doctest.testfile('../apps/ExplorePrior/ExplorePrior.py')
#doctest.testfile('../apps/FixedParams/FixedParams.py')
#doctest.testfile('../apps/GelfandGhosh/GelfandGhosh.py')
#doctest.testfile('../apps/LikelihoodTest/LikelihoodTest.py')
doctest.testfile('../apps/MCMCSimple/MCMCSimple.py')
doctest.testfile('../apps/Phycas/Phycas.py')
#doctest.testfile('../apps/Polytomies/Polytomies.py')
#doctest.testfile('../apps/Simulator/Simulator.py')

raw_input('Press return to quit')
