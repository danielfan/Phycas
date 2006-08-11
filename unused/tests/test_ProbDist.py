from ProbDist import *
import unittest,math,sys

def pause():
    raw_input('\nPress return to continue...')

def almostEqual(x, y, tol=1.e-8):
    diff = math.fabs(x-y)
    if diff <= tol:
        return True
    else:
        return False

class ProbDistTest(unittest.TestCase):
    def testBernoulli(self):
        b = BernoulliDist(.2)
        self.assertEqual(b.getMean(), 0.2)
        self.assertEqual(almostEqual(b.getVar(), 0.16), True)
    def testExponential(self):
        d = ExponentialDist(1.0)
        self.assertEqual(almostEqual(b.getMean(), 1.0), True)

def exerciseDist(dist, lot, test_values):
    print '\nExercising %s:' % dist.getDistName()
    #if dist.isDiscrete():
    print '  Description ', dist.getDistDescription()
    print '  Mean        ', dist.getMean()
    print '  Variance    ', dist.getVar()
    print '  Std. dev.   ', dist.getStdDev()
    x = [dist.sample(lot) for i in range(1000)]
    s = sum(x)
    ss = sum([y*y for y in x])
    m = float(s)/1000.0
    v = (ss - 1000.0*m*m)/1000.0
    print '  Mean of sample of 1000 values:',float(s)/1000.0
    print '  Variance of same sample:',v
    for v in test_values:
        if dist.isDiscrete():
            print '  Pr(%.0f) = %f (cum. prob. = %f)' % (v,math.exp(dist.getLnPDF(v)),dist.getCDF(v))
        else:
            print '  density at',v,'=',math.exp(dist.getLnPDF(v))

if __name__=='__main__':
    lot = Lot()
    print 'Starting pseudorandom number seed was',lot.getInitSeed()
    exerciseDist(BernoulliDist(0.2), lot, (0, 1))
    exerciseDist(BinomialDist(10, 0.6), lot, (0, 1, 5))
    exerciseDist(UniformDist(0.0, 2.0), lot, (0.1, 1.0, 1.9))
    exerciseDist(GammaDist(0.2, 5.0), lot, (0.1, 1.0, 10.0))
    exerciseDist(ExponentialDist(0.1), lot, (0.1, 1.0, 10.0))
    try:
        print "\nTrying to create an InverseGammaDist with shape < 2..."
        exerciseDist(InverseGammaDist(1.0, 0.5), lot, ())
    except ValueError, e:
        print 'Oops!:',e
    exerciseDist(InverseGammaDist(3.0, 0.5), lot, (0.1, 1.0, 10.0))
    pause()
    #unittest.main()
    
