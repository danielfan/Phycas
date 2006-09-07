from _ProbDist import *

class NormalDist(NormalDistBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate normal probability distribution.

    """
    def __init__(self, mean, stddev):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the mean and standard deviation when creating a new NormalDist
        object. e.g.,

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.getMean()
        2.0
        >>> print d.getStdDev()
        3.0

        """
        NormalDistBase.__init__(self, mean, stddev)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the normal distribution is continuous.
        
        """
        return NormalDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Normal'
        
        """
        return NormalDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another NormalDist
        object identical to this one. e.g.,

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.__str__()
        NormalDist(2.00000, 3.00000)
        
        """
        return NormalDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another NormalDist
        object identical to this one. e.g.,

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.__repr__()
        NormalDist(2.00000, 3.00000)
        
        """
        return NormalDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = NormalDist(0.0, 1.0)
        >>> d2 = NormalDist(10.0, 1.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        -2.303677281009
        >>> print "%.12f" % d2.sample()
        9.991924830267
        
        In this example, both d1 and d2 use the same random number generator
        when their sample functions are called. Generating all random numbers
        from a single random number generator has the advantage of allowing
        an entire analysis to be recreated by simply specifying the original
        seed (here, 1357) to the master Lot object. If setLot is not used,
        each distribution object maintains its own random number generator,
        which is initialized using the system clock at the time the object is
        created. This behavior is convenient, but makes it very hard to
        replicate results.

        """
        return NormalDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        -7.221565598551
        >>> print "%.12f" % d.sample()
        4.095823211335
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        -7.221565598551

        """
        return NormalDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return NormalDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.getMean()
        2.0

        """
        return NormalDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.getVar()
        9.0

        """
        return NormalDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.getStdDev()
        3.0

        """
        return NormalDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 
        negative infinity up to the value x.

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(0, 3)
        >>> print d.getCDF(0.0)
        0.5

        """
        return NormalDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the normal distribution specified
        by this NormalDist object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(0, 3)
        >>> d.setSeed(97531)
        >>> samples = [d.sample() for i in range(3)]
        >>> for x in samples:
        ...    print round(x,6)
        2.151006
        -6.76991
        4.733003

        """
        return NormalDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(0, 3)
        >>> ln_pdf = d.getLnPDF(1.0)
        >>> print round(ln_pdf, 9)
        -1.961995266

        The normal density is

              1           / (x - mean)^2 \
        ------------- exp| -------------- |
        sd*sqrt(2*pi)     \    2*sd^2    /
        
        If mean = 0.0, sd = 3.0 and x = 1.0. The density is

            1        / 1.0  \
        -------- exp| -----  | = 0.140578
        7.519884     \ 18.0 /

        The natural log (i.e., logarithm base e) density is

        (x - mean)^2
        ------------ - log(sd) - 0.5*log(2*pi)
           2*sd^2

        and the log density for the above example is thus

        1.0/18.0 - log(3.0) - 0.5*log(2*3.1415927)
        = 0.055556 - 1.098612 - 0.918939
        = -1.961995
        
        """
        return NormalDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability density function at the supplied
        value x. Returns the natural logarithm of the relative density at x.
        Use this function if speed is important but normalization is not, say
        in MCMC calculations. Use getLnPDF instead if you need to have a
        correctly normalized density value (i.e. from a density function that
        integrates to 1.0)

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(0, 3)
        >>> print d.getRelativeLnPDF(1.5)
        0.125

        The normal density is

              1           / (x - mean)^2 \
        ------------- exp| -------------- |
        sd*sqrt(2*pi)     \    2*sd^2    /
        
        If mean = 0.0, sd = 3.0 and x = 1.5. The density is

            1        / 2.25 \
        -------- exp| -----  | = 0.150687
        7.519884     \ 18.0 /

        The natural log (i.e., logarithm base e) density is

        (x - mean)^2
        ------------ - log(sd) - 0.5*log(2*pi)
           2*sd^2

        For the above example we have
        2.25/18.0 - log(3) - 0.5*log(2*3.1414927)
        = 0.125 - 1.098612 - 0.918939 = -1.892551

        This function returns only the first term (0.125) because the second
        and third terms involve only the standard deviation (a constant) and
        not the value x.
        
        """
        return NormalDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. 

        >>> from phypy.ProbDist import *
        >>> d = NormalDist(2, 3)
        >>> print d.getMean()
        2.0
        >>> print d.getVar()
        9.0
        >>> d.setMeanAndVariance(5, 4)
        >>> print d.getMean()
        5.0
        >>> print d.getVar()
        4.0
        
        """
        return NormalDistBase.setMeanAndVariance(self, mean, var)
    
