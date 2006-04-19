from _LikelihoodBase import *

class JCModel(JCModelBase):    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|    """    Encapsulates the Jukes and Cantor (1969) substitution model, which
    assumes base frequencies are equal and all types of substitutions
    occur at the same rate.
        Literature Cited:
    
    Jukes, T. H., and C. R. Cantor. 1969. Evolution of protein molecules.
    Pages 21-132 in Mammalian Protein Metabolism (H. N. Munro, ed.)
    Academic Press, New York.    

    """
    def getNStates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the number of states (always 4 for this model).
        
        """
        return JCModelBase.getNStates(self)
    def getStateFreqs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a tuple comprising the 4 state frequencies. Always
        (0.25, 0.25, 0.25, 0.25) for this model.

        >>> import Likelihood
        >>> model = Likelihood.JCModel()
        >>> print model.getStateFreqs()
        (0.25, 0.25, 0.25, 0.25)
        
        """
        return JCModelBase.getStateFreqs(self)
    def setAllFreqsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all four state frequencies to 0.25. Superfluous for this model,
        but included for conformity
        
        """
        return JCModelBase.setAllFreqsEqual(self)
    def getNGammaRates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the current number of relative rate categories.
        
        """
        return JCModelBase.getNGammaRates(self)
    def setNGammaRates(self, n):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return JCModelBase.setNGammaRates(self, n)
    def getRateProbs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return JCModelBase.getRateProbs(self)
    def setAllRateProbsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return JCModelBase.setAllRateProbsEqual(self)

    def setPriorOnShapeInverse(self, invert):            #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        If True is specified, then the gamma shape parameter will actually
        update the inverse of the shape parameter rather than the shape
        parameter itself.
        
        """
        JCModelBase.setPriorOnShapeInverse(self, invert)

class HKYModel(HKYModelBase):    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|    """    Encapsulates the Hasegawa-Kishino-Yano (1985) substitution model,
    which allows unequal base frequencies and transition-type
    substitutions to occur at a different rate than transversion-type
    substitutions. The transition/transversion rate ratio is termed kappa,
    whereas the probability of any transition divided by the probability
    of any transversion is known as tratio. This model can be made equal
    to the Felsenstein (1981) model by setting kappa to 1.0, and to the
    Kimura (1980) 2-parameter model by making the base frequencies all
    0.25. This model is identical to the Jukes-Cantor (1969) model if
    the base frequencies are equal and kappa is 1.0.

    Literature Cited:
    
    Felsenstein, J. 1981. Evolutionary trees from DNA sequences:  a
    maximum likelihood approach. Journal of Molecular Evolution
    17: 368-376.

    Hasegawa, M., H. Kishino, and T. Yano. 1985. Dating of the human-ape
    splitting by a molecular clock of mitochondrial DNA. Journal of
    Molecular Evolution 22: 160-174.

    Jukes, T. H., and C. R. Cantor. 1969. Evolution of protein molecules.
    Pages 21-132 in Mammalian Protein Metabolism (H. N. Munro, ed.)
    Academic Press, New York.

    Kimura, M. 1980. A simple method for estimating evolutionary rate of
    base substitutions through comparative studies of nucleotide sequences.
    Journal of Molecular Evolution 16:111-120.    
        """
    def getNStates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the number of states (always 4 for this model).
        
        """
        return HKYModelBase.getNStates(self)
    def getStateFreqs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a tuple comprising the 4 state frequencies.
        
        """
        return HKYModelBase.getStateFreqs(self)
    def setStateFreqParam(self, i, value):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets frequency parameter for state i to value. use i = 0 for A,
        i = 1 for C, i = 2 for G and i = 3 for T. Note that value can be any
        non-negative number; there is no need to ensure that it is between
        0.0 and 1.0 (although there is nothing wrong with providing normalized
        frequencies). The four frequency parameters are normalized for use in
        all calculations involving base frequencies. Thus, specifying 1, 2, 3,
        and 4 for the four frequency parameters will result in the relative
        base frequencies being set to 0.1, 0.2, 0.3 and 0.4.
        
        """
        return HKYModelBase.setStateFreqParam(self, i, value)
    def setAllFreqsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all four state frequencies to 0.25.
        
        """
        return HKYModelBase.setAllFreqsEqual(self)
    def setNucleotideFreqs(self, freqA, freqC, freqG, freqT):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the four state frequencies to the values provided, which should
        all be greater than or equal to 0.0.
        
        """
        assert freqA >= 0.0 and freqC >= 0.0 and freqG >= 0.0 and freqT >= 0.0
        return HKYModelBase.setNucleotideFreqs(self, freqA, freqC, freqG, freqT)
    def getKappa(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the current value of the transition-transversion rate
        ratio kappa.
        
        """
        return HKYModelBase.getKappa(self)
    def setKappa(self, k):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the value of the transition-transversion rate ratio kappa.
        
        """
        return HKYModelBase.setKappa(self, k)
    def setKappaFromTRatio(self, tratio):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the value of the transition-transversion rate ratio kappa given
        a supplied value of tratio, the transition-transversion ratio. The
        rate ratio kappa is related to tratio as follows (where piA means the
        frequency of base A, piC means the frequency of base C, etc.
                tratio (piA + piG) (piC + piT)
        kappa = -------------------------------
                    (piA piG + piC piT)

        Thus, if piA = piC = piG = piT = 0.25, kappa is twice the tratio
        because there are twice as many kinds of transversion-type
        substitutions compared to transition-type substitutions.
        """
        return HKYModelBase.setKappaFromTRatio(self, tratio)
    def calcTRatio(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Calculates the transition/transversion ratio (tratio) given the
        transition/transversion rate ratio (kappa) and the relative base
        frequencies. Here are the details of the calculation (for brevity,
        piA symbolizes the frquency of base A, piC the frequency of base C,
        etc.
        
        Parameters: b = transversion rate, k = kappa, dt = infinitesimal time
        The probability that a transition from base A to base G occurs over
        time dt is the probability of starting in state A (piA) times the
        probability of a transition to base G (piG k b dt).
        
        Pr(any transition | dt) = (piA piG k b dt) + (piC piT k b dt)
          + (piG piA k b dt) + (piT piC k b dt)
          = 2 k b dt (piA piG + piC piT)
        
        Pr(any transversion | dt) = (piA piC b dt) + (piA piT b dt)
          + (piC piA b dt) + (piC piG b dt) + (piG piC b dt) + (piG piT b dt)
          + (piT piA b dt) + (piT piG b dt)
          = 2 b dt (piA + piG) (piC + piT)

                 2 k b dt (piA piG + piC piT)     k (piA piG + piC piT)
        tratio = ------------------------------ = -----------------------
                 2 b dt (piA + piG) (piC + piT)   (piA + piG) (piC + piT)
        
        """
        return HKYModelBase.calcTRatio(self)
    def getNGammaRates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the current number of relative rate categories.
        
        """
        return HKYModelBase.getNGammaRates(self)
    def setNGammaRates(self, n):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return HKYModelBase.setNGammaRates(self, n)
    def getRateProbs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return HKYModelBase.getRateProbs(self)
    def setAllRateProbsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return HKYModelBase.setAllRateProbsEqual(self)
    def setPriorOnShapeInverse(self, invert):            #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        If True is specified, then the gamma shape parameter will actually
        update the inverse of the shape parameter rather than the shape
        parameter itself.
        
        """
        HKYModelBase.setPriorOnShapeInverse(self, invert)

class GTRModel(GTRModelBase):    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|    """    Encapsulates the General Time Reversible substitution model, first
    described by Rodriguez et al. (1990). This model allows unequal base
    frequencies as well as six different relative rates corresponding to
    the substitution classes A <-> C, A <-> G, A <-> T, C <-> G, C <-> T
    and G <-> T. Constrained versions of this model can be made equal to
    the HKY85, K80, F81 and JC models.

    Literature Cited:

    Rodrigues, F., J. L. Oliver, A. Marin, and J. R. Medina. 1990. The
    general stochastic model of nucleotide substitution. Journal of
    Theoretical Biology 142: 485-501.
    
    """
    def getNStates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the number of states (always 4 for this model).
        
        """
        return GTRModelBase.getNStates(self)
    def getStateFreqs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a tuple comprising the 4 state frequencies.
        
        """
        return GTRModelBase.getStateFreqs(self)
    def setStateFreqParam(self, i, value):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets frequency parameter for state i to value. use i = 0 for A,
        i = 1 for C, i = 2 for G and i = 3 for T. Note that value can be any
        non-negative number; there is no need to ensure that it is between
        0.0 and 1.0 (although there is nothing wrong with providing normalized
        frequencies). The four frequency parameters are normalized for use in
        all calculations involving base frequencies. Thus, specifying 1, 2, 3,
        and 4 for the four frequency parameters will result in the relative
        base frequencies being set to 0.1, 0.2, 0.3 and 0.4.
        
        """
        return GTRModelBase.setStateFreqParam(self, i, value)
    def setAllFreqsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all four state frequencies to 0.25.
        
        """
        return GTRModelBase.setAllFreqsEqual(self)
    def setNucleotideFreqs(self, freqA, freqC, freqG, freqT):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the four state frequencies to the values provided, which should
        all be greater than or equal to 0.0.
        
        """
        assert freqA >= 0.0 and freqC >= 0.0 and freqG >= 0.0 and freqT >= 0.0
        return GTRModelBase.setNucleotideFreqs(self, freqA, freqC, freqG, freqT)
    def getRelRates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a list comprising the six relative rates.
        
        """
        return GTRModelBase.getRelRates(self)
    def setRelRates(self, rr):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all six relative rates.
        
        """
        return GTRModelBase.setRelRates(self, rr)
    def calcTRatio(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Calculates the transition/transversion ratio (tratio) given the six
        relative rates and the relative base frequencies. Here are the
        details of the calculation.

        Symbols:
            rAC = rate of A <-> C   piA = frequency of base A
            rAG = rate of A <-> G   piC = frequency of base C
            rAT = rate of A <-> T   piG = frequency of base G
            rCG = rate of C <-> G   piT = frequency of base T
            rCT = rate of C <-> T
            rGT = rate of G <-> T   dt = infinitesimal time
            
        The probability that a transition from base A to base G occurs over
        time dt is the probability of starting in state A (piA) times the
        probability of a transition to base G (piG rAG dt).
        
        Pr(any transition | dt) = (piA piG rAG dt) + (piC piT rCT dt)
          + (piG piA rAG dt) + (piT piC rCT dt)
          = 2 dt (piA piG rAG + piC piT rCT)
        
        Pr(any transversion | dt) = (piA piC rAC dt) + (piA piT rAT dt)
          + (piC piA rAC dt) + (piC piG rCG dt) + (piG piC rCG dt)
          + (piG piT rGT dt) + (piT piA rAT dt) + (piT piG rGT dt)
          = 2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)

                            2 dt (piA piG rAG + piC piT rCT)
        tratio = ------------------------------------------------------------
                 2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)
        
                            piA piG rAG + piC piT rCT
               = -----------------------------------------------------
                 piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT

        Example:
        >>> import Likelihood
        >>> model = Likelihood.GTRModel()
        >>> model.setRelRates([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
        >>> model.setNucleotideFreqs(0.25, 0.25, 0.25, 0.25)
        >>> print model.calcTRatio()
        2.0
        
        """
        return GTRModelBase.calcTRatio(self)
    def getNGammaRates(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns the current number of relative rate categories.
        
        """
        return GTRModelBase.getNGammaRates(self)
    def setNGammaRates(self, n):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return GTRModelBase.setNGammaRates(self, n)
    def getRateProbs(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return GTRModelBase.getRateProbs(self)
    def setAllRateProbsEqual(self):        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return GTRModelBase.setAllRateProbsEqual(self)
    def setPriorOnShapeInverse(self, invert):            #---+----|----+----|----+----|----+----|----+----|----+----|----+----|        """
        If True is specified, then the gamma shape parameter will actually
        update the inverse of the shape parameter rather than the shape
        parameter itself.
        
        """
        GTRModelBase.setPriorOnShapeInverse(self, invert)
