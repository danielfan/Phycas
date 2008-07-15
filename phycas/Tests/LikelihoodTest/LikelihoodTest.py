# This example is designed to check the likelihood calculation under most models
# supported by Phycas. A data set is simulated under the most complex model, and
# analyzed under a spectrum of simpler models. The data set is saved as a nexus file
# complete with PAUP blocks that allow verification of Phycas's likelihood
# calculations by PAUP. A second sweep of models is done for a real data set
# (nyldna4.nex) and a paup command file is written to allow verification of Phycas'
# likelihood calculations.

from phycas import *

def tryAllModels(fn):
    # Create string containing PAUP commands that will be added to the end of the
    # file fn to check results - we will add more commands to this string as we go
    paup_commands = []
    #paup_commands.append('\n[!\n***** HKY+G+I (estimate everything) *****]')
    paup_commands.append('\n[!\n***** GTR+G+I (estimate everything) *****]')
    #paup_commands.append('lset nst=2 variant=hky basefreq=estimate tratio=estimate rates=gamma shape=estimate pinvar=estimate;')
    paup_commands.append('lset nst=6 basefreq=estimate rmatrix=estimate rates=gamma shape=estimate pinvar=estimate;')
    paup_commands.append('lscores 1 / userbrlen;')

    print
    print '************* Testing GTRModel *******************'

    # Compute likelihood using the GTR+G+I model
    print '\nGTR+G+I model'
    phycas.model = Likelihood.GTRModel()
    phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    ref_lnL = lnL
    print 'lnL = %.5f (this is the reference lnL)' % (lnL)
    paup_commands.append('\n[!\n***** GTR+G+I (using GTRModel) *****]')
    paup_commands.append('lset nst=6 basefreq=(0.1 0.2 0.3) rmatrix=(1.8 4.0 1.5 1.2 5.0) rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas GTR+G+I lnL = %.5f]' % lnL)

    # Compute likelihood using the GTR+I model
    print '\nGTR+I model'
    phycas.model = Likelihood.GTRModel()
    phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    #ref_lnL = lnL
    print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** GTR+I (using GTRModel) *****]')
    paup_commands.append('lset nst=6 basefreq=(0.1 0.2 0.3) rmatrix=(1.8 4.0 1.5 1.2 5.0) rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas GTR+I lnL = %.5f]' % lnL)

    # Compute likelihood using the GTR+G model
    print '\nGTR+G model'
    phycas.model = Likelihood.GTRModel()
    phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** GTR+G (using GTRModel) *****]')
    paup_commands.append('lset nst=6 basefreq=(0.1 0.2 0.3) rmatrix=(1.8 4.0 1.5 1.2 5.0) rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas GTR+G lnL = %.5f]' % lnL)

    # ************* temporary below here **************
    #raw_input('debug stop: just about to compute GTR+psr')

    #print '\nGTR+psr model'
    #phycas.model = Likelihood.GTRModel()
    #phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    #phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    #phycas.model.setNGammaRates(4)
    #phycas.model.setShape(1.2)
    #phycas.model.setNotPinvarModel()
    #phycas.likelihood.usePatternSpecificRates() 
    #phycas.likelihood.replaceModel(phycas.model)
    #phycas.likelihood.prepareForLikelihood(phycas.tree)
    #lnL = phycas.likelihood.calcLnL(phycas.tree)
    #print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    #paup_commands.append('\n[!\n***** GTR+psr (actually, using GTR+G since no way to do psr in PAUP*) *****]')
    #paup_commands.append('lset nst=6 basefreq=(0.1 0.2 0.3) rmatrix=(1.8 4.0 1.5 1.2 5.0) rates=gamma shape=1.2 pinvar=0.0;')
    #paup_commands.append('lscores 1 / userbrlen;')
    #paup_commands.append('[!Phycas GTR+G lnL = %.5f]' % lnL)
    
    #phycas.likelihood.doNotUsePatternSpecificRates()
    # ************* temporary above here **************

    # Compute likelihood using the GTR model
    print '\nGTR model'
    phycas.model = Likelihood.GTRModel()
    phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** GTR (using GTRModel) *****]')
    paup_commands.append('lset nst=6 basefreq=(0.1 0.2 0.3) rmatrix=(1.8 4.0 1.5 1.2 5.0) rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas GTR lnL = %.5f]' % lnL)

    print
    print '************* Testing HKYModel *******************'

    # Compute likelihood using the HKY+G+I model
    print '\nHKY+G+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    ref_lnL = lnL
    print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** HKY+G+I (using HKYModel) *****]')
    paup_commands.append('lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas HKY+G+I lnL = %.5f]' % lnL)

    # Compute likelihood using the HKY+I model
    print '\nHKY+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    ref_lnL = lnL
    print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** HKY+I (using HKYModel) *****]')
    paup_commands.append('lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas HKY+I lnL = %.5f]' % lnL)

    # Compute likelihood using the HKY+G model
    print '\nHKY+G model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than the reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** HKY+G (using HKYModel) *****]')
    paup_commands.append('lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas HKY+G lnL = %.5f]' % lnL)

    # Compute likelihood using the HKY model
    print '\nHKY model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('\n[!\n***** HKY (using HKYModel) *****]')
    paup_commands.append('lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas HKY lnL = %.5f]' % lnL)

    # Compute likelihood using the F81+G+I model
    print '\nF81+G+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** F81+G+I (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=(0.1 0.2 0.3) rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas F81+G+I lnL = %.5f]' % lnL)

    # Compute likelihood using the F81+I model
    print '\nF81+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** F81+I (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=(0.1 0.2 0.3) rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas F81+I lnL = %.5f]' % lnL)

    # Compute likelihood using the F81+G model
    print '\nF81+G model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** F81+G (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=(0.1 0.2 0.3) rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas F81+G lnL = %.5f]' % lnL)

    # Compute likelihood using the F81 model
    print '\nF81 model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** F81 (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=(0.1 0.2 0.3) rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas F81 lnL = %.5f]' % lnL)

    # Compute likelihood using the K80+G+I model
    print '\nK80+G+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** K80+G+I (using HKYModel) *****]')
    paup_commands.append('lset nst=2 basefreq=equal tratio=2.0 rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas K80+G+I lnL = %.5f]' % lnL)

    # Compute likelihood using the K80+I model
    print '\nK80+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** K80+I (using HKYModel) *****]')
    paup_commands.append('lset nst=2 basefreq=equal tratio=2.0 rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas K80+I lnL = %.5f]' % lnL)

    # Compute likelihood using the K80+G model
    print '\nK80+G model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** K80+G (using HKYModel) *****]')
    paup_commands.append('lset nst=2 basefreq=equal tratio=2.0 rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas K80+G lnL = %.5f]' % lnL)

    # Compute likelihood using the K80 model
    print '\nK80 model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** K80 (using HKYModel) *****]')
    paup_commands.append('lset nst=2 basefreq=equal tratio=2.0 rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas K80 lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+G+I model
    print '\nJC+G+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+G+I (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+G+I lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+I model
    print '\nJC+I model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+I (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+I lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+G model
    print '\nJC+G model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+G (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+G lnL = %.5f]' % lnL)

    # Compute likelihood using the JC model
    print '\nJC model'
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setAllFreqsEqual()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC (using HKYModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC lnL = %.5f]' % lnL)

    print
    print '************** Testing JCModel *******************'

    # Compute likelihood using the JC model
    print '\nJC model'
    phycas.model = Likelihood.JCModel()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC (using JCModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=equal pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+G model
    print '\nJC+G model'
    phycas.model = Likelihood.JCModel()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNotPinvarModel()
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+G (using JCModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=gamma shape=1.2 pinvar=0.0;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+G lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+I model
    print '\nJC+I model'
    phycas.model = Likelihood.JCModel()
    phycas.model.setNGammaRates(1)
    #phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+I (using JCModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=equal pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+I lnL = %.5f]' % lnL)

    # Compute likelihood using the JC+G+I model
    print '\nJC+G+I model'
    phycas.model = Likelihood.JCModel()
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)
    phycas.likelihood.replaceModel(phycas.model)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    lnL = phycas.likelihood.calcLnL(phycas.tree)
    print 'lnL = %.5f (%.5f worse than reference lnL)' % (lnL, ref_lnL - lnL)
    paup_commands.append('[!\n***** JC+G+I (using JCModel) *****]')
    paup_commands.append('lset nst=1 basefreq=equal rates=gamma shape=1.2 pinvar=0.3;')
    paup_commands.append('lscores 1 / userbrlen;')
    paup_commands.append('[!Phycas JC+G+I lnL = %.5f]' % lnL)

    # Add a PAUP block to the file named fn to make it easy to check the results
    f = file(fn, 'a')
    f.write('\n')
    f.write('\nbegin paup;')
    f.write('\n  set criterion=likelihood storebrlen;')
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin trees;')
    f.write('\n  translate')
    for i,nm in enumerate(phycas.taxon_names):
        if nm.count(' ') > 0:
            f.write("\n    %d '%s'" % (i, nm))
        else:
            f.write("\n    %d %s" % (i, nm))
        if i < len(phycas.taxon_names) - 1:
            f.write(',')
        else:
            f.write(';')
    f.write('\n  utree simtree = %s;' % model_tree)
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin paup;\n')
    f.write('\n'.join(paup_commands))
    f.write('\nend;')
    f.write('\n')
    f.close()

def simulateData(fn):
    # Define the names of the taxa to use when the simulated data set is saved to a file
    phycas.taxon_names = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

    # Create a simulation model
    #phycas.model = Likelihood.HKYModel()
    phycas.model = Likelihood.GTRModel()
    #phycas.model.setKappa(4.0)
    phycas.model.setRelRates([1.8, 4.0, 1.5, 1.2, 5.0, 1.0])
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(1.2)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.3)

    # Create a likelihood object to orchestrate both simulations and likelihood calculations
    phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)

    # Prepare the tree for simulation (i.e. equip nodes with transition matrices)
    phycas.likelihood.prepareForSimulation(phycas.tree)

    # Simulation settings
    phycas.r.setSeed(13579)
    phycas.sim_nreps = 1 # ignored at present
    phycas.sim_outfile = 'simout.nex'
    #num_sites = 5000
    num_sites = 100000

    # Create a SimData object to hold the simulated data
    sim_data = Likelihood.SimData()

    # Simulate num_sites of data and store in sim_data
    # Use the function simulateFirst (rather than just simulate) in order
    # to force calculation of transition probabilities
    phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)

    # Save simulated data to a NEXUS file using taxon_names, datatype=dna and
    # using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
    sim_data.saveToNexusFile('simulated.nex', phycas.taxon_names, 'dna', ('a','c','g','t'))

    # Copy the simulated data from sim_data to phycas.likelihood so that
    # we can compute the likelihood for the simulated data
    phycas.likelihood.copyDataFromSimData(sim_data)    

def readData(fn):
    phycas.reader = ReadNexus.NexusReader()
    phycas.reader.readFile(fn)
    phycas.data_matrix = phycas.reader.getLastDiscreteMatrix(True)
    phycas.ntax = phycas.data_matrix.getNTax()
    phycas.nchar = phycas.data_matrix.getNChar()
    phycas.taxon_names = phycas.reader.getTaxLabels()
    #print 'ntax =',phycas.ntax
    #print 'nchar =',phycas.nchar
    #print 'taxon names =',phycas.taxon_names
    phycas.model = Likelihood.JCModel()
    phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)
    phycas.likelihood.copyDataFromDiscreteMatrix(phycas.data_matrix)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

def createCommandFile(fn, dataf):
    outf = file(fn, 'w')
    outf.write('#nexus\n\n')
    outf.write('begin paup;\n')
    outf.write("  set nowarnroot;\n")
    outf.write("  exe '%s';\n" % dataf)
    outf.write('end;\n')
    outf.close()    

if __name__ == '__main__':    
    phycas = Phycas()

    # Create a model tree
    phycas.tree = Phylogeny.Tree()
    model_tree = '(0:0.1,1:0.15,(2:0.025,3:0.15):0.05)'
    phycas.tree.buildFromString(model_tree, True)   # 2nd arg. is zero_based_tips

    #m = Lock()
    #tv = TreeViewer(tree=phycas.tree, mutex=m)
    #tv.start()

    print
    print '+------------------------------------------------+'
    print '|           Analyzing nyldna4.nex                |'
    print '+------------------------------------------------+'
    dataf = '../Data/nyldna4-compressed.nex'
    #dataf = '../Data/nyldna4.nex'
    readData(dataf)
    createCommandFile('check.nex', dataf)
    tryAllModels('check.nex')

    print
    print '+------------------------------------------------+'
    print '|          Analyzing Simulated Data              |'
    print '+------------------------------------------------+'

    simulateData('simulated.nex')
    tryAllModels('simulated.nex')

