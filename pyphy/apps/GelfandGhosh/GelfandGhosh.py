# This file demonstrates simulation of data, analysis under various models, and
# assessment of the performance of each model on the simulated data using the
# posterior predictive loss approach advocated by Gelfand and Ghosh:
#
# Gelfand, A. E., and S. Ghosh. 1998. Model choice: a minimum posterior predictive
# loss approach. Biometrika 85:1-11.
#
# Start near the bottom of this file with the __name__ == "__main___" conditional
# Note that Phycas.py is used to most of the actual work.

import Conversions
from ProbDist import Lot
from Phycas import *

def commonSetup():
    # Set up MCMC
    phycas.ncycles = 4000
    phycassample_every = 20
    phycas.adapt_first = 50
    phycas.verbose = True
    phycas.metropolis_weight = 300
    phycas.slice_weight = 1
    phycas.gg_do = True
    phycas.gg_kvect = [0.5, 1.0, 2.0]

def runJC(rnseed):
    global jc_p, jc_g, jc_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'JC analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    

    # Create JC model and specify priors
    phycas.model = Likelihood.JCModel()
    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open the parameter and tree files
    phycas.param_file_name = 'analJC.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analJC.nex.t'
    phycas.treeFileOpen(taxon_names)

    commonSetup()
    phycas.run()

    jc_p = phycas.gg_Pm
    jc_g = phycas.gg_Gm
    jc_d = phycas.gg_Dm

def runHKY(rnseed):
    global hky_p, hky_g, hky_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    

    # Change to HKY model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(1.0, 1.0, 1.0, 1.0)
    phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))

    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open new parameter and tree files
    phycas.param_file_name = 'analHKY.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analHKY.nex.t'
    phycas.treeFileOpen(taxon_names)

    commonSetup()
    phycas.run()

    hky_p = phycas.gg_Pm
    hky_g = phycas.gg_Gm
    hky_d = phycas.gg_Dm

def runHKYg(rnseed):
    global hkyg_p, hkyg_g, hkyg_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY+Gamma analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    

    # Change to HKY+Gamma model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(1.0, 1.0, 1.0, 1.0)
    phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(0.5)
    phycas.model.setDiscreteGammaShapePrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNotPinvarModel()

    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open new parameter and tree files
    phycas.param_file_name = 'analHKYg.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analHKYg.nex.t'
    phycas.treeFileOpen(taxon_names)

    commonSetup()
    phycas.run()

    hkyg_p = phycas.gg_Pm
    hkyg_g = phycas.gg_Gm
    hkyg_d = phycas.gg_Dm

def runHKYgi(rnseed):
    global hkygi_p, hkygi_g, hkygi_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY+G+I analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    

    # Change to HKY+G+I model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(1.0, 1.0, 1.0, 1.0)
    phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNGammaRates(4)
    phycas.model.setShape(0.5)
    phycas.model.setDiscreteGammaShapePrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setPinvarModel()
    phycas.model.setPinvar(0.5)
    phycas.model.setPinvarPrior(ProbDist.BetaDist(1.0, 1.0))

    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open new parameter and tree files
    phycas.param_file_name = 'analHKYgi.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analHKYgi.nex.t'
    phycas.treeFileOpen(taxon_names)

    commonSetup()
    phycas.run()

    hkygi_p = phycas.gg_Pm
    hkygi_g = phycas.gg_Gm
    hkygi_d = phycas.gg_Dm

def runGTR(rnseed):
    global gtr_p, gtr_g, gtr_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~~'
    print 'GTR analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    

    # Change to GTR model
    #phycas.using_gtr = True
    phycas.model = Likelihood.GTRModel()
    phycas.model.setRelRates([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    phycas.model.setRelRatePrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNGammaRates(1)
    phycas.model.setNotPinvarModel()

    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open new parameter and tree files
    phycas.param_file_name = 'analGTR.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analGTR.nex.t'
    phycas.treeFileOpen(taxon_names)

    commonSetup()
    phycas.run()

    gtr_p = phycas.gg_Pm
    gtr_g = phycas.gg_Gm
    gtr_d = phycas.gg_Dm

if __name__ == "__main__":
    phycas = Phycas()

    # Define the names of the taxa to use when the simulated data set is saved to a file
    taxon_names = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

    # Create a model tree
    phycas.ntax = 4
    phycas.tree = Phylogeny.Tree()
    phycas.tree_topology = '(1:0.1,2:0.1,(3:0.1,4:0.1):0.1)'
    phycas.tree.buildFromString(phycas.tree_topology)

    # Create a model
    phycas.model = Likelihood.JCModel()
    sim_model = phycas.model.getModelName()

    # Create a TreeLikelihood object to orchestrate simulations
    phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)

    # Prepare the tree for simulation (i.e. equip nodes with transition matrices)
    phycas.likelihood.prepareForSimulation(phycas.tree)

    # Simulation settings
    master_seed = 13579
    phycas.r.setSeed(master_seed)
    phycas.sim_nreps = 1 # ignored at present
    num_sites = 1000

    # Create a SimData object to hold the simulated data
    sim_data = Likelihood.SimData()

    # Simulate num_sites of data and store in sim_data
    # Use the function simulateFirst (rather than just simulate) in order
    # to force calculation of transition probabilities
    phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)

    # Save simulated data to a NEXUS file using taxon_names, datatype=dna and
    # using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
    sim_data.saveToNexusFile('simJC.nex', taxon_names, 'dna', ('a','c','g','t'))

    # Copy the simulated data from sim_data to phycas.likelihood so that
    # we can run MCMC analyses on the simulated data
    phycas.likelihood.copyDataFromSimData(sim_data)
    phycas.nchar = num_sites # this should be set by copyDataFromSimData

    runJC(master_seed)
    runHKY(master_seed)
    runHKYg(master_seed)
    runHKYgi(master_seed)
    runGTR(master_seed)
    
    # Output results
    outf = file('ggout.txt','w')

    outf.write('Simulation conditions:\n')
    outf.write('  no. sites   : %d\n' % num_sites)
    outf.write('  model tree  : %s\n' % phycas.tree_topology)
    outf.write('  subst. model: %s\n' % sim_model)
    outf.write('\n')

    outf.write('Gelfand-Ghosh calculation for JC analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, jc_g, jc_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, jc_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for HKY analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hky_g, hky_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hky_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for HKY+G analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hkyg_g, hkyg_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hkyg_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for HKY+G+I analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hkygi_g, hkygi_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hkygi_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for GTR analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, gtr_g, gtr_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, gtr_p, G, D))

    outf.write('\nPenalty rankings (best listed first):\n')
    plist = [('JC',jc_p), ('HKY', hky_p), ('HKYg',hkyg_p), ('HKYgi',hkygi_p), ('GTR', gtr_p)]
    plist.sort(lambda x,y: cmp(x[1], y[1]))
    for name,value in plist:
        outf.write('%12s %12.5f\n' % (name, value))

    for i,k in enumerate(phycas.gg_kvect):
        outf.write('\nFor k = %.1f:\n' % k)
        
        outf.write('\n  Goodness-of-fit rankings (best listed first):\n')
        glist = [('JC',jc_g[i]), ('HKY', hky_g[i]), ('HKYg',hkyg_g[i]), ('HKYgi',hkygi_g[i]), ('GTR', gtr_g[i])]
        glist.sort(lambda x,y: cmp(x[1], y[1]))
        for name,value in glist:
            outf.write('%12s %12.5f\n' % (name, value))

        outf.write('\n  Overall rankings (best listed first):\n')
        dlist = [('JC',jc_d[i]), ('HKY', hky_d[i]), ('HKYg',hkyg_d[i]), ('HKYgi',hkygi_d[i]), ('GTR', gtr_d[i])]
        dlist.sort(lambda x,y: cmp(x[1], y[1]))
        for name,value in dlist:
            outf.write('%12s %12.5f\n' % (name, value))

    outf.close()

