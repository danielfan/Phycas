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
    phycas.sample_every = 20
    phycas.adapt_first = 2 # was 10
    phycas.verbose = True
    phycas.metropolis_weight = 300
    phycas.slice_weight = 1
    phycas.gg_do = True
    phycas.gg_nreps = 5 # was 1
    phycas.gg_kvect = [0.5, 1.0, 2.0]
    phycas.slice_max_units = 0
    phycas.use_inverse_shape = False

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
    phycas.treeFileOpen()

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
    commonSetup()

    # Change to HKY+Gamma model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(1.0, 1.0, 1.0, 1.0)
    phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNGammaRates(4)
    phycas.model.setPriorOnShapeInverse(phycas.use_inverse_shape)
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
    phycas.treeFileOpen()

    #commonSetup()
    phycas.run()

    hkyg_p = phycas.gg_Pm
    hkyg_g = phycas.gg_Gm
    hkyg_d = phycas.gg_Dm

def runHKYFLEX(rnseed):
    global hkyflex_p, hkyflex_g, hkyflex_d
    print
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY+FLEX analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    phycas.r.setSeed(rnseed)    
    commonSetup()

    # Change to HKY+FLEX model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(1.0)
    phycas.model.setNucleotideFreqs(1.0, 1.0, 1.0, 1.0)
    phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNGammaRates(4)
    phycas.model.setFlexModel()
    phycas.model.setFLEXProbParamPrior(ProbDist.ExponentialDist(1.0))
    phycas.model.setNotPinvarModel()

    phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))
    phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))
    phycas.likelihood.replaceModel(phycas.model)

    # Start with same tree we used for the simulation
    phycas.tree.buildFromString(phycas.tree_topology)
    phycas.likelihood.prepareForLikelihood(phycas.tree)

    # Open new parameter and tree files
    phycas.param_file_name = 'analHKYFLEX.nex.p'
    phycas.paramFileOpen()
    phycas.tree_file_name = 'analHKYFLEX.nex.t'
    phycas.treeFileOpen()

    phycas.run()

    hkyflex_p = phycas.gg_Pm
    hkyflex_g = phycas.gg_Gm
    hkyflex_d = phycas.gg_Dm

if __name__ == "__main__":
    phycas = Phycas()

    # Simulation settings
    master_seed = 15397  # was 13579
    phycas.r.setSeed(master_seed)
    num_sites = 2000    # was 1000

    # Define the names of the taxa to use when the simulated data set is saved to a file
    #POLPY_OLDWAY
    #phycas.taxon_labels = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']
    #POLPY_NEWWAY
    phycas.taxon_labels = ['long_1', 'short_2', 'long_3', 'short_4']

    # Create a model tree
    phycas.ntax = 4
    phycas.tree = Phylogeny.Tree()
    
    # Build a random tree
    #edge_len_dist = ProbDist.ExponentialDist(1.0)
    #edge_len_dist.setLot(phycas.r)
    #Phylogeny.TreeManip(phycas.tree).randomTree(
    #    phycas.ntax,     # number of tips
    #    phycas.r,        # pseudorandom number generator
    #    edge_len_dist,   # distribution from which to draw edge lengths
    #    False)           # Yule tree if True, edge lengths independent if False
    #self.starting_tree = self.tree.makeNewick()
    
    phycas.tree_topology = '(1:0.2,2:0.02,(3:0.2,4:0.02):0.02)'
    #phycas.tree_topology = '(1:0.1,2:0.1,(3:0.1,4:0.1):0.1)'
    #phycas.tree_topology = '(1:0.142857143,2:0.071428571,(3:0.142857143,4:0.071428571):0.071428571)'
    #phycas.tree_topology = '(1:0.192307692,2:0.038461538,(3:0.192307692,4:0.038461538):0.038461538)'
    #phycas.tree_topology = '(1:0.217391304,2:0.02173913,(3:0.217391304,4:0.02173913):0.02173913)'
    #phycas.tree_topology = '(1:0.23255814,2:0.011627907,(3:0.23255814,4:0.011627907):0.011627907)'

    phycas.tree.buildFromString(phycas.tree_topology)

    # Create a model
    phycas.model = Likelihood.HKYModel()
    sim_kappa = 4.0
    phycas.model.setKappa(sim_kappa)
    phycas.model.setNGammaRates(4)
    phycas.model.setNotPinvarModel()
    sim_shape = 0.2
    phycas.model.setShape(sim_shape)
    sim_piA = 0.2
    sim_piC = 0.3
    sim_piG = 0.3 
    sim_piT = 0.2
    phycas.model.setNucleotideFreqs(sim_piA, sim_piC, sim_piG, sim_piT)
    sim_model = phycas.model.getModelName()

    # Create a TreeLikelihood object to orchestrate simulations
    phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)

    # Prepare the tree for simulation (i.e. equip nodes with transition matrices)
    phycas.likelihood.prepareForSimulation(phycas.tree)

    # Create a SimData object to hold the simulated data
    sim_data = Likelihood.SimData()

    # Simulate num_sites of data and store in sim_data
    # Use the function simulateFirst (rather than just simulate) in order
    # to force calculation of transition probabilities
    phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)

    # Save simulated data to a NEXUS file using taxon_labels, datatype=dna and
    # using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
    sim_data.saveToNexusFile('simHKYg.nex', phycas.taxon_labels, 'dna', ('a','c','g','t'))

    # Copy the simulated data from sim_data to phycas.likelihood so that
    # we can run MCMC analyses on the simulated data
    phycas.likelihood.copyDataFromSimData(sim_data)
    phycas.nchar = num_sites # this should be set by copyDataFromSimData

    runHKY(master_seed)
    runHKYg(master_seed)
    runHKYFLEX(master_seed)
    
    # Output results
    outf = file('ggout.txt','w')

    outf.write('Simulation conditions:\n')
    outf.write('  no. sites   : %d\n' % num_sites)
    outf.write('  model tree  : %s\n' % phycas.tree_topology)
    outf.write('  subst. model: %s\n' % sim_model)
    outf.write('  kappa       : %.1f\n' % sim_kappa)
    outf.write('  shape       : %.1f\n' % sim_shape)
    outf.write('  piA         : %.1f\n' % sim_piA)
    outf.write('  piC         : %.1f\n' % sim_piC)
    outf.write('  piG         : %.1f\n' % sim_piG)
    outf.write('  piT         : %.1f\n' % sim_piT)
    outf.write('\n')

    outf.write('\nGelfand-Ghosh calculation for HKY analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hky_g, hky_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hky_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for HKY+G analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hkyg_g, hkyg_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hkyg_p, G, D))

    outf.write('\nGelfand-Ghosh calculation for HKY+FLEX analysis:\n')
    outf.write('  %6s %12s %12s %12s\n' % ('k','Pm','Gm','Dm'))
    for k,G,D in zip(phycas.gg_kvect, hkyflex_g, hkyflex_d):
        outf.write('  %6.1f %12.5f %12.5f %12.5f\n' % (k, hkyflex_p, G, D))

    outf.write('\nPenalty rankings (best listed first):\n')
    plist = [('HKY', hky_p), ('HKYg',hkyg_p), ('HKYFLEX',hkyflex_p)]
    plist.sort(lambda x,y: cmp(x[1], y[1]))
    for name,value in plist:
        outf.write('%12s %12.5f\n' % (name, value))

    for i,k in enumerate(phycas.gg_kvect):
        outf.write('\nFor k = %.1f:\n' % k)
        
        outf.write('\n  Goodness-of-fit rankings (best listed first):\n')
        glist = [('HKY', hky_g[i]), ('HKYg',hkyg_g[i]), ('HKYFLEX',hkyflex_g[i])]
        glist.sort(lambda x,y: cmp(x[1], y[1]))
        for name,value in glist:
            outf.write('%12s %12.5f\n' % (name, value))

        outf.write('\n  Overall rankings (best listed first):\n')
        dlist = [('HKY', hky_d[i]), ('HKYg',hkyg_d[i]), ('HKYFLEX',hkyflex_d[i])]
        dlist.sort(lambda x,y: cmp(x[1], y[1]))
        for name,value in dlist:
            outf.write('%12s %12.5f\n' % (name, value))

    outf.close()

