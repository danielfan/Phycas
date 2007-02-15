# This file demonstrates simulation of data, analysis under various models, and
# assessment of the performance of each model on the simulated data using the
# posterior predictive loss approach advocated by Gelfand and Ghosh:
#
# Gelfand, A. E., and S. Ghosh. 1998. Model choice: a minimum posterior predictive
# loss approach. Biometrika 85:1-11.
#
# Start near the bottom of this file with the __name__ == "__main___" conditional
# Note that Phycas.py is used to most of the actual work.

from phycas import *

def commonSetup():
    # Set up MCMC
    phycas.ncycles = 200
    phycas.sample_every = 2 
    phycas.adapt_first = 10
    phycas.verbose = True
    phycas.ls_move_weight = 100
    phycas.slice_weight = 1
    phycas.gg_do = True
    phycas.gg_nreps = 5 # was 1
    phycas.gg_kvect = [0.5, 1.0, 2.0]
    phycas.slice_max_units = 0
    phycas.use_inverse_shape = False
    phycas.using_hyperprior = True
    phycas.gg_outfile = None
    phycas.starting_tree_source = 'usertree'
    phycas.starting_tree = '(1:0.2,2:0.02,(3:0.2,4:0.02):0.02)'
    phycas.default_model = 'hky'
    phycas.data_source = 'memory'
    phycas.estimate_pinvar = False
    phycas.use_flex_model = False

def runHKY(rnseed):
    print
    print '~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~'

    commonSetup()

    phycas.random_seed = rnseed
    phycas.data_file_name = 'analHKY.nex'

    phycas.setup()
    phycas.run()

    global hky_p, hky_g, hky_d
    hky_p = phycas.gg_Pm
    hky_g = phycas.gg_Gm
    hky_d = phycas.gg_Dm

def runHKYg(rnseed):
    print
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY+Gamma analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    commonSetup()

    phycas.random_seed = rnseed
    phycas.data_file_name = 'analHKYg.nex'
    phycas.num_rates = 4

    phycas.setup()
    phycas.run()

    global hkyg_p, hkyg_g, hkyg_d
    hkyg_p = phycas.gg_Pm
    hkyg_g = phycas.gg_Gm
    hkyg_d = phycas.gg_Dm

def runHKYFLEX(rnseed):
    print
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'HKY+FLEX analysis beginning'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    commonSetup()

    phycas.random_seed            = rnseed
    phycas.data_file_name         = 'analHKYflex.nex'
    phycas.num_rates              = 1
    phycas.use_flex_model         = True
    phycas.flex_ncat_move_weight  = 1
    phycas.flex_num_spacers       = 1
    phycas.flex_phi               = 0.25
    phycas.flex_L                 = 1.0
    phycas.flex_lambda            = 1.0 
    phycas.flex_prob_param_prior  = ProbDist.ExponentialDist(1.0)

    phycas.setup()
    phycas.run()

    global hkyflex_p, hkyflex_g, hkyflex_d
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
    phycas.taxon_labels = ['long_1', 'short_2', 'long_3', 'short_4']

    # Create a model tree
    phycas.ntax = 4
    phycas.tree = Phylogeny.Tree()
    phycas.tree_topology = '(1:0.2,2:0.02,(3:0.2,4:0.02):0.02)'
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

