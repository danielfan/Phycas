# This example program simulates data on a 5-taxon tree containing one polytomy
# This data set is then analyzed using an MCMC analysis that allows polytomies
# to demonstrate that the polytomous true tree can be recovered

import Conversions
from Phycas import *

phycas = Phycas()

print
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Simulating a data set'
print '~~~~~~~~~~~~~~~~~~~~~'

# Define the names of the taxa to use when the simulated data set is saved to a file
taxon_names = ['P._fimbriata', 'P._parksii', 'P._articulata', 'P._gracilis', 'P._macrophylla']

# Create a model tree containing one polytomy
phycas.ntax = 5
phycas.tree = Phylogeny.Tree()
phycas.tree_topology = '(1:0.1,2:0.1,(3:0.1,4:0.1,5:0.1):0.1)'
phycas.tree.buildFromString(phycas.tree_topology)

# Create a model
phycas.model_type = 'hky'
phycas.model = Likelihood.HKYModel()
phycas.model.setKappa(4.0)
phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))

# Create a likelihood object to orchestrate both simulations and likelihood calculations
phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)

# Prepare the tree for simulation (i.e. equip nodes with transition matrices)
phycas.likelihood.prepareForSimulation(phycas.tree)

# Simulation settings
phycas.r.setSeed(461019)
phycas.sim_nreps = 1 # ignored at present
num_sites = 5000

# Create a SimData object to hold the simulated data
sim_data = Likelihood.SimData()

# Simulate num_sites of data and store in sim_data
# Use the function simulateFirst (rather than just simulate) in order
# to force calculation of transition probabilities
phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)

# Save simulated data to a NEXUS file using taxon_names, datatype=dna and
# using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
sim_data.saveToNexusFile('simHKY.nex', taxon_names, 'dna', ('a','c','g','t'))

# Add a MrBayes block to make it easier to summariz trees later
# A temporary measure: this functionality should be built into Phycas
dataf = file('simHKY.nex', 'a')
dataf.write('\n\nbegin mrbayes;')
dataf.write('\n  sumt file=analHKY.nex nrun=1;')
dataf.write('\nend;\n')
dataf.close()

# Copy the simulated data from sim_data to phycas.likelihood so that
# we can run MCMC analyses on the simulated data
phycas.likelihood.copyDataFromSimData(sim_data)
phycas.nchar = num_sites # this should be set by copyDataFromSimData

print
print '~~~~~~~~~~~~~~~~~~~~~~'
print 'HKY analysis beginning'
print '~~~~~~~~~~~~~~~~~~~~~~'

# Finish setting up the model for MCMC
d = ProbDist.ExponentialDist(10.0)
d.setLot(phycas.r)
phycas.model.setEdgeLenPrior(d)
d = ProbDist.InverseGammaDist(2.1, 0.9090909)
d.setLot(phycas.r)
phycas.model.setEdgeLenHyperPrior(d)
phycas.likelihood.replaceModel(phycas.model)

# Tell phycas we want to allow polytomies
phycas.allow_polytomies = True
phycas.polytomy_prior   = True
phycas.topo_prior_C     = 2.0

# Open new parameter and tree files
phycas.param_file_name = 'analHKY.nex.p'
phycas.paramFileOpen()
phycas.tree_file_name = 'analHKY.nex.t'
phycas.treeFileOpen(taxon_names)

# Build a random starting tree for an MCMC analysis and prepare it
# for likelihood calculations (i.e. equip nodes with both transition
# matrices and conditional likelihood arrays)
phycas.starting_tree_source = 'random'
phycas.setupTree()
phycas.likelihood.prepareForLikelihood(phycas.tree)

# start MCMC
#raw_input('about to run')
phycas.run()

