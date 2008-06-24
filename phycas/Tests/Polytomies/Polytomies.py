# This example program simulates data on a 5-taxon tree containing one polytomy
# This data set is then analyzed using an MCMC analysis that allows polytomies
# to demonstrate that the polytomous true tree can be recovered

from phycas import *

phycas = Phycas()

print
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Simulating a data set'
print '~~~~~~~~~~~~~~~~~~~~~'

# Define the names of the taxa to use when the simulated data set is saved to a file
phycas.taxon_labels = ['P._fimbriata', 'P._parksii', 'P._articulata', 'P._gracilis', 'P._macrophylla']

# Create a model tree containing one polytomy
phycas.ntax = 5
phycas.tree_topology = Newick('(1:0.1,2:0.1,(3:0.1,4:0.1,5:0.1):0.1)')
phycas.tree = phycas.tree_topology.buildTree()

# Create a model
phycas.default_model = 'hky'
phycas.model = Likelihood.HKYModel()
phycas.model.setKappa(4.0)
phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))
phycas.model.setStateFreqParamPrior(ProbDist.ExponentialDist(1.0))

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

# Save simulated data to a NEXUS file using phycas.taxon_labels, datatype=dna and
# using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
sim_data.saveToNexusFile('simHKY.nex', phycas.taxon_labels, 'dna', ('a','c','g','t'))

# Add a MrBayes block to make it easier to summarize trees later
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

# Tell Phycas that the data is already in memory (it will not obtain data
# from a file)
phycas.data_source = 'memory'

print
print '~~~~~~~~~~~~~~~~~~~~~~'
print 'HKY analysis beginning'
print '~~~~~~~~~~~~~~~~~~~~~~'

# Tell phycas we want to allow polytomies
phycas.allow_polytomies = True
phycas.polytomy_prior   = True
phycas.topo_prior_C     = 2.0

# Start with a random tree
phycas.starting_tree_source = 'random'
phycas.outfile_prefix = 'HKYpolytomy'

# By default, Phycas assumes slice_max_units = 1000, but here we will let the maximum
# number of slice sampler units be the largest possible unsigned int. To do this, set
# slice_max_units to 0
phycas.slice_max_units = 0

phycas.ncycles = 200
phycas.sample_every = 10
phycas.adapt_first = 100

phycas.mcmc()

