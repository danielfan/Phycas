# This example program simulates data from a four-taxon tree using the JC model,
# drawing edge lengths in the model tree from an exponential distribution with mean
# 0.02. These data are then analyzed by a model having a different edge length
# prior mean and the Gelfand-Ghosh measure is computed so that the model used 
# for the analysis can be compared to other models.
#
# Warning: this is a more advanced example that involves more than simply changing
# default values and saying go. If this is the first example you have looked at,
# you might want to consider looking at Simplest or Paradox first.

from phycas import *
from math import exp

# Here are some numbers we will use later on in this script. They are defined
# here to make it easy to find later should we wish to run the script again
# with different values.
rnseed = 12345
num_sites = 2000
true_edgelen_mean = 0.02
analysis_prior_mean = 0.00002

# Create a Phycas object (always the first step)
phycas = Phycas()

#############################################################
#################### Simulate a data set ####################
#############################################################

# Create a pseudorandom number generator
rng = ProbDist.Lot()
rng.setSeed(rnseed)

# Set up the substitution model that will be used for both simulation and the MCMC analysis
jcmodel = Likelihood.JCModel() # use the Jukes-Cantor (1969) model
jcmodel.setNGammaRates(1)      # specifying one rate means *no* rate heterogeneity

# Create a tree
true_tree = Phylogeny.Tree()
true_tree.buildFromString('(1,2,(3,4))')

# Apply edge lengths randomly
true_edgelen_dist = ProbDist.ExponentialDist(1.0/true_edgelen_mean)
true_edgelen_dist.setLot(rng)
tm = Phylogeny.TreeManip(true_tree)
tm.setRandomEdgeLengths(true_edgelen_dist)

# Create a likelihood object to orchestrate both simulations and likelihood calculations
phycas.likelihood = Likelihood.TreeLikelihood(jcmodel)

# Prepare the tree for simulation (i.e. equip nodes with transition matrices)
phycas.likelihood.prepareForSimulation(true_tree)

# Create a SimData object to hold the simulated data
sim_data = Likelihood.SimData()

# Simulate num_sites of data and store in sim_data
# Use the function simulateFirst (rather than just simulate) in order
# to force calculation of transition probabilities
phycas.likelihood.simulateFirst(sim_data, true_tree, rng, num_sites)

# Define the names of the taxa to use when the simulated data set is saved to a file
taxon_labels = ['A', 'B', 'C', 'D']

# Save simulated data to a NEXUS file using taxon_labels, datatype=dna and
# using the symbols a, c, g, and t to represent the nucleotides
sim_data.saveToNexusFile('simulated.nex', taxon_labels, 'dna', ('a','c','g','t'))

####################################################################
#################### Analyze the simulated data ####################
####################################################################

# Set various options before starting the MCMC run. To see all available options, you
# can go to the source: take a look at the beginning part (the __init__ function)
# of the Phycas.py file
#   phycas/Examples/GelfandGhosh/GelfandGhosh.py  <-- you are here
#   phycas/Phycas/Phycas.py                       <-- here is Phycas.py

# Use the same random number seed in the analysis that we used for the simulation
phycas.random_seed = rnseed

# Specify the model to be used for the analysis
phycas.default_model = 'jc' # use the Jukes-Cantor (1969) model
phycas.num_rates = 1        # specifying 1 rate means no rate heterogeneity

# The JC model has only edge length parameters, so the only priors we need to worry
# about are the edge length priors. The master_edgelen_dist specifies the prior used
# for all edge lengths. Tell Phycas to not use a hyperprior for edge lengths, otherwise
# it will ignore your master_edgelen_dist specification.
phycas.using_hyperprior = False 
phycas.master_edgelen_dist = ProbDist.ExponentialDist(1.0/analysis_prior_mean)

# Don't allow polytomies
phycas.allow_polytomies = False

# Specify the data file
phycas.data_file_name = 'simulated.nex'

# Start with a random tree
phycas.starting_tree_source = 'random'

# Tell phycas that we want to run the MCMC analysis for 20000 cycles.
phycas.ncycles = 20000
phycas.sample_every = 10                 # save tree and parameters every 10 cycles
phycas.report_every = phycas.ncycles/20  # report progress only 20 times during the run

# Tell Phycas to calculate the Gelfand-Ghosh measure as the MCMC analysis runs
phycas.gg_do = True

# Finally, call the mcmc method, which prepares and then runs the MCMC sampler
phycas.mcmc()
