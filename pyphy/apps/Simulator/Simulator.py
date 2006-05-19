import Conversions
from Likelihood import SimData
from Phycas import *

phycas = Phycas()

# Define the names of the taxa to use when the simulated data set is saved to a file
phycas.taxon_labels = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

# Create a model tree
phycas.tree = Phylogeny.Tree()
model_tree = '(0:0.1,1:0.15,(2:0.025,3:0.15):0.05)'
phycas.tree.buildFromString(model_tree)

# Create a model
phycas.model = Likelihood.HKYModel()
phycas.model.setKappa(4.0)
phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)

# Create a likelihood object to orchestrate both simulations and likelihood calculations
phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)

# Prepare the tree for simulation (i.e. equip nodes with transition matrices)
phycas.likelihood.prepareForSimulation(phycas.tree)

# Simulation settings
phycas.r.setSeed(13579)
phycas.sim_nreps = 1 # ignored at present
phycas.sim_outfile = 'simout.nex'
num_sites = 100000

# Create a SimData object to hold the simulated data
sim_data = SimData()

# Simulate num_sites of data and store in sim_data
# Use the function simulateFirst (rather than just simulate) in order
# to force calculation of transition probabilities
phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)

# Save simulated data to a NEXUS file using phycas.taxon_labels, datatype=dna and
# using the symbols a, c, g, and t for state codes 0, 1, 2, and 3, respectively
sim_data.saveToNexusFile('simulated.nex', phycas.taxon_labels, 'dna', ('a','c','g','t'))

# Copy the simulated data from sim_data to phycas.likelihood so that
# we can compute the likelihood for the simulated data
phycas.likelihood.copyDataFromSimData(sim_data)

# Prepare the tree for calculating the likelihood (i.e. equip nodes with
# both transition matrices and conditional likelihood arrays)
phycas.likelihood.prepareForLikelihood(phycas.tree)

# Compute likelihood for simulated data
print 'lnL =',phycas.likelihood.calcLnL(phycas.tree.getFirstPreorder())

# Add a PAUP block to the simulated.nex file to make it easy to check the results
f = file('simulated.nex', 'a')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  set criterion=likelihood storebrlen;')
f.write('\nend;')
f.write('\n')
f.write('\nbegin trees;')
f.write('\n  translate')
for i,nm in enumerate(phycas.taxon_labels):
    if nm.count(' ') > 0:
        f.write("\n    %d '%s'" % (i, nm))
    else:
        f.write("\n    %d %s" % (i, nm))
    if i < len(phycas.taxon_labels) - 1:
        f.write(',')
    else:
        f.write(';')
f.write('\n  utree simtree = %s;' % model_tree)
f.write('\nend;')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=equal;')
f.write('\n  lscores 1 / userbrlen;')
f.write('\n  lset basefreq=estimate tratio=estimate;')
f.write('\n  lscores 1 / nouserbrlen;')
f.write('\nend;')
f.write('\n')
f.close()
