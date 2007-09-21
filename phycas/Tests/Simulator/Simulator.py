from phycas import *

phycas = Phycas()

# Set data_source to None because we will not be reading data from a file
phycas.data_source = None

# Define the names of the taxa to use when the simulated data set is saved to a file
phycas.sim_taxon_labels = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

# Create a model tree
phycas.starting_tree_source = 'usertree'
phycas.tree_topology = '(0:0.1,1:0.15,(2:0.025,3:0.15):0.05)'

# Create a model
phycas.default_model = 'hky'
phycas.starting_kappa = 4.0
phycas.starting_freqs = [0.1, 0.2, 0.3, 0.4]

# Simulation settings
phycas.random_seed = 13579
phycas.sim_file_name = 'simulated.nex'
phycas.sim_nchar = 100000

# Simulate data
phycas.simulateDNA()

# Now compute the likelihood of the model tree
phycas.data_source = 'file'
phycas.data_file_name = 'simulated.nex'
lnL = phycas.likelihood()
print 'lnL =',lnL

# Add a PAUP block to the simulated.nex file to make it easy to check the results
f = file('simulated.nex', 'a')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  set criterion=likelihood storebrlen;')
f.write('\nend;')
f.write('\n')
f.write('\nbegin trees;')
f.write('\n  translate')
for i,nm in enumerate(phycas.sim_taxon_labels):
    if nm.count(' ') > 0:
        f.write("\n    %d '%s'" % (i, nm))
    else:
        f.write("\n    %d %s" % (i, nm))
    if i < len(phycas.sim_taxon_labels) - 1:
        f.write(',')
    else:
        f.write(';')
f.write('\n  utree simtree = %s;' % phycas.tree_topology)
f.write('\nend;')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  lset nst=2 variant=hky basefreq=(0.1 0.2 0.3) tratio=1.8333333 rates=equal;')
f.write('\n  lscores 1 / userbrlen;')
f.write('\n  lset basefreq=estimate tratio=estimate;')
f.write('\n  lscores 1 / nouserbrlen;')
f.write('\n  describe 1 / brlens;')
f.write('\nend;')
f.write('\n')
f.close()
