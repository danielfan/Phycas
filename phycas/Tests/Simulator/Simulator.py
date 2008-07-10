from phycas import *

# Define the names of the taxa to use when the simulated data set is saved to a file
sim.taxon_labels = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

# Create a model tree
sim.tree_source = 'usertree'
sim.tree_topology = Newick('(1:0.1,2:0.15,(3:0.025,4:0.15):0.05)')

# Create a model
sim.default_model = 'hky'
sim.kappa = 4.0
sim.base_freqs = [0.1, 0.2, 0.3, 0.4]
#sim.pinvar_model = False

# Simulation settings
sim.random_seed = 13579
sim.file_name = 'simulated.nex'
sim.nchar = 100000

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('attach debugger to python_d process now')

# Simulate data
sim()

# Now compute the likelihood of the model tree
#phycas.data_source = 'file'
#phycas.data_file_name = 'simulated.nex'
#lnL = phycas.likelihood()
#print 'lnL =',lnL

# Add a PAUP block to the simulated.nex file to make it easy to check the results
f = file('simulated.nex', 'a')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  set criterion=likelihood storebrlen;')
f.write('\nend;')
f.write('\n')
f.write('\nbegin trees;')
f.write('\n  translate')
for i,nm in enumerate(sim.taxon_labels):
    if nm.count(' ') > 0:
        f.write("\n    %d '%s'" % (i + 1, nm))
    else:
        f.write("\n    %d %s" % (i + 1, nm))
    if i < len(sim.taxon_labels) - 1:
        f.write(',')
    else:
        f.write(';')
f.write('\n  utree simtree = %s;' % sim.tree_topology)
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
