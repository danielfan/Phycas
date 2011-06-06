import sys
from phycas import *

setMasterSeed(25215436)

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

step_tag = sys.argv[1]
mcmc.data_source = '../../rokas8T1000C.nex'
mcmc.out.log = 'output.multitree.' + step_tag + 'steps.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees.multitree' + step_tag + 'steps'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params.multitree' + step_tag + 'steps'
mcmc.out.params.mode = REPLACE
mcmc.ncycles = int(step_tag)
mcmc.sample_every = max(1, mcmc.ncycles//1000)
print 'running', mcmc.ncycles, 'cycles.  Sampling every', mcmc.sample_every
mcmc.report_every = 100

model.type = 'gtr'
model.state_freqs = [0.25, 0.25, 0.25, 0.25]
model.state_freq_prior = Dirichlet((1.0,1.0,1.0,1.0))
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.relrates = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)
model.pinvar_model = False
model.pinvar = 0.5
model.pinvar_prior = Beta(1.0, 1.0)

mcmc.fix_topology = False
mcmc.allow_polytomies = False

newick = '(1:0.07057,2:0.03471,(((5:0.13251,4:0.05314):0.03344,((8:1.14146,7:0.40024):0.15627,6:0.36899):0.20003):0.04305,3:0.10839):0.03644)'
mcmc.starting_tree_source = TreeCollection(newick=newick)

ss.override_fixed_topology_restriction = True
ss.refdist_definition_file = 'refdist.txt'
ss.nbetavals = 25
ss.xcycles = -(mcmc.ncycles - 1)
ss()

sump.out.log = 'ss.sump.multitree' + step_tag + 'steps.txt'
sump.out.log.mode = REPLACE
sump.file = 'params.multitree' + step_tag + 'steps.p'
sump()
