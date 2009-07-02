# This file does a short MCMC analysis under the GTR+I+G model

from phycas import *
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

rng = ProbDist.Lot()
rng.setSeed(13579)

model.type = 'codon'
model.num_rates = 1
#model.gamma_shape = 2.0
model.pinvar_model = False
#model.use_inverse_shape = True
model.use_flex_model = False
model.edgelen_prior = Exponential(1.0)
model.update_freqs_separately = False
model.state_freqs = [1.0]*61
model.state_freq_prior = Dirichlet([1.0]*61)
model.state_freq_param_prior = Exponential(1.0)
#model.update_relrates_separately = True

mcmc.out.log = 'output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params'
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.ncycles = 100
mcmc.sample_every = 5
mcmc.report_every = 20
mcmc.adapt_first = 10
mcmc.verbose = True
mcmc.ls_move_weight = 10
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
mcmc.rng = rng
mcmc.data_source = blob.characters

if False:
	import sys,os
	if os.path.basename(sys.executable) == 'python_d.exe':
		raw_input('debug stop')

mcmc()
