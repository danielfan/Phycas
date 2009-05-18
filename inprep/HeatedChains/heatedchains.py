from phycas import *

#filename = getPhycasTestData('green.nex')
filename = 'green.nex'
blob = readFile(filename)

myrng = ProbDist.Lot()
myrng.setSeed(13579)

model.type = 'jc'
model.num_rates = 1
model.gamma_shape = 2.0
model.pinvar_model = False
model.use_inverse_shape = False
model.use_flex_model = False
model.edgelen_prior = Exponential(1.0)
#model.update_freqs_separately = False
#model.update_relrates_separately = False
#model.edgelen_hyperprior = None

mcmc.out.log = 'output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params'
mcmc.out.params.mode = REPLACE
mcmc.nchains = 8
mcmc.min_heat_power = 1.0/2.6
#mcmc.heat_vector = [1.0]
mcmc.ncycles = 10000
mcmc.sample_every = 5
mcmc.report_every = 20
mcmc.adapt_first = 10
mcmc.verbose = True
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
n = len(blob.taxon_labels)
print 'ntaxa =',n
mcmc.starting_tree_source = randomtree(n_taxa=n, rng = myrng)
mcmc.rng = myrng
mcmc.data_source = blob.characters

if False:
	import sys,os
	if os.path.basename(sys.executable) == 'python_d.exe':
		raw_input('debug stop')

mcmc.debugging = False
mcmc()

#min_heat = 0.6
#n_levels = 5

#def lamda(z, n):
#    return (1-z)/(z*(n-1))

#t = [1/(1+i*lamda(min_heat, n_levels)) for i in range(n_levels)]
#print t
