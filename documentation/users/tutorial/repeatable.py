from phycas import *

r = ProbDist.Lot()
r.setSeed(98765)
randomtree.rng = r
mcmc.rng = r

mcmc.data_file_name = 'green.nex'
mcmc.out.log = 'basic.log'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'green'
mcmc.out.params.prefix = 'green'
mcmc.ncycles = 2000
mcmc.sample_every = 10
mcmc()
