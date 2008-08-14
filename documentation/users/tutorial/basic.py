from phycas import *
mcmc.data_source = 'green.nex'
mcmc.out.log = 'basic.log'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'green'
mcmc.out.params.prefix = 'green'
mcmc.ncycles = 2000
mcmc.sample_every = 10
mcmc()
