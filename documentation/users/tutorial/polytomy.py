from phycas import *

mcmc.allow_polytomies = True
mcmc.polytomy_prior = False
mcmc.topo_prior_C = 1.0

mcmc.data_source = None
mcmc.ntax = 5

mcmc.ncycles = 20000
mcmc.sample_every = 1

mcmc()

sumt.trees = 'trees.t'
sumt.burnin = 1
sumt.tree_credible_prob = 1.0
sumt()


