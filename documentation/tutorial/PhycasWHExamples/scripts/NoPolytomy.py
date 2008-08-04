#!/usr/bin/env python

from phycas import *

file_contents = readFile("ShoupLewis.nex")
mcmc.data_source = file_contents.characters
model.help()
model.type = "hky"
model.num_rates = 4

model.current()

model.gamma_shape_prior = Exponential(2.0)
model.gamma_shape_prior.getMean()
model.kappa_prior = Exponential(0.25)
model.base_freq_param_prior = Gamma(1.0, 1.0)
model.edgelen_hyperprior = InverseGamma(2.10000, 0.90909)

mcmc.allow_polytomies = False

rng = Lot()
rng.setSeed(13957)
mcmc.rng = rng
mcmc.starting_tree_source = randomtree(rng=rng)

mcmc.ncycles = 2000
mcmc.sample_every = 10

mcmc.out.trees.prefix = "no_p_trees" 
mcmc.out.trees.mode = REPLACE 
mcmc.out.params.prefix = "no_p_params" 
mcmc.out.params.mode = REPLACE 
mcmc.out.log.prefix = "no_p_output" 
mcmc.out.log.mode = REPLACE 
mcmc.current()

mcmc()

sumt.outgroup_taxon = "Oedogonium cardiacum" 
sumt.trees = "no_p_trees.t" 
sumt.burnin = 101 
sumt.out.log.prefix = "no_p_output" 
sumt.out.log.mode = APPEND 
sumt.out.trees.prefix = "no_p_sumt_trees" 
sumt.out.trees.mode = REPLACE 
sumt.out.splits.prefix = "no_p_sumt_splits" 
sumt.out.splits.mode = REPLACE 

sumt()
