import sys
from phycas import *

setMasterSeed(13579)

model.type = 'jc'
model.num_rates = 1
model.pinvar_model = False

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

mcmc.data_source = 'rokas6first1000.nex'
mcmc.out.log = 'output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params'
mcmc.out.params.mode = REPLACE
mcmc.ncycles = 10000
mcmc.sample_every = 10
mcmc.report_every = 1000

mcmc.fix_topology = False
mcmc.allow_polytomies = False

mcmc.starting_tree_source = TreeCollection(newick='(1:0.064188,2:0.037378,(3:0.088169,(4:0.087427,(5:0.081579,6:0.836677):0.042230):0.042388):0.035303)')
mcmc()

sump.file = 'params.p'
sump.burnin = 1
sump()

sumt.trees = 'trees.t'
sumt.burnin = 1
sumt.out.refdistfile.prefix = 'sumt_ref_dist'
sumt.out.refdistfile.mode = REPLACE
sumt()