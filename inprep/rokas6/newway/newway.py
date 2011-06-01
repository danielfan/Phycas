import sys
from phycas import *

setMasterSeed(13759)

model.type = 'jc'
model.num_rates = 1
model.pinvar_model = False

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

mcmc.data_source = 'rokas6first1000.nex'
mcmc.out.log = '_output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = '_trees'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = '_params'
mcmc.out.params.mode = REPLACE
mcmc.ncycles = 10000
mcmc.sample_every = 100
mcmc.report_every = 1000

mcmc.fix_topology = False
mcmc.allow_polytomies = False

mcmc.starting_tree_source = TreeCollection(newick='(1:0.06304,2:0.03328,(3:0.09271,(4:0.05231,(6:0.28321,5:0.09495):0.03792):0.04309):0.03581)')

ss.override_fixed_topology_restriction = True
ss.refdist_definition_file = 'sumt_ref_dist.txt'
ss.nbetavals = 21
ss.xcycles = -(mcmc.ncycles - 1)
ss()

sump.out.log = '_ss.sump.txt'
sump.out.log.mode = REPLACE
sump.file = '_params.p'
sump()
