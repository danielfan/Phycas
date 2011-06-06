import sys
from phycas import *

setMasterSeed(13759)

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

mcmc.data_source = sys.argv[1]
mcmc.out.log = 'output.pilot.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees.pilot'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params.pilot'
mcmc.out.params.mode = REPLACE
mcmc.ncycles = 50000
mcmc.sample_every = 10
mcmc.report_every = 1000

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

rerunning_for_sump = True
if rerunning_for_sump:
    sump.out.refdistfile = 'sump.refdist.txt'
    sump.out.log = 'sump.pilot.txt'
    sump.out.log.mode = REPLACE
    sump.file = 'params.pilot.p'
    sump()
    sys.exit(0)

mcmc()

sumt.burnin = 5
sumt.out.refdistfile.prefix = "refdist"
sumt.out.refdistfile.mode = REPLACE
sumt.trees = 'trees.pilot.t'
sumt.out.log = 'sumt.pilot.txt'
sumt.out.log.mode = REPLACE
sumt.out.splits.prefix = 'sumt.pilot.splits'
sumt.out.splits.mode = REPLACE
sumt.out.trees.prefix = 'sumt.pilot.trees'
sumt.out.trees.mode = REPLACE
sumt()
