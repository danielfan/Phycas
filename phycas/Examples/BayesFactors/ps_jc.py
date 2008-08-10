from phycas import *
from phycas.ProbDist import Exponential
import shutil

burnin     = 1000
ncycles    = 1000
samplefreq = 10
lsweight   = 100
rnseed     = 11193

rng = ProbDist.Lot()
rng.setSeed(rnseed)

model.type                 = 'jc'
model.num_rates            = 4
model.pinvar_model         = False
model.edgelen_prior        = Exponential(10.0)
model.edgelen_hyperprior   = None 
model.use_flex_model       = False

mcmc.data_source           = 'green.nex'
mcmc.adapt_first           = 10
mcmc.verbose               = True
mcmc.tree_scaler_weight    = 1
mcmc.slice_weight          = 1
mcmc.burnin                = burnin
mcmc.ncycles               = ncycles
mcmc.sample_every          = samplefreq
mcmc.report_every          = ncycles//20
mcmc.ls_move_weight        = lsweight
mcmc.allow_polytomies      = False
mcmc.rng                   = rng
mcmc.out.log               = 'ps_jc.log'
mcmc.out.trees             = 'ps_jc.t'
mcmc.out.params            = 'ps_jc.p'

mcmc.draw_directly_from_prior = False

ps.maxbeta                 = 1.0
ps.minbeta                 = 0.0
ps.nbetavals               = 101
ps()

