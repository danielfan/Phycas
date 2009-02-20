from phycas import *
from phycas.ProbDist import Exponential,Gamma
import shutil

burnin     = 1000
ncycles    = 1000
samplefreq = 10
lsweight   = 100
rnseed     = 28767

rng = ProbDist.Lot()
rng.setSeed(rnseed)

model.type                 = 'hky'
model.num_rates            = 4
model.pinvar_model         = False
model.edgelen_prior        = Exponential(10.0)
model.edgelen_hyperprior   = None 
model.kappa_prior          = Exponential(0.2)
model.state_freqs           = [0.25, 0.25, 0.25, 0.25]
model.fix_freqs            = True
model.use_flex_model       = False

randomtree.rng             = rng

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
mcmc.out.log               = 'ps_k80.log'
mcmc.out.trees             = 'ps_k80.t'
mcmc.out.params            = 'ps_k80.p'

mcmc.draw_directly_from_prior = False
mcmc.ndecimals             = 8

ps.maxbeta                 = 1.0
ps.minbeta                 = 0.0
ps.nbetavals               = 101
ps()

