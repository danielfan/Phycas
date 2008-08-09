from phycas import *
from phycas.ProbDist import Exponential, Beta
import shutil

ncycles    = 10000
samplefreq = 100
lsweight   = 100
rnseed     = 75391

rng = ProbDist.Lot()
rng.setSeed(rnseed)

model.type                  = 'hky'
model.num_rates             = 4
model.pinvar_model          = True
model.edgelen_prior         = Exponential(10.0)
model.edgelen_hyperprior    = None 
model.kappa_prior           = Exponential(0.2)
model.base_freq_param_prior = Exponential(1.0)
model.gamma_shape_prior     = Exponential(1.0)
model.pinvar_prior          = Beta(1.0, 1.0)
model.use_flex_model        = False

mcmc.data_source           = 'green6.nex'
#mcmc.starting_tree_source  = randomtree(n_taxa=6, rng=rng)
mcmc.adapt_first           = 10
mcmc.verbose               = True
mcmc.tree_scaler_weight    = 1
mcmc.slice_weight          = 1
mcmc.burnin                = ncycles//10
mcmc.ncycles               = 10000
mcmc.sample_every          = samplefreq
mcmc.report_every          = ncycles//20
mcmc.ls_move_weight        = lsweight
mcmc.allow_polytomies      = False
mcmc.rng                   = rng

ps.maxbeta                 = 1.0
ps.minbeta                 = 0.0
ps.nbetavals               = 101
ps()

