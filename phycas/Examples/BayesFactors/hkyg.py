from phycas import *
from ProbDist import Exponential
import shutil

data       = getPhycasTestData('green.nex')
ncycles    = 10000
samplefreq = 100
lsweight   = 100
rnseed     = 75391

rng = ProbDist.Lot()
rng.setSeed(rnseed)

model.type                 = 'hky'
model.num_rates            = 4
model.pinvar_model         = False
model.edgelen_prior        = Exponential(10.0)
model.edgelen_hyperprior   = None 
model.use_flex_model = False

mcmc.data_source           = data.characters
mcmc.starting_tree_source  = randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
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
mcmc.random_seed           = rnseed
mcmc.out                   = 'hkyg.log'

ps.maxbeta                 = 1.0
ps.minbeta                 = 0.0
ps.nbetavals               = 101
ps()

